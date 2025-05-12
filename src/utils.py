import datetime
import os
import sys
import gzip
import shutil
import time
import pandas as pd
import numpy as np
import requests
import gcsfs
import hail as hl

# --- Global GCS FileSystem Object ---
# Initialized when first needed or by explicit call
FS = None
GCS_INIT_ATTEMPTS = 3

# --- Caching Configuration ---
# This will be relative to the execution directory of the script (Nextflow work_dir for the process)
CACHE_DIR = "local_cache"

def get_cache_dir():
    """Ensures the local cache directory exists."""
    os.makedirs(CACHE_DIR, exist_ok=True)

def get_gcs_fs():
    """Initializes and returns the GCSFileSystem object with requester_pays."""
    global FS
    if FS is None:
        print("Initializing GCSFileSystem...")
        for attempt in range(GCS_INIT_ATTEMPTS):
            try:
                fs_instance = gcsfs.GCSFileSystem(requester_pays=True)
                print("GCSFileSystem initialized successfully (with requester_pays=True).")
                FS = fs_instance
                break
            except Exception as e:
                print(f"Error initializing/testing GCSFileSystem (Attempt {attempt + 1}/{GCS_INIT_ATTEMPTS}): {e}")
                if attempt < GCS_INIT_ATTEMPTS - 1:
                    time.sleep(5)
                else:
                    print("FATAL: Could not initialize GCSFileSystem after multiple attempts.")
                    sys.exit(1)
    return FS

def init_hail(gcs_hail_temp_dir, run_timestamp):
    """Initializes Hail."""
    print("Configuring Hail environment...")
    hail_init_attempts = 3
    for attempt in range(hail_init_attempts):
        try:
            log_file_name = f'hail_prs_log_{run_timestamp}_{os.getpid()}.txt' # unique log name per task
            hl.init(tmp_dir=gcs_hail_temp_dir, log=log_file_name)
            hl.default_reference = 'GRCh38'
            print(f"Hail initialized successfully. Log: {log_file_name}\n")
            # Confirm GCS access through Hail
            try:
                test_gcs_path = f"{gcs_hail_temp_dir}/hail_gcs_access_test_{run_timestamp}_{os.getpid()}.txt"
                hl.utils.hadoop_write(test_gcs_path, "test")
                if hl.hadoop_exists(test_gcs_path):
                    hl.hadoop_remove(test_gcs_path)
                    print("Hail GCS access confirmed via hadoop commands.")
                else:
                    print("WARN: Hail GCS access test file not found after write attempt.")
            except Exception as he:
                print(f"WARN: Hail GCS access test failed: {he}. This might indicate issues with Hail's GCS connector or permissions.")

            break
        except Exception as e:
            print(f"Hail initialization failed (Attempt {attempt + 1}/{hail_init_attempts}): {e}")
            if attempt < hail_init_attempts - 1:
                print("Retrying Hail initialization...")
                time.sleep(10 + attempt * 5) # Increasing sleep
            else:
                print("FATAL ERROR: Hail initialization failed after multiple attempts.")
                sys.exit(1)

def cache_path(name):
    """Generates local cache file path."""
    get_cache_dir()
    return os.path.join(CACHE_DIR, f"{name}.pkl")

def load_from_cache(name):
    """Loads pandas DataFrame from local pickle cache."""
    path = cache_path(name)
    if os.path.exists(path):
        try:
            if os.path.getsize(path) > 0:
                return pd.read_pickle(path)
            else:
                print(f"[CACHE WARNING] Cache file '{path}' is empty. Will recompute.")
                os.remove(path)
                return None
        except Exception as e:
            print(f"[CACHE WARNING] Failed to load '{name}' from cache: {e}. Will recompute.")
            try:
                os.remove(path)
            except OSError:
                pass
            return None
    return None

def cache_result(name_prefix):
    """Decorator for caching pandas DataFrame results to local files."""
    def decorator(func):
        def wrapper(*args, **kwargs):
            dynamic_name_parts = [name_prefix]
            if args:
                for arg in args:
                    if isinstance(arg, str) and arg.startswith("PGS"): # Often prs_id
                        dynamic_name_parts.append(arg)
                    elif isinstance(arg, list) and len(arg) > 0 and isinstance(arg[0], (int,str)): # like concept_ids
                        dynamic_name_parts.append("cids"+str(len(arg)))

            prs_id_kwarg = kwargs.get('prs_id')
            if prs_id_kwarg:
                 dynamic_name_parts.append(str(prs_id_kwarg))

            phenotype_name_kwarg = kwargs.get('phenotype_name')
            if phenotype_name_kwarg:
                dynamic_name_parts.append(str(phenotype_name_kwarg).replace(" ", "_"))

            dynamic_name = "_".join(dynamic_name_parts)

            cached = load_from_cache(dynamic_name)
            if cached is not None:
                print(f"[CACHE HIT] Loading '{dynamic_name}' from local cache: {cache_path(dynamic_name)}")
                return cached

            print(f"[CACHE MISS] Running step to generate '{dynamic_name}'...")
            result = func(*args, **kwargs)

            if isinstance(result, pd.DataFrame):
                if not result.empty:
                    try:
                        pd.to_pickle(result, cache_path(dynamic_name))
                        print(f"[CACHE SAVE] Saved '{dynamic_name}' to local cache: {cache_path(dynamic_name)}")
                    except Exception as e:
                        print(f"[CACHE WARNING] Failed to save '{dynamic_name}' to cache: {e}")
                else:
                    print(f"[CACHE SKIP] Skipping local caching for empty DataFrame result: '{dynamic_name}'")
            elif result is not None:
                print(f"[CACHE SKIP] Skipping local caching for non-DataFrame/non-None result type: {type(result)}")
            return result
        return wrapper
    return decorator

def gcs_path_exists(gcs_path):
    """Checks if a path (file or directory) exists on GCS using gcsfs."""
    fs = get_gcs_fs()
    if fs is None:
        print("[GCS CHECK ERROR] GCSFileSystem not initialized.")
        return False
    try:
        return fs.exists(gcs_path)
    except Exception as e:
        print(f"[GCS CHECK WARNING] Failed to check existence of {gcs_path}: {e}")
        return False

def hail_path_exists(hail_path):
    """Checks if a Hail table/matrix table path _directory_ exists."""
    try:
        return hl.hadoop_exists(hail_path)
    except Exception as e:
        print(f"[HAIL CHECK WARNING] Failed to check existence of {hail_path} using hl.hadoop_exists: {e}")
        if hail_path.startswith("gs://"):
            print(f"[HAIL CHECK INFO] Falling back to gcsfs check for {hail_path}")
            return gcs_path_exists(hail_path)
        return False

def delete_gcs_path(gcs_path, recursive=True):
    """Deletes a GCS path using gcsfs with retries."""
    fs = get_gcs_fs()
    if fs is None:
        print(f"[DELETE ERROR] GCSFileSystem not initialized. Cannot delete {gcs_path}")
        return False
    print(f"Attempting to delete GCS path using gcsfs: {gcs_path}")
    try:
        if fs.exists(gcs_path):
            delete_attempts = 3
            for i in range(delete_attempts):
                try:
                    fs.rm(gcs_path, recursive=recursive)
                    print(f"Successfully initiated delete for GCS path: {gcs_path} (Attempt {i+1}/{delete_attempts})")
                    time.sleep(2 + i * 2) # Increasing delay for GCS to reflect deletion
                    if not fs.exists(gcs_path):
                        print(f"Deletion confirmed for GCS path: {gcs_path}")
                        return True
                    else:
                        print(f"[WARNING] GCS path {gcs_path} still exists after deletion attempt {i+1}.")
                        if i == delete_attempts - 1:
                            print(f"[ERROR] Failed to confirm deletion of {gcs_path} after {delete_attempts} attempts.")
                            return False
                except Exception as attempt_err:
                    print(f"[ERROR] Attempt {i+1}/{delete_attempts} to delete {gcs_path} failed: {attempt_err}")
                    if "rate" in str(attempt_err).lower() and i < delete_attempts - 1:
                        print(f"Rate limit likely hit, retrying delete in {5 + i*5} seconds...")
                        time.sleep(5 + i*5)
                    elif i == delete_attempts - 1:
                        # raise attempt_err # Re-raise the last error if all attempts fail
                        print(f"[ERROR] All attempts to delete {gcs_path} failed. Last error: {attempt_err}")
                        return False # Do not raise, allow script to decide if this is fatal
            return False # Fallback, should be caught by checks above
        else:
            print(f"GCS path {gcs_path} does not exist. No deletion needed.")
            return True
    except Exception as delete_err:
        print(f"[ERROR] Overall exception during gcsfs deletion process for {gcs_path}: {delete_err}")
        return False

def download_file(url, local_filename_gz, local_filename_txt, prs_id=""):
    """Downloads and decompresses a file."""
    # Prefix prs_id for log messages
    log_prefix = f"[{prs_id}] " if prs_id else ""

    os.makedirs(os.path.dirname(local_filename_gz), exist_ok=True)

    if os.path.exists(local_filename_txt):
        print(f"{log_prefix}Using existing decompressed file: {local_filename_txt}")
        return local_filename_txt

    print(f"{log_prefix}Downloading file from: {url} to {local_filename_gz}")
    try:
        response = requests.get(url, stream=True, timeout=300) # 5 min timeout
        response.raise_for_status()
        with open(local_filename_gz, 'wb') as out_file:
            for chunk in response.iter_content(chunk_size=8192):
                out_file.write(chunk)
        print(f"{log_prefix}Downloaded to {local_filename_gz}.")

        print(f"{log_prefix}Decompressing {local_filename_gz}...")
        with gzip.open(local_filename_gz, 'rb') as f_in:
            with open(local_filename_txt, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        os.remove(local_filename_gz)
        print(f"{log_prefix}Decompressed to {local_filename_txt}.\n")
        return local_filename_txt
    except requests.exceptions.RequestException as e:
        print(f"{log_prefix}ERROR: Failed to download {url}. Error: {e}")
        return None
    except Exception as e:
        print(f"{log_prefix}ERROR: Failed to decompress {local_filename_gz}. Error: {e}")
        if os.path.exists(local_filename_gz): os.remove(local_filename_gz)
        if os.path.exists(local_filename_txt): os.remove(local_filename_txt)
        return None
