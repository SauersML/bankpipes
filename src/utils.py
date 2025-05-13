"""
Utility functions for the PRS Nextflow pipeline.

Includes helpers for:
- GCS File System interaction (requester pays)
- Hail initialization
- GCS path checking and deletion
- File downloading and decompression
- Local disk caching for Pandas DataFrames
"""

import datetime
import os
import sys
import gzip
import shutil
import time
import pandas as pd
# import numpy as np
import requests
import gcsfs
import hail as hl

# --- Globals ---
_FS = None # Global GCSFileSystem object, prefixed to avoid accidental direct use
_GCS_INIT_ATTEMPTS = 3 # Retries for initializing GCS connection
_HAIL_INIT_ATTEMPTS = 3 # Retries for initializing Hail
_LOCAL_CACHE_DIR_NAME = "local_script_cache" # Name of the local cache directory

# --- GCS Interaction ---

def get_gcs_fs():
    """
    Initializes and returns a GCSFileSystem object with requester_pays enabled.
    Uses a global variable _FS to avoid re-initialization. Includes retry logic.
    Exits fatally if initialization fails.
    """
    global _FS
    if _FS is None:
        print("Initializing GCSFileSystem (requester_pays=True)...")
        for attempt in range(_GCS_INIT_ATTEMPTS):
            try:
                fs_instance = gcsfs.GCSFileSystem(requester_pays=True)
                print("GCSFileSystem initialized successfully.")
                _FS = fs_instance
                break # Success
            except Exception as e:
                print(f"Error initializing GCSFileSystem (Attempt {attempt + 1}/{_GCS_INIT_ATTEMPTS}): {e}")
                if attempt < _GCS_INIT_ATTEMPTS - 1:
                    print("Retrying GCS initialization...")
                    time.sleep(5 * (attempt + 1))
                else:
                    print("FATAL: Could not initialize GCSFileSystem after multiple attempts.")
                    print("Check GCS credentials, permissions, and network connectivity.")
                    sys.exit(1)
    return _FS

def gcs_path_exists(gcs_path):
    """Checks if a path (file or directory) exists on GCS using gcsfs."""
    fs = get_gcs_fs()
    if not gcs_path:
        return False
    try:
        return fs.exists(gcs_path)
    except Exception as e:
        print(f"[GCS CHECK WARNING] Failed to check existence of {gcs_path}: {e}")
        return False

def delete_gcs_path(gcs_path, recursive=True):
    """
    Deletes a GCS path using gcsfs with retries.
    Returns True if deletion is successful or path doesn't exist, False otherwise.
    """
    fs = get_gcs_fs()
    if not gcs_path_exists(gcs_path): # Use our own checker first
        print(f"GCS path {gcs_path} does not exist. No deletion needed.")
        return True

    print(f"Attempting to delete GCS path: {gcs_path}")
    try:
        delete_attempts = 3
        for i in range(delete_attempts):
            try:
                fs.rm(gcs_path, recursive=recursive)
                time.sleep(2 + i * 2) # Allow for GCS consistency
                if not fs.exists(gcs_path): # Verify deletion
                    print(f"Deletion confirmed for GCS path: {gcs_path}")
                    return True
                else:
                    print(f"[WARNING] GCS path {gcs_path} still exists after delete attempt {i+1}.")
                    if i == delete_attempts - 1:
                        print(f"[ERROR] Failed to confirm deletion of {gcs_path} after {delete_attempts} attempts.")
                        return False
            except Exception as attempt_err:
                print(f"[ERROR] Attempt {i+1}/{delete_attempts} to delete {gcs_path} failed: {attempt_err}")
                if "rate" in str(attempt_err).lower() and i < delete_attempts - 1:
                    wait_time = 5 + i*5
                    print(f"Rate limit likely hit, retrying delete in {wait_time} seconds...")
                    time.sleep(wait_time)
                elif i == delete_attempts - 1:
                    print(f"[ERROR] All attempts to delete {gcs_path} failed. Last error: {attempt_err}")
                    return False
        return False # Should be unreachable if logic is correct
    except Exception as delete_err:
        print(f"[ERROR] Overall exception during gcsfs deletion process for {gcs_path}: {delete_err}")
        return False

# --- Hail Interaction ---

def init_hail(gcs_hail_temp_dir, log_suffix="task"):
    """
    Initializes Hail with retry logic. Uses a unique log file name.
    Exits fatally if initialization fails.
    """
    print("Configuring Hail environment...")
    # Ensure Hail temp dir is valid gs:// path
    if not gcs_hail_temp_dir or not gcs_hail_temp_dir.startswith("gs://"):
        print(f"FATAL ERROR: Invalid GCS path for Hail temp directory: {gcs_hail_temp_dir}")
        sys.exit(1)

    for attempt in range(_HAIL_INIT_ATTEMPTS):
        try:
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            # Ensure log file is written to current working directory, which Nextflow manages
            log_file_name = f'hail_{timestamp}_{log_suffix}_{os.getpid()}.log'
            print(f"Attempting Hail initialization (Attempt {attempt + 1}/{_HAIL_INIT_ATTEMPTS}). Log: ./{log_file_name}")

            hl.init(tmp_dir=gcs_hail_temp_dir, log=log_file_name, default_reference='GRCh38')
            print(f"Hail initialized successfully. Log file: ./{log_file_name}")
            return
        except Exception as e:
            print(f"Hail initialization failed (Attempt {attempt + 1}/{_HAIL_INIT_ATTEMPTS}): {e}")
            if attempt < _HAIL_INIT_ATTEMPTS - 1:
                print("Retrying Hail initialization...")
                # Stop any existing SparkContext to allow re-initialization
                if hl.utils.java.Env.spark_context():
                    hl.stop()
                time.sleep(10 * (attempt + 1))
            else:
                print("FATAL ERROR: Hail initialization failed after multiple attempts.")
                sys.exit(1)

def hail_path_exists(hail_path):
    """
    Checks if a Hail table/matrix table path directory exists using hl.hadoop_exists.
    Falls back to gcs_path_exists for gs:// paths if Hail check fails unexpectedly.
    """
    if not hail_path:
        return False
    try:
        if hl.utils.java.Env.backend() is None:
            print("[HAIL CHECK WARNING] Hail not initialized. Attempting fallback for GCS path check.")
            if hail_path.startswith("gs://"):
                return gcs_path_exists(hail_path)
            return False # Cannot check non-GCS path without Hail
        return hl.hadoop_exists(hail_path)
    except Exception as e:
        print(f"[HAIL CHECK WARNING] Failed to check existence of {hail_path} using hl.hadoop_exists: {e}")
        if hail_path.startswith("gs://"):
            print(f"[HAIL CHECK INFO] Falling back to gcsfs check for {hail_path}")
            return gcs_path_exists(hail_path)
        return False

# --- Local DataFrame Caching (Mirroring original notebook) ---

def get_cache_dir(cache_dir_name=_LOCAL_CACHE_DIR_NAME):
    """Creates and returns the path to the local cache directory."""
    # Use current working directory for Nextflow tasks
    # Nextflow manages task-specific work directories.
    cache_path = os.path.join(os.getcwd(), cache_dir_name)
    os.makedirs(cache_path, exist_ok=True)
    return cache_path

def _local_cache_path(name_prefix, prs_id=None, cache_dir_name=_LOCAL_CACHE_DIR_NAME):
    """Generates local cache file path, incorporating prs_id if provided."""
    cache_dir = get_cache_dir(cache_dir_name)
    filename_parts = [name_prefix]
    if prs_id:
        filename_parts.append(str(prs_id))
    filename = "_".join(filename_parts) + ".pkl"
    return os.path.join(cache_dir, filename)

def _load_from_local_cache(name_prefix, prs_id=None):
    """Loads pandas DataFrame from local pickle cache."""
    path = _local_cache_path(name_prefix, prs_id)
    if os.path.exists(path):
        try:
            if os.path.getsize(path) > 0:
                return pd.read_pickle(path)
            else:
                print(f"[LOCAL CACHE WARNING] Cache file '{path}' is empty. Will recompute.")
                os.remove(path)
                return None
        except Exception as e:
            print(f"[LOCAL CACHE WARNING] Failed to load '{name_prefix}' (ID: {prs_id}) from cache '{path}': {e}. Will recompute.")
            try:
                os.remove(path)
            except OSError:
                pass
            return None
    return None

def cache_result(name_prefix: str):
    """
    Decorator for caching pandas DataFrame results to local pickle files.
    Mirrors the caching behavior of the original notebook.
    The decorated function can optionally take a 'prs_id' kwarg (or as first arg if string)
    to make the cache file name specific.
    """
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Determine if prs_id is passed for cache naming specificity
            prs_id_arg = kwargs.get('prs_id')
            if not prs_id_arg and args and isinstance(args[0], str) and args[0].startswith("PGS"):
                prs_id_arg = args[0]
            dynamic_name = name_prefix
            if prs_id_arg:
                dynamic_name_with_id = f"{name_prefix}_{prs_id_arg}" # Used for messages and actual load/save
            else:
                dynamic_name_with_id = name_prefix # For messages if no ID, and actual load/save

            cached_df = _load_from_local_cache(name_prefix, prs_id_arg) # Pass potentially None prs_id_arg
            if cached_df is not None:
                print(f"[LOCAL CACHE HIT] Loading '{dynamic_name_with_id}' from local cache.")
                return cached_df

            print(f"[LOCAL CACHE MISS] Running step '{dynamic_name_with_id}' for function '{func.__name__}'...")
            result_df = func(*args, **kwargs)

            if isinstance(result_df, pd.DataFrame):
                if not result_df.empty:
                    try:
                        # Use the same logic for path generation as in _load_from_local_cache
                        save_path = _local_cache_path(name_prefix, prs_id_arg)
                        result_df.to_pickle(save_path)
                        print(f"[LOCAL CACHE SAVE] Saved '{dynamic_name_with_id}' to local cache: {save_path}")
                    except Exception as e:
                        print(f"[LOCAL CACHE WARNING] Failed to save '{dynamic_name_with_id}' to cache: {e}")
                else:
                    print(f"[LOCAL CACHE SKIP] Skipping local caching for empty DataFrame result: '{dynamic_name_with_id}'")
            elif result_df is not None:
                print(f"[LOCAL CACHE SKIP] Result is {type(result_df)}, not DataFrame. Skipping caching for '{dynamic_name_with_id}'.")
            return result_df
        return wrapper
    return decorator


# --- File Downloading ---
def download_file(url: str, local_gz_path: str, local_txt_path: str, log_id: str = "") -> str | None:
    """
    Downloads a file (potentially gzipped) from a URL and decompresses if needed.
    Checks if the final decompressed file already exists locally first.
    Local paths are relative to the current working directory.
    """
    log_prefix = f"[{log_id}] " if log_id else ""

    # Ensure local directories for these paths exist if they are nested
    os.makedirs(os.path.dirname(local_gz_path), exist_ok=True)
    os.makedirs(os.path.dirname(local_txt_path), exist_ok=True)

    if os.path.exists(local_txt_path):
        print(f"{log_prefix}Using existing decompressed file: {local_txt_path}")
        return local_txt_path

    print(f"{log_prefix}Downloading file from: {url} to {local_gz_path}")
    try:
        response = requests.get(url, stream=True, timeout=600) # 10 min timeout
        response.raise_for_status()
        with open(local_gz_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024*1024): # 1MB chunks
                f.write(chunk)
        print(f"{log_prefix}Download complete: {local_gz_path}")

        final_path_to_return = local_txt_path
        if url.endswith('.gz') or local_gz_path.endswith('.gz'):
            print(f"{log_prefix}Decompressing {local_gz_path} to {local_txt_path}...")
            with gzip.open(local_gz_path, 'rb') as f_in, open(local_txt_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(local_gz_path)
            print(f"{log_prefix}Decompression complete. Final file: {local_txt_path}")
        else:
            # If not gzipped, the downloaded file is the final file.
            # If local_gz_path is different from local_txt_path, move it.
            if local_gz_path != local_txt_path:
                shutil.move(local_gz_path, local_txt_path)
            else: # They are the same, gz_path is actually the final txt_path
                 pass
            print(f"{log_prefix}File is not gzipped. Final file: {local_txt_path}")

        return final_path_to_return

    except requests.exceptions.Timeout:
        print(f"{log_prefix}ERROR: Timeout occurred while downloading {url}")
    except requests.exceptions.RequestException as e:
        print(f"{log_prefix}ERROR: Failed to download {url}. Error: {e}")
    except Exception as e:
        print(f"{log_prefix}ERROR: Failed during download/decompression for {url}. Error: {e}")
    finally:
        # Cleanup attempt for gz file if txt file wasn't created or if it was a non-gz file that failed
        if os.path.exists(local_gz_path) and not os.path.exists(local_txt_path):
            try:
                os.remove(local_gz_path)
            except OSError:
                pass
    return None
