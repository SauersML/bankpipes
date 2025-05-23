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
import json
import pandas as pd
import requests
import gcsfs
import hail as hl

# --- Globals ---
_FS = None 
_FS_PROJECT_ID = None # Project ID used to initialize _FS
_GCS_INIT_ATTEMPTS = 3
_HAIL_INIT_ATTEMPTS = 3
_LOCAL_CACHE_DIR_NAME = "local_script_cache"

# --- GCS Interaction ---

def get_gcs_fs(project_id_for_billing=None):
    global _FS, _FS_PROJECT_ID

    # Re-initialize if FS is not set, or if a project ID is provided and it's different from the cached one
    if _FS is None or (project_id_for_billing is not None and _FS_PROJECT_ID != project_id_for_billing):
        if project_id_for_billing:
            print(f"Initializing GCSFileSystem (requester_pays=True, project={project_id_for_billing})...")
        else:
            print("Initializing GCSFileSystem (requester_pays=True, project not specified, relying on ADC default for billing)...")
        
        current_attempt = 0
        while current_attempt < _GCS_INIT_ATTEMPTS:
            current_attempt += 1
            try:
                fs_instance = gcsfs.GCSFileSystem(project=project_id_for_billing, requester_pays=True)
                # This helps catch auth/permission issues early.
                if project_id_for_billing:
                    pass
                print("GCSFileSystem initialized successfully.")
                _FS = fs_instance
                _FS_PROJECT_ID = project_id_for_billing # Cache the project ID used
                break 
            except Exception as e:
                print(f"Error initializing GCSFileSystem (Attempt {current_attempt}/{_GCS_INIT_ATTEMPTS}, project: {project_id_for_billing}): {e}")
                if current_attempt < _GCS_INIT_ATTEMPTS:
                    print("Retrying GCS initialization...")
                    time.sleep(5 * current_attempt)
                else:
                    print("FATAL: Could not initialize GCSFileSystem after multiple attempts.")
                    sys.exit(1)
    elif _FS is not None and project_id_for_billing is None and _FS_PROJECT_ID is not None:
        print(f"Returning GCSFileSystem previously initialized with project: {_FS_PROJECT_ID}")
        pass


    return _FS

def gcs_path_exists(gcs_path, project_id_for_billing=None):
    fs = get_gcs_fs(project_id_for_billing=project_id_for_billing)
    if not gcs_path:
        return False
    try:
        return fs.exists(gcs_path)
    except Exception as e:
        print(f"[GCS CHECK WARNING] Failed to check existence of {gcs_path}: {e}")
        return False

def delete_gcs_path(gcs_path, project_id_for_billing=None, recursive=True):
    fs = get_gcs_fs(project_id_for_billing=project_id_for_billing)
    if not gcs_path_exists(gcs_path, project_id_for_billing=project_id_for_billing):
        print(f"GCS path {gcs_path} does not exist. No deletion needed.")
        return True

    print(f"Attempting to delete GCS path: {gcs_path}")
    try:
        delete_attempts = 3
        for i in range(delete_attempts):
            try:
                fs.rm(gcs_path, recursive=recursive)
                time.sleep(2 + i * 2) 
                if not fs.exists(gcs_path): 
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
        return False 
    except Exception as delete_err: 
        print(f"[ERROR] Overall exception during gcsfs deletion process for {gcs_path}: {delete_err}")
        return False

# --- Hail Interaction ---

def init_hail(gcs_hail_temp_dir, log_suffix="task", spark_configurations_json_str=None):
    print("Configuring Hail environment...")
    if not gcs_hail_temp_dir or not gcs_hail_temp_dir.startswith("gs://"):
        print(f"FATAL ERROR: Invalid GCS path for Hail temp directory: {gcs_hail_temp_dir}")
        sys.exit(1)

    spark_conf_dict = {}
    if spark_configurations_json_str:
        try:
            config_list_of_strings = json.loads(spark_configurations_json_str)
            for conf_item_str in config_list_of_strings:
                if '=' in conf_item_str:
                    key, value = conf_item_str.split('=', 1)
                    spark_conf_dict[key.strip()] = value.strip()
                else:
                    print(f"[WARN] Hail Init: Skipping malformed Spark configuration string: {conf_item_str}")
            if not spark_conf_dict:
                 print("[INFO] Hail Init: No valid Spark configurations found in provided JSON string.")
        except json.JSONDecodeError as e:
            print(f"[WARN] Hail Init: Could not parse spark_configurations_json_str: {e}. Proceeding without these extra Spark confs.")
            spark_conf_dict = {} # Default to empty if parsing fails
    
    if not spark_conf_dict:
        print("[INFO] Hail Init: No explicit Spark configurations provided via JSON. Hail will use defaults and environment settings.")

    for attempt in range(_HAIL_INIT_ATTEMPTS):
        try:
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            log_file_name = f'hail_{timestamp}_{log_suffix}_{os.getpid()}.log'
            print(f"Attempting Hail initialization (Attempt {attempt + 1}/{_HAIL_INIT_ATTEMPTS}). Log: ./{log_file_name}")
            print(f"Passing Spark configurations to hl.init: {spark_conf_dict if spark_conf_dict else 'None'}")

            hl.init(
                tmp_dir=gcs_hail_temp_dir,
                log=log_file_name,
                default_reference='GRCh38',
                spark_conf=spark_conf_dict if spark_conf_dict else None,
                master='local[*]'
            )
            print(f"Hail initialized successfully. Log file: ./{log_file_name}")
            print(f"Spark context master: {hl.spark_context().master if hl.spark_context() else 'Not available'}")
            if hl.spark_context():
                 print(f"Spark default parallelism: {hl.spark_context().defaultParallelism}")
            return
        except Exception as e:
            print(f"Hail initialization failed (Attempt {attempt + 1}/{_HAIL_INIT_ATTEMPTS}): {e}")
            if attempt < _HAIL_INIT_ATTEMPTS - 1:
                print("Retrying Hail initialization...")
                if hl.utils.java.Env.spark_context():
                    hl.stop()
                time.sleep(10 * (attempt + 1))
            else:
                print("FATAL ERROR: Hail initialization failed after multiple attempts.")
                sys.exit(1)

def hail_path_exists(hail_path, project_id_for_billing=None):
    if not hail_path:
        return False
    try:
        if hl.utils.java.Env.backend() is None:
            print("[HAIL CHECK WARNING] Hail not initialized. Cannot use hl.hadoop_exists for this check.")
            if hail_path.startswith("gs://"):
                print(f"[HAIL CHECK INFO] Falling back to gcsfs check for {hail_path}")
                return gcs_path_exists(hail_path, project_id_for_billing=project_id_for_billing)
            return False
        return hl.hadoop_exists(hail_path)
    except Exception as e:
        print(f"[HAIL CHECK WARNING] Failed to check existence of {hail_path} using hl.hadoop_exists: {e}")
        if hail_path.startswith("gs://"):
            print(f"[HAIL CHECK INFO] Falling back to gcsfs check for {hail_path}")
            return gcs_path_exists(hail_path, project_id_for_billing=project_id_for_billing)
        return False

# --- Local DataFrame Caching ---

def get_cache_dir(cache_dir_name=_LOCAL_CACHE_DIR_NAME):
    cache_path = os.path.join(os.getcwd(), cache_dir_name)
    os.makedirs(cache_path, exist_ok=True)
    return cache_path

def _local_cache_path(name_prefix, prs_id=None, cache_dir_name=_LOCAL_CACHE_DIR_NAME):
    cache_dir = get_cache_dir(cache_dir_name)
    filename_parts = [name_prefix]
    if prs_id:
        filename_parts.append(str(prs_id))
    filename = "_".join(filename_parts) + ".pkl"
    return os.path.join(cache_dir, filename)

def _load_from_local_cache(name_prefix, prs_id=None):
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
    def decorator(func):
        def wrapper(*args, **kwargs):
            prs_id_arg = kwargs.get('prs_id')
            if not prs_id_arg and args and isinstance(args[0], str) and args[0].startswith("PGS"):
                prs_id_arg = args[0]
            
            dynamic_name_with_id = f"{name_prefix}_{prs_id_arg}" if prs_id_arg else name_prefix

            cached_df = _load_from_local_cache(name_prefix, prs_id_arg)
            if cached_df is not None:
                print(f"[LOCAL CACHE HIT] Loading '{dynamic_name_with_id}' from local cache.")
                return cached_df

            print(f"[LOCAL CACHE MISS] Running step '{dynamic_name_with_id}' for function '{func.__name__}'...")
            result_df = func(*args, **kwargs)

            if isinstance(result_df, pd.DataFrame):
                if not result_df.empty:
                    try:
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
    log_prefix = f"[{log_id}] " if log_id else ""
    
    abs_local_gz_path = os.path.join(os.getcwd(), local_gz_path)
    abs_local_txt_path = os.path.join(os.getcwd(), local_txt_path)
    
    os.makedirs(os.path.dirname(abs_local_gz_path), exist_ok=True)
    os.makedirs(os.path.dirname(abs_local_txt_path), exist_ok=True)

    if os.path.exists(abs_local_txt_path):
        print(f"{log_prefix}Using existing decompressed file: {abs_local_txt_path}")
        return abs_local_txt_path

    # If input URL is GCS, copy directly using gcsfs
    if url.startswith("gs://"):
        print(f"{log_prefix}Copying GCS file from: {url}")
        try:
            fs = get_gcs_fs()

            # Determine if download target should be gz or txt based on GCS URI
            target_local_path_for_gcs_copy = abs_local_gz_path if url.endswith(".gz") else abs_local_txt_path
            
            fs.get(url, target_local_path_for_gcs_copy)
            print(f"{log_prefix}Copied GCS file to: {target_local_path_for_gcs_copy}")

            if url.endswith(".gz"): # Decompress if it was a .gz file from GCS
                print(f"{log_prefix}Decompressing {target_local_path_for_gcs_copy} to {abs_local_txt_path}...")
                with gzip.open(target_local_path_for_gcs_copy, 'rb') as f_in, open(abs_local_txt_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                if target_local_path_for_gcs_copy != abs_local_txt_path: # if gz path was used for download
                    os.remove(target_local_path_for_gcs_copy)
                print(f"{log_prefix}Decompression complete. Final file: {abs_local_txt_path}")
                return abs_local_txt_path
            else: # Was not a .gz file from GCS
                return target_local_path_for_gcs_copy # which is abs_local_txt_path in this case

        except Exception as e:
            print(f"{log_prefix}ERROR: Failed to copy/decompress GCS file {url}. Error: {e}")
            return None # Fall through to general error handling
            
    # Original HTTP/FTP download logic
    print(f"{log_prefix}Downloading file from: {url} to {abs_local_gz_path}")
    try:
        response = requests.get(url, stream=True, timeout=600) 
        response.raise_for_status()
        with open(abs_local_gz_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=1024*1024): 
                f.write(chunk)
        print(f"{log_prefix}Download complete: {abs_local_gz_path}")

        final_path_to_return = abs_local_txt_path
        if url.endswith('.gz') or abs_local_gz_path.endswith('.gz'): # Check original URL and local path
            print(f"{log_prefix}Decompressing {abs_local_gz_path} to {abs_local_txt_path}...")
            with gzip.open(abs_local_gz_path, 'rb') as f_in, open(abs_local_txt_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            if abs_local_gz_path != abs_local_txt_path: # If .gz was downloaded to a different name than the final .txt
                os.remove(abs_local_gz_path)
            print(f"{log_prefix}Decompression complete. Final file: {abs_local_txt_path}")
        else:
            if abs_local_gz_path != abs_local_txt_path:
                shutil.move(abs_local_gz_path, abs_local_txt_path)
            print(f"{log_prefix}File is not gzipped. Final file: {abs_local_txt_path}")
        
        return final_path_to_return

    except requests.exceptions.Timeout:
        print(f"{log_prefix}ERROR: Timeout occurred while downloading {url}")
    except requests.exceptions.RequestException as e:
        print(f"{log_prefix}ERROR: Failed to download {url}. Error: {e}")
    except Exception as e:
        print(f"{log_prefix}ERROR: Failed during download/decompression for {url}. Error: {e}")
    finally:
        # Cleanup attempt for gz file if txt file wasn't created or if it was a non-gz file that failed
        if os.path.exists(abs_local_gz_path) and not os.path.exists(abs_local_txt_path):
            try:
                os.remove(abs_local_gz_path)
            except OSError:
                pass
    return None
