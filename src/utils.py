"""
Utility functions for the PRS Nextflow pipeline.

Includes helpers for:
- GCS File System interaction (requester pays)
- Hail initialization
- GCS path checking and deletion
- File downloading and decompression
- GCS-based caching for Pandas DataFrames
"""

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

# --- Globals ---
FS = None # Global GCSFileSystem object
GCS_INIT_ATTEMPTS = 3 # Retries for initializing GCS connection
HAIL_INIT_ATTEMPTS = 3 # Retries for initializing Hail

# --- GCS Interaction ---

def get_gcs_fs():
    """
    Initializes and returns a GCSFileSystem object with requester_pays enabled.
    Uses a global variable FS to avoid re-initialization. Includes retry logic.
    Exits fatally if initialization fails.
    """
    global FS
    if FS is None:
        print("Initializing GCSFileSystem (requester_pays=True)...")
        for attempt in range(GCS_INIT_ATTEMPTS):
            try:
                fs_instance = gcsfs.GCSFileSystem(requester_pays=True)
                print("GCSFileSystem initialized successfully.")
                FS = fs_instance
                break # Success
            except Exception as e:
                print(f"Error initializing GCSFileSystem (Attempt {attempt + 1}/{GCS_INIT_ATTEMPTS}): {e}")
                if attempt < GCS_INIT_ATTEMPTS - 1:
                    print("Retrying GCS initialization...")
                    time.sleep(5 * (attempt + 1)) # Exponential backoff might be better
                else:
                    print("FATAL: Could not initialize GCSFileSystem after multiple attempts.")
                    print("Check GCS credentials, permissions, and network connectivity.")
                    sys.exit(1) # Exit if GCS cannot be accessed
    return FS

def gcs_path_exists(gcs_path):
    """Checks if a path (file or directory) exists on GCS using gcsfs."""
    fs = get_gcs_fs()
    if not gcs_path: # Handle empty path strings
        return False
    try:
        return fs.exists(gcs_path)
    except Exception as e:
        print(f"[GCS CHECK WARNING] Failed to check existence of {gcs_path}: {e}")
        return False # Assume it doesn't exist if check fails

def delete_gcs_path(gcs_path, recursive=True):
    """
    Deletes a GCS path using gcsfs with retries.
    Returns True if deletion is successful or path doesn't exist, False otherwise.
    """
    fs = get_gcs_fs()
    print(f"Attempting to delete GCS path: {gcs_path}")
    try:
        if fs.exists(gcs_path):
            delete_attempts = 3
            for i in range(delete_attempts):
                try:
                    fs.rm(gcs_path, recursive=recursive)
                    # Give GCS some time to become consistent
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
                    # Basic retry logic for potential rate limiting
                    if "rate" in str(attempt_err).lower() and i < delete_attempts - 1:
                        wait_time = 5 + i*5
                        print(f"Rate limit likely hit, retrying delete in {wait_time} seconds...")
                        time.sleep(wait_time)
                    elif i == delete_attempts - 1:
                        print(f"[ERROR] All attempts to delete {gcs_path} failed. Last error: {attempt_err}")
                        return False # Do not raise, let caller decide if fatal
            return False # Fallback if loop finishes without success
        else:
            print(f"GCS path {gcs_path} does not exist. No deletion needed.")
            return True
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
    for attempt in range(HAIL_INIT_ATTEMPTS):
        try:
            # Generate a unique log file name using timestamp, pid and suffix
            timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            log_file_name = f'hail_{timestamp}_{log_suffix}_{os.getpid()}.log'
            print(f"Attempting Hail initialization (Attempt {attempt + 1}/{HAIL_INIT_ATTEMPTS}). Log: {log_file_name}")

            # Ensure tmp_dir path is valid
            if not gcs_hail_temp_dir or not gcs_hail_temp_dir.startswith("gs://"):
                 raise ValueError(f"Invalid GCS path provided for Hail temp directory: {gcs_hail_temp_dir}")

            hl.init(tmp_dir=gcs_hail_temp_dir, log=log_file_name, default_reference='GRCh38')

            print(f"Hail initialized successfully. Log file: {log_file_name}")
            return # Success
        except Exception as e:
            print(f"Hail initialization failed (Attempt {attempt + 1}/{HAIL_INIT_ATTEMPTS}): {e}")
            if attempt < HAIL_INIT_ATTEMPTS - 1:
                print("Retrying Hail initialization...")
                time.sleep(10 * (attempt + 1)) # Increasing sleep time
            else:
                print("FATAL ERROR: Hail initialization failed after multiple attempts.")
                sys.exit(1) # Exit if Hail cannot be initialized

def hail_path_exists(hail_path):
    """
    Checks if a Hail table/matrix table path directory exists using hl.hadoop_exists.
    Falls back to gcs_path_exists for gs:// paths if Hail check fails unexpectedly.
    """
    if not hail_path: # Handle empty path strings
        return False
    try:
        # Ensure Hail is initialized before using hadoop commands
        if hl.utils.java.Env.backend() is None:
             print("WARN: Hail not initialized, cannot check path existence via Hail. Attempting gcsfs fallback if applicable.")
             if hail_path.startswith("gs://"): return gcs_path_exists(hail_path)
             else: return False # Cannot check non-GCS paths without Hail backend

        return hl.hadoop_exists(hail_path)
    except Exception as e:
        # Log specific Hail exceptions if possible, e.g., related to auth or backend
        print(f"[HAIL CHECK WARNING] Failed to check existence of {hail_path} using hl.hadoop_exists: {e}")
        # Fallback for GCS paths
        if hail_path.startswith("gs://"):
            print(f"[HAIL CHECK INFO] Falling back to gcsfs check for {hail_path}")
            return gcs_path_exists(hail_path)
        return False # Assume non-existent if check fails for non-GCS paths

# --- GCS DataFrame Caching ---

def cache_result_to_gcs(name_prefix: str, cache_file_type: str = 'csv'):
    """
    Decorator for caching pandas DataFrame results to Google Cloud Storage.

    Checks for a CSV or Parquet file on GCS. If found, loads it.
    Otherwise, executes the decorated function, and saves the resulting DataFrame to GCS.

    Requires the decorated function to accept 'gcs_cache_dir' (str) in its kwargs.
    Optionally accepts 'fs' (gcsfs.GCSFileSystem) in kwargs, otherwise uses global FS.

    Args:
        name_prefix: A prefix for the cache file name.
        cache_file_type: 'csv' (default) or 'parquet'. Parquet is generally preferred
                         for preserving data types and better performance.
    """
    if cache_file_type not in ['csv', 'parquet']:
        raise ValueError("cache_file_type must be 'csv' or 'parquet'")

    def decorator(func):
        def wrapper(*args, **kwargs):
            # --- Get GCS FS and Cache Directory ---
            fs = kwargs.get('fs', get_gcs_fs()) # Allow passing fs or use global
            gcs_cache_dir = kwargs.get('gcs_cache_dir')
            if not gcs_cache_dir or not gcs_cache_dir.startswith("gs://"):
                print(f"ERROR in @cache_result_to_gcs for '{func.__name__}':")
                print(" Decorated function must be called with a valid 'gcs_cache_dir' argument (gs://...).")
                # Decide how to handle: raise error or proceed without caching?
                # Raising error is safer as caching is expected.
                raise ValueError(f"Missing or invalid 'gcs_cache_dir' kwarg for cached function {func.__name__}")

            # --- Generate Dynamic Cache File Name ---
            # Creates a unique name based on prefix and function arguments/kwargs
            # to avoid collisions for different calls.
            dynamic_name_parts = [name_prefix]
            # Incorporate relevant args (avoiding large objects)
            for arg in args:
                 if isinstance(arg, (str, int, float, bool)):
                      dynamic_name_parts.append(str(arg).replace(" ", "_").replace("/", "_"))
                 elif isinstance(arg, list) and len(arg) > 0 and isinstance(arg[0], (int, str)):
                      # Example: handle list of concept IDs
                      dynamic_name_parts.append(f"list{len(arg)}_{hash(tuple(arg))}") # Use hash for lists

            # Incorporate relevant kwargs
            for key, value in kwargs.items():
                # Exclude kwargs used by the decorator itself
                if key not in ['fs', 'gcs_cache_dir']:
                    if isinstance(value, (str, int, float, bool)):
                        dynamic_name_parts.append(f"{key}_{str(value).replace(' ', '_').replace('/', '_')}")

            dynamic_name = "_".join(dynamic_name_parts)
            # Basic sanitization for GCS object names (replace common problematic chars)
            dynamic_name = dynamic_name.replace(':', '_').replace('*','_').replace('?','_')
            cache_filename = f"{dynamic_name}.{cache_file_type}"
            gcs_cache_path = os.path.join(gcs_cache_dir, cache_filename)

            # --- Check GCS Cache ---
            if gcs_path_exists(gcs_cache_path): # Uses helper from this file
                print(f"[GCS CACHE HIT] Found cache file: {gcs_cache_path}")
                try:
                    print(f"Attempting to load DataFrame from {gcs_cache_path}...")
                    with fs.open(gcs_cache_path, 'rb') as f: # Read in binary mode for parquet
                        if cache_file_type == 'csv':
                             df = pd.read_csv(f)
                        elif cache_file_type == 'parquet':
                             df = pd.read_parquet(f)
                    print(f"Successfully loaded '{dynamic_name}' from GCS cache.")
                    return df
                except Exception as e:
                    print(f"[GCS CACHE WARNING] Failed to load '{dynamic_name}' from GCS cache file {gcs_cache_path}: {e}")
                    print("Will recompute and attempt to overwrite corrupted cache file.")
                    # Attempt to delete the corrupted file before recomputing
                    delete_gcs_path(gcs_cache_path, recursive=False)

            # --- Cache Miss or Failed Load: Execute Function ---
            print(f"[GCS CACHE MISS] Running function '{func.__name__}' to generate '{dynamic_name}'...")
            result = func(*args, **kwargs) # Execute the actual function (e.g., BQ query)

            # --- Save Result to GCS Cache ---
            if isinstance(result, pd.DataFrame):
                if not result.empty:
                    print(f"[GCS CACHE SAVE] Saving result to GCS cache: {gcs_cache_path}")
                    try:
                        # Ensure parent GCS directory exists
                        parent_dir = os.path.dirname(gcs_cache_path)
                        if not gcs_path_exists(parent_dir):
                            print(f"Creating GCS cache directory: {parent_dir}")
                            fs.mkdirs(parent_dir, exist_ok=True)

                        # Write the file
                        with fs.open(gcs_cache_path, 'wb') as f: # Write in binary mode
                            if cache_file_type == 'csv':
                                result.to_csv(f, index=False, na_rep='NA')
                            elif cache_file_type == 'parquet':
                                # Requires 'pyarrow' or 'fastparquet' installed
                                result.to_parquet(f, index=False)
                        print(f"Successfully saved '{dynamic_name}' to GCS cache.")
                    except Exception as e:
                        print(f"[GCS CACHE WARNING] Failed to save '{dynamic_name}' to GCS cache at {gcs_cache_path}: {e}")
                        # Should we delete the failed write attempt? Maybe not, could be partial.
                else:
                    print(f"[GCS CACHE SKIP] Skipping GCS caching for empty DataFrame result: '{dynamic_name}'")
            elif result is not None:
                print(f"[GCS CACHE SKIP] Result type is {type(result)}, not a non-empty DataFrame. Skipping GCS caching.")

            return result
        return wrapper
    return decorator

# --- File Downloading ---

def download_file(url: str, local_gz_path: str, local_txt_path: str, log_id: str = "") -> str | None:
    """
    Downloads a file (potentially gzipped) from a URL and decompresses if needed.
    Checks if the final decompressed file already exists locally first.

    Args:
        url: The URL to download from.
        local_gz_path: The local path to save the downloaded (potentially gzipped) file.
        local_txt_path: The local path for the final (decompressed) file.
        log_id: An identifier (e.g., PRS ID) for log messages.

    Returns:
        The path to the final local (decompressed) file if successful, otherwise None.
    """
    log_prefix = f"[{log_id}] " if log_id else ""

    # Ensure local directory exists
    os.makedirs(os.path.dirname(local_gz_path), exist_ok=True)
    os.makedirs(os.path.dirname(local_txt_path), exist_ok=True)

    # Check if final target file already exists
    if os.path.exists(local_txt_path):
        print(f"{log_prefix}Using existing decompressed file: {local_txt_path}")
        return local_txt_path

    # If target doesn't exist, proceed with download
    print(f"{log_prefix}Downloading file from: {url} to {local_gz_path}")
    try:
        # Use stream=True for potentially large files
        response = requests.get(url, stream=True, timeout=600) # Increased timeout to 10 mins
        response.raise_for_status() # Raises HTTPError for bad responses (4xx or 5xx)

        with open(local_gz_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192*4): # Larger chunk size
                f.write(chunk)
        print(f"{log_prefix}Download complete: {local_gz_path}")

        # Decompress if it's a gzipped file (check extension)
        if url.endswith('.gz') or local_gz_path.endswith('.gz'):
            print(f"{log_prefix}Decompressing {local_gz_path} to {local_txt_path}...")
            with gzip.open(local_gz_path, 'rb') as f_in:
                with open(local_txt_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(local_gz_path) # Remove the gzipped file after successful decompression
            print(f"{log_prefix}Decompression complete. Final file: {local_txt_path}")
            return local_txt_path
        else:
            # If not gzipped, rename the downloaded file to the target txt path
            print(f"{log_prefix}File is not gzipped. Renaming {local_gz_path} to {local_txt_path}")
            shutil.move(local_gz_path, local_txt_path)
            return local_txt_path

    except requests.exceptions.Timeout:
        print(f"{log_prefix}ERROR: Timeout occurred while downloading {url}")
        return None
    except requests.exceptions.RequestException as e:
        print(f"{log_prefix}ERROR: Failed to download {url}. Error: {e}")
        # Clean up potentially incomplete download
        if os.path.exists(local_gz_path): os.remove(local_gz_path)
        return None
    except Exception as e:
        print(f"{log_prefix}ERROR: Failed during download/decompression for {url}. Error: {e}")
        # Clean up local files on any other error
        if os.path.exists(local_gz_path): os.remove(local_gz_path)
        if os.path.exists(local_txt_path): os.remove(local_txt_path)
        return None
