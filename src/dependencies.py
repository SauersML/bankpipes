# src/dependencies.py
import sys
import os
import subprocess
import importlib
import logging

# Configure basic logging for this script
LOG_FMT = "%(asctime)s [%(levelname)s] %(filename)s:%(lineno)d – %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FMT, handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

# --- Configuration ---
PYSPARK_VERSION_TO_INSTALL = "3.5.0"
NEST_ASYNCIO_SPEC = "nest-asyncio>=1.5.8,<2"
PYTHON_EXECUTABLE = sys.executable # pip uses the same python

# --- Functions ---
def install_package(package_name_with_version: str) -> bool:
    """Installs a single package using pip and returns True on success."""
    log.info(f"Attempting to install: {package_name_with_version}...")
    
    install_command = [
        PYTHON_EXECUTABLE, "-m", "pip", "install", 
        "--no-cache-dir", # Avoid using cached packages that might be problematic
        "--disable-pip-version-check", # Suppress pip version warnings
        package_name_with_version
    ]
    log.info(f"Executing: {' '.join(install_command)}")

    try:
        process_result = subprocess.run(
            install_command,
            capture_output=True,
            text=True,
            check=False # We check returncode manually for better logging
        )

        log.info(f"--- pip install output for {package_name_with_version} ---")
        if process_result.stdout:
            log.info("Pip STDOUT:\n%s", process_result.stdout.strip())
        if process_result.stderr:
            # Pip often uses stderr for normal download/install progress too,
            # so log as INFO unless returncode is non-zero.
            if process_result.returncode != 0:
                log.error("Pip STDERR:\n%s", process_result.stderr.strip())
            else:
                log.info("Pip STDERR (might contain progress/warnings):\n%s", process_result.stderr.strip())
        
        if process_result.returncode == 0:
            log.info(f"Pip command for {package_name_with_version} completed with exit code 0 (assumed success).")
            return True
        else:
            log.error(f"ERROR: Pip command for {package_name_with_version} FAILED with exit code {process_result.returncode}.")
            return False
            
    except FileNotFoundError:
        log.error(f"ERROR: Python executable '{PYTHON_EXECUTABLE}' or pip not found. Cannot install packages.")
        return False
    except Exception as e:
        log.error(f"ERROR: An unexpected error occurred during installation of {package_name_with_version}: {e}", exc_info=True)
        return False

def verify_pyspark(expected_version: str) -> bool:
    """Verifies PySpark installation and checks its version."""
    log.info("Verifying PySpark installation...")
    try:
        # Invalidate import caches to help pick up newly installed versions
        importlib.invalidate_caches()
        
        # Attempt to import pyspark
        log.info("Attempting to import 'pyspark' module...")
        import pyspark
        
        pyspark_location = os.path.dirname(pyspark.__file__)
        pyspark_version_actual = pyspark.__version__
        
        log.info(f"Successfully imported pyspark from: {pyspark_location}")
        log.info(f"Detected PySpark version: {pyspark_version_actual}")

        if pyspark_version_actual == expected_version:
            log.info(f"PySpark version {pyspark_version_actual} MATCHES target {expected_version}.")
            return True
        else:
            log.warning(f"WARNING: Detected PySpark version '{pyspark_version_actual}' DOES NOT MATCH target '{expected_version}'.")
            log.warning("This might be due to existing conflicting installations, pip resolver behavior, or issues with the install process.")
            # Depending on strictness, this could be treated as a failure.
            # For this critical dependency, we'll consider it a failure for the setup script's success.
            return False
            
    except ModuleNotFoundError:
        log.error("CRITICAL ERROR: PySpark module NOT FOUND after installation attempt.")
        return False
    except Exception as e:
        log.error(f"CRITICAL ERROR: An unexpected error occurred during PySpark verification: {e}", exc_info=True)
        return False

def verify_nest_asyncio() -> bool:
    """Verifies nest_asyncio installation."""
    log.info("Verifying nest_asyncio installation...")
    try:
        importlib.invalidate_caches()
        log.info("Attempting to import 'nest_asyncio' module...")
        import nest_asyncio
        nest_asyncio_location = os.path.dirname(nest_asyncio.__file__)
        log.info(f"Successfully imported nest_asyncio from: {nest_asyncio_location}")
        # Version check for nest_asyncio is less critical than exact pyspark version for this problem
        return True
    except ModuleNotFoundError:
        log.warning("WARNING: nest_asyncio module not found after installation attempt. This might cause issues for some Hail/GCS operations.")
        return False # Consider this non-fatal for the overall setup success, but log it
    except Exception as e:
        log.error(f"ERROR: An error occurred during nest_asyncio verification: {e}", exc_info=True)
        return False

# --- Main Execution ---
if __name__ == "__main__":
    log.info(f"--- Starting Dependency Setup using Python: {PYTHON_EXECUTABLE} ---")

    pyspark_installed_correctly = install_package(f"pyspark=={PYSPARK_VERSION_TO_INSTALL}")
    nest_asyncio_installed_correctly = install_package(NEST_ASYNCIO_SPEC)
    
    pyspark_verified = False
    if pyspark_installed_correctly:
        pyspark_verified = verify_pyspark(PYSPARK_VERSION_TO_INSTALL)
    else:
        log.error("Skipping PySpark verification because installation command failed or indicated failure.")

    # We can verify nest_asyncio even if its install_package call returned True but pip had other non-fatal warnings
    nest_asyncio_verified = verify_nest_asyncio()

    log.info("--- Dependency Setup Script Finished ---")

    if not pyspark_verified: # PySpark verification is the most critical part
        log.error("CRITICAL FAILURE: PySpark was not installed or verified to the target version.")
        sys.exit(1) # Exit with error code if PySpark isn't correctly set up
    
    if not nest_asyncio_verified: # nest_asyncio is important but perhaps not immediately fatal if missing
        log.warning("Nest-asyncio was not successfully verified. This might lead to issues later.")
        # Not exiting with error for nest_asyncio failure to allow pipeline to proceed if pyspark is okay

    log.info("Dependency setup deemed successful (PySpark verified to target version).")
    sys.exit(0)
