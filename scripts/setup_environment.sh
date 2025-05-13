#!/bin/bash
#
# setup_environment.sh
#
# This script sets up the Python environment for Nextflow tasks.
# It expects arguments defining the python executable, the task label (if any),
# and all required package versions.

# Exit immediately if a command exits with a non-zero status.
set -e

# --- Argument Parsing ---
if [ "$#" -lt 18 ]; then
    echo "ERROR: Not enough arguments provided to setup_environment.sh"
    echo "Usage: $0 <python_exe> <task_label> <pandas_ver> <numpy_ver> <pyarrow_ver> <requests_ver> <gcsfs_ver> <gcloud_bq_ver> <db_dtypes_ver> <pyspark_ver> <hail_ver> <nest_asyncio_ver> <sklearn_ver> <scipy_ver> <matplotlib_ver> <seaborn_ver>"
    exit 1
fi

received_python_exe="$1"
task_label="$2" # Can be empty string '' if no specific label setup needed
version_pandas="$3"
version_numpy="$4"
version_pyarrow="$5"
version_requests="$6"
version_gcsfs="$7"
version_google_cloud_bigquery="$8"
version_db_dtypes="$9"
version_pyspark="${10}"
version_hail="${11}"
version_nest_asyncio="${12}"
version_scikit_learn="${13}"
version_scipy="${14}"
version_matplotlib="${15}"
version_seaborn="${16}"

# --- Environment Initialization ---
echo "INFO: Initializing Python environment (PID: $$)"
PYTHON_EXE="${received_python_exe}" # Use the python executable passed as an argument

# Validate Python Version
if ! ${PYTHON_EXE} -V 2>&1 | grep -q "Python 3"; then
    echo "ERROR: PYTHON_EXE (set to '${PYTHON_EXE}') does not appear to be Python 3."
    exit 1
fi
echo "INFO: Using Python: $(${PYTHON_EXE} --version)"
echo "INFO: Python executable path: $(command -v ${PYTHON_EXE} || echo 'Not in PATH')"
echo "INFO: PIP version: $(${PYTHON_EXE} -m pip --version)"

# Ensure ~/.local/bin is in PATH
if [[ ":$PATH:" != *":${HOME}/.local/bin:"* ]]; then
    export PATH="${HOME}/.local/bin:${PATH}"
    echo "INFO: Added ${HOME}/.local/bin to PATH for task."
fi

# --- Helper Functions for Installation ---
check_install_exact() {
    local package_name=$1
    local exact_version=$2
    local installed_version
    installed_version=$(${PYTHON_EXE} -m pip show ${package_name} 2>/dev/null | grep '^Version:' | awk '{print $2}')

    if [[ "$installed_version" == "$exact_version" ]]; then
        echo "INFO: Verified ${package_name} version ${installed_version} matches required ${exact_version}."
    else
        if [[ -n "$installed_version" ]]; then
             echo "WARN: Found ${package_name} version ${installed_version}, but require ${exact_version}. Reinstalling..."
        else
             echo "INFO: ${package_name} not found. Installing required version ${exact_version}..."
        fi
        if ${PYTHON_EXE} -m pip install --user --no-cache-dir --upgrade --force-reinstall "${package_name}==${exact_version}"; then
             echo "INFO: ${package_name}==${exact_version} install/update successful."
             local new_version
             new_version=$(${PYTHON_EXE} -m pip show ${package_name} 2>/dev/null | grep '^Version:' | awk '{print $2}')
             if [[ "$new_version" == "$exact_version" ]]; then
                 echo "INFO: Confirmed installed version: ${new_version}"
             else
                 echo "ERROR: Could not confirm ${package_name} version ${exact_version} (${new_version} found). Exiting."
                 exit 1
             fi
        else
             echo "ERROR: Failed to install ${package_name}==${exact_version}. Exiting."
             exit 1
        fi
    fi
}

check_install_other() {
    local package_spec=$1
    local pkg_name_for_check=$(echo "$package_spec" | sed -E 's/[<>=!].*//')

    if ${PYTHON_EXE} -m pip show ${pkg_name_for_check} > /dev/null 2>&1; then
        echo "INFO: Package ${pkg_name_for_check} found. Ensuring version spec '${package_spec}' is met..."
        if ! ${PYTHON_EXE} -m pip install --user --no-cache-dir --upgrade "${package_spec}"; then
            echo "ERROR: Failed updating package spec '${package_spec}'. Exiting."
            exit 1
        fi
         # No further action needed if update successful
    else
        echo "INFO: Package ${pkg_name_for_check} not found. Installing '${package_spec}'..."
        if ${PYTHON_EXE} -m pip install --user --no-cache-dir --upgrade "${package_spec}"; then
             echo "INFO: Install successful for '${package_spec}'"
        else
             echo "ERROR: Failed to install package spec '${package_spec}'. Exiting."
             exit 1
        fi
    fi
}

# --- Environment Variable Checks ---
echo "INFO: Verifying required environment variables WORKSPACE_BUCKET and WORKSPACE_CDR..."
# These are expected to be set by the Nextflow process environment
if [[ -z "${WORKSPACE_BUCKET:-}" ]]; then
    echo "ERROR: Environment variable WORKSPACE_BUCKET is not set or is empty."
    exit 1
fi
if [[ -z "${WORKSPACE_CDR:-}" ]]; then
    echo "ERROR: Environment variable WORKSPACE_CDR is not set or is empty."
    exit 1
fi
echo "INFO: Confirmed WORKSPACE_BUCKET=${WORKSPACE_BUCKET}"
echo "INFO: Confirmed WORKSPACE_CDR=${WORKSPACE_CDR}"
echo "----------------------------------------"

# --- Label-Specific Package Installation ---
# This section uses the task_label ($2) passed to the script

if [[ "$task_label" == "python_bq_env" ]]; then
    echo "INFO: Setting up 'python_bq_env' packages..."
    check_install_other "pandas==${version_pandas}"
    check_install_other "numpy==${version_numpy}"
    check_install_other "pyarrow==${version_pyarrow}"
    check_install_other "requests==${version_requests}"
    check_install_other "gcsfs==${version_gcsfs}"
    check_install_other "google-cloud-bigquery==${version_google_cloud_bigquery}"
    check_install_other "db-dtypes==${version_db_dtypes}"
    echo "INFO: 'python_bq_env' package setup complete."
    echo "----------------------------------------"

elif [[ "$task_label" == "pyhail_env" ]]; then
    echo "INFO: Setting up 'pyhail_env' (Spark/Hail + analysis tools) packages..."
    # Core Spark/Hail
    check_install_exact "pyspark" "${version_pyspark}"
    check_install_exact "hail" "${version_hail}"
    check_install_other "nest-asyncio>=${version_nest_asyncio},<2"
    # Common Data Science & GCP
    check_install_other "pandas==${version_pandas}"
    check_install_other "numpy==${version_numpy}"
    check_install_other "pyarrow==${version_pyarrow}"
    check_install_other "requests==${version_requests}"
    check_install_other "gcsfs==${version_gcsfs}"
    check_install_other "google-cloud-bigquery==${version_google_cloud_bigquery}"
    check_install_other "db-dtypes==${version_db_dtypes}"
    # Analysis & Viz
    check_install_other "scikit-learn==${version_scikit_learn}"
    check_install_other "scipy==${version_scipy}"
    check_install_other "matplotlib==${version_matplotlib}"
    check_install_other "seaborn==${version_seaborn}"
    echo "INFO: 'pyhail_env' package setup complete."
    echo "Python executable: $(${PYTHON_EXE} -V)"
    ${PYTHON_EXE} -m pip list | grep -E 'pyspark|hail|nest-asyncio|pandas|numpy|pyarrow|requests|gcsfs|google-cloud-bigquery|db-dtypes|scikit-learn|scipy|matplotlib|seaborn' || echo "WARN: Some key packages might not be listed by grep."
    echo "----------------------------------------"

else
    # Default case: No specific label, or a label we don't have special setup for.
    # The base setup (Python checks, PATH, env var checks) already ran.
    echo "INFO: No specific package installations requested for label: '${task_label}' (or label is empty)."
    echo "----------------------------------------"
fi

echo "INFO: Environment setup script finished successfully."
exit 0
