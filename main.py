import argparse
import datetime
import os
import subprocess
import sys
import yaml
import pandas as pd
import json
import shutil
import time
import psutil
import logging

# --- Global Constants and Default Settings ---

# Default filename for the CSV listing PRS models, in the same directory as main.py
MODELS_CSV_FILENAME = "models.csv"

# Default Python executable for running src scripts
PYTHON_EXECUTABLE = "/opt/conda/bin/python3"

# Suffix for the flagged samples file relative to the CDR_STORAGE_PATH environment variable
FLAGGED_SAMPLES_FILE_SUFFIX = "/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv"

# GCS directory suffixes (will be prepended with gs://<WORKSPACE_BUCKET_NAME>/ by PipelineConfig)
GCS_REUSABLE_INTERMEDIATES_SUFFIX = "pgs_pipeline_python_orchestrator/intermediates"
GCS_HAIL_TEMP_RUN_SPECIFIC_SUFFIX = "pgs_pipeline_python_orchestrator/hail_temp" # run_timestamp will be appended by PipelineConfig
GCS_RUN_OUTPUTS_SUFFIX = "pgs_pipeline_python_orchestrator/runs" # run_timestamp will be appended by PipelineConfig

# VDS Preparation Default Settings
VDS_PREP_ENABLE_DOWNSAMPLING = True
VDS_PREP_N_CASES_DOWNSAMPLE = 500
VDS_PREP_N_CONTROLS_DOWNSAMPLE = 500
VDS_PREP_DOWNSAMPLING_RANDOM_STATE = 2025

# --- Logging Setup ---
# Configures a logger to output messages to stdout.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] - %(module)s - %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


class PipelineConfig:
    """
    Manages pipeline configuration by loading from a user-provided YAML file,
    deriving values from environment variables, and applying hardcoded defaults.
    """
    def __init__(self, config_yaml_path):
        logger.info(f"Initializing pipeline configuration from: {config_yaml_path}")
        try:
            with open(config_yaml_path, 'r') as f:
                user_config = yaml.safe_load(f)
        except FileNotFoundError:
            logger.error(f"FATAL: Configuration file not found at {config_yaml_path}")
            sys.exit(1)
        except yaml.YAMLError as e:
            logger.error(f"FATAL: Error parsing YAML configuration file {config_yaml_path}: {e}")
            sys.exit(1)

        if not user_config or 'phenotype_definition' not in user_config:
            logger.error("FATAL: 'phenotype_definition' section is missing in the configuration file.")
            sys.exit(1)

        # --- User-specified phenotype definition from YAML ---
        self.phenotype_target_name = user_config['phenotype_definition'].get('target_name')
        self.phenotype_concept_ids = user_config['phenotype_definition'].get('concept_ids')

        if not self.phenotype_target_name:
            logger.error("FATAL: 'target_name' is missing under 'phenotype_definition' in the configuration file.")
            sys.exit(1)
        if not self.phenotype_concept_ids:
            logger.error("FATAL: 'concept_ids' are missing under 'phenotype_definition' in the configuration file.")
            sys.exit(1)
        if not isinstance(self.phenotype_concept_ids, list):
            logger.error("FATAL: 'concept_ids' must be a list in the configuration file.")
            sys.exit(1)
        logger.info(f"  Phenotype Target Name: {self.phenotype_target_name}")
        logger.info(f"  Phenotype Concept IDs: {self.phenotype_concept_ids}")

        # --- Load from Environment Variables (Mandatory for core functionality) ---
        self.google_billing_project = os.getenv('GOOGLE_PROJECT')
        self.workspace_bucket_name_raw = os.getenv('WORKSPACE_BUCKET') # Expected format: gs://bucket-name or bucket-name
        self.aou_wgs_vds_path = os.getenv('WGS_VDS_PATH')
        self.aou_workspace_cdr = os.getenv('WORKSPACE_CDR')
        self.cdr_storage_path_base = os.getenv('CDR_STORAGE_PATH')

        # Validate mandatory environment variables
        env_vars_to_check = {
            "GOOGLE_PROJECT": self.google_billing_project,
            "WORKSPACE_BUCKET": self.workspace_bucket_name_raw,
            "WGS_VDS_PATH": self.aou_wgs_vds_path,
            "WORKSPACE_CDR": self.aou_workspace_cdr,
            "CDR_STORAGE_PATH": self.cdr_storage_path_base
        }
        for var_name, value in env_vars_to_check.items():
            if not value:
                logger.error(f"FATAL: Required environment variable '{var_name}' is not set.")
                sys.exit(1)
        
        # Process workspace bucket name to remove 'gs://' prefix if present
        if self.workspace_bucket_name_raw.startswith("gs://"):
            self.workspace_bucket_name = self.workspace_bucket_name_raw[5:]
        else:
            self.workspace_bucket_name = self.workspace_bucket_name_raw

        # --- Apply Hardcoded Defaults & Construct Derived Paths ---
        self.run_timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
        self.script_dir = os.path.dirname(os.path.abspath(__file__)) # Directory where main.py resides
        
        self.models_csv_path = os.path.join(self.script_dir, MODELS_CSV_FILENAME)
        if not os.path.isfile(self.models_csv_path): # Check if it's a file
            logger.warning(f"Models CSV file not found at the default location: {self.models_csv_path}. Processing will fail if models are expected.")

        self.python_executable = PYTHON_EXECUTABLE
        self.aou_flagged_samples_gcs_path = f"{self.cdr_storage_path_base.rstrip('/')}{FLAGGED_SAMPLES_FILE_SUFFIX}"

        # Base GCS path for all pipeline operations, constructed from the workspace bucket name
        self.base_gcs_pipeline_path = f"gs://{self.workspace_bucket_name}"

        self.gcs_reusable_intermediates_dir = f"{self.base_gcs_pipeline_path}/{GCS_REUSABLE_INTERMEDIATES_SUFFIX.strip('/')}"
        self.gcs_run_hail_temp_dir = f"{self.base_gcs_pipeline_path}/{GCS_HAIL_TEMP_RUN_SPECIFIC_SUFFIX.strip('/')}/{self.run_timestamp}"
        self.gcs_run_outputs_dir = f"{self.base_gcs_pipeline_path}/{GCS_RUN_OUTPUTS_SUFFIX.strip('/')}/{self.run_timestamp}"

        # VDS Preparation Settings from global constants
        self.vds_prep_enable_downsampling = VDS_PREP_ENABLE_DOWNSAMPLING
        self.vds_prep_n_cases = VDS_PREP_N_CASES_DOWNSAMPLE
        self.vds_prep_n_controls = VDS_PREP_N_CONTROLS_DOWNSAMPLE
        self.vds_prep_random_state = VDS_PREP_DOWNSAMPLING_RANDOM_STATE

        # Dataproc specific configurations
        self.dataproc_region = os.getenv("DATAPROC_REGION", "us-central1") # Default Dataproc region
        self.dataproc_cluster_name = f"prs-pipeline-{self.run_timestamp}" # Unique cluster name for this run

        # Prepare Spark configurations JSON string, replacing internal template marker
        spark_conf_list_template = [
            # --- GCS Requester Pays and Project ID  ---
            "spark.hadoop.fs.gs.requester.pays.mode=AUTO",
            "spark.hadoop.fs.gs.requester.pays.project.id={{GOOGLE_BILLING_PROJECT}}",
            "spark.hadoop.fs.gs.project.id={{GOOGLE_BILLING_PROJECT}}", # General GCS project ID

            # --- GCS Connector Performance ---
            # Larger block size for GCS operations, perhaps improves read/write performance for large files.
            "spark.hadoop.fs.gs.block.size=134217728", # 128MB

            # --- Spark Dynamic Allocation ---
            # They will be active when jobs are submitted to Dataproc.
            "spark.dynamicAllocation.enabled=true",

            # --- Memory Management and Serialization ---
            # Kryo is generally faster than Java serialization for Spark.
            "spark.serializer=org.apache.spark.serializer.KryoSerializer",

            # --- Garbage Collection (More relevant for executors on a cluster) ---
            # Using G1GC. InitiatingHeapOccupancyPercent can tune when G1GC kicks in.
            # These apply to the JVMs of the executors.
            "spark.executor.defaultJavaOptions=-XX:+UseG1GC -XX:InitiatingHeapOccupancyPercent=35 -XX:+PrintGCDetails -XX:+PrintGCTimeStamps"
        ]
        processed_spark_list = [
            s.replace("{{GOOGLE_BILLING_PROJECT}}", self.google_billing_project) for s in spark_conf_list_template
        ]
        dynamic_shuffle_partitions = max(200, (os.cpu_count() or 1) * 4)
        processed_spark_list.append(f"spark.sql.shuffle.partitions={dynamic_shuffle_partitions}")
        self.spark_conf_json_str = json.dumps(processed_spark_list)

        logger.info("Pipeline configuration loaded and processed:")
        logger.info(f"  Run Timestamp: {self.run_timestamp}")
        logger.info(f"  Google Billing Project: {self.google_billing_project}")
        logger.info(f"  Workspace Bucket: gs://{self.workspace_bucket_name}")
        logger.info(f"  Reusable Intermediates GCS Dir: {self.gcs_reusable_intermediates_dir}")
        logger.info(f"  Run Outputs GCS Dir: {self.gcs_run_outputs_dir}")
        logger.info(f"  Run Hail Temp GCS Dir: {self.gcs_run_hail_temp_dir}")
        logger.info(f"  Models CSV Path: {self.models_csv_path}")

def execute_script(script_name_in_src, command_args, config, log_dir_local, step_log_identifier, cwd_for_script=None):
    """
    Executes a Python script located in the 'src' directory using a subprocess.
    Manages environment variables (PYSPARK_PYTHON, PYTHONPATH) for the subprocess.
    Logs stdout and stderr to a specified file.
    Provides periodic status updates including new log lines or system resource usage.

    Args:
        script_name_in_src (str): Filename of the script in the 'src' directory (e.g., "fetch_phenotypes.py").
        command_args (list): A list of command-line arguments to pass to the script.
        config (PipelineConfig): The pipeline's configuration object.
        log_dir_local (str): The local directory where log files for script executions will be stored.
        step_log_identifier (str): A unique string to identify this execution step in log filenames.
        cwd_for_script (str, optional): The working directory from which to run the script.
                                       Defaults to the current working directory of main.py.

    Returns:
        str: The full path to the log file generated for this script execution.

    Raises:
        SystemExit: If the script is not found or if the script execution fails (non-zero exit code).
    """
    script_full_path = os.path.join(config.script_dir, "src", script_name_in_src)
    if not os.path.isfile(script_full_path):
        logger.error(f"FATAL: [{step_log_identifier}] Script file not found at {script_full_path}")
        sys.exit(1)

    command = [config.python_executable, script_full_path] + command_args
    log_file_name = f"{step_log_identifier}_{config.run_timestamp}.log"
    log_file_path = os.path.join(log_dir_local, log_file_name)

    logger.info(f"[{step_log_identifier}] Executing command: {' '.join(command)}")
    logger.info(f"[{step_log_identifier}] Log file for this step: {log_file_path}")
    if cwd_for_script:
        logger.info(f"[{step_log_identifier}] Setting working directory for script execution: {cwd_for_script}")

    process_env = os.environ.copy()
    process_env["PYSPARK_PYTHON"] = config.python_executable
    src_directory_path = os.path.join(config.script_dir, "src")
    current_pythonpath = process_env.get('PYTHONPATH', '')
    if current_pythonpath:
        process_env['PYTHONPATH'] = f"{src_directory_path}{os.pathsep}{current_pythonpath}"
    else:
        process_env['PYTHONPATH'] = src_directory_path
    logger.debug(f"[{step_log_identifier}] Subprocess PYTHONPATH set to: {process_env['PYTHONPATH']}")

    process = None
    last_log_line_count = 0
    monitoring_interval_seconds = 15

    try:
        # log directory exists before opening log file for subprocess
        os.makedirs(os.path.dirname(log_file_path), exist_ok=True)

        with open(log_file_path, 'w') as log_file_handle:
            process = subprocess.Popen(
                command,
                env=process_env,
                stdout=log_file_handle,
                stderr=subprocess.STDOUT,
                text=True,
                cwd=cwd_for_script,
            )
        logger.info(f"[{step_log_identifier}] Subprocess '{script_name_in_src}' started (PID: {process.pid}). Monitoring progress every {monitoring_interval_seconds} seconds...")

        while process.poll() is None: # While subprocess is running
            time.sleep(monitoring_interval_seconds)

            new_lines_printed_this_interval = False
            try:
                if os.path.exists(log_file_path):
                    with open(log_file_path, 'r') as current_log_read_handle:
                        all_lines = current_log_read_handle.readlines()
                    current_line_count = len(all_lines)

                    if current_line_count > last_log_line_count:
                        for i in range(last_log_line_count, current_line_count):
                            logger.info(f"[{step_log_identifier} LOG]: {all_lines[i].strip()}")
                        last_log_line_count = current_line_count
                        new_lines_printed_this_interval = True
                else:
                    # Log file might not have been created yet if script is very fast or stdout buffer not flushed
                    pass
            except Exception as log_read_e:
                logger.warning(f"[{step_log_identifier}] Error reading log file {log_file_path} for progress: {log_read_e}")

            if not new_lines_printed_this_interval:
                try:
                    # Get system-wide CPU and RAM usage.
                    cpu_usage = psutil.cpu_percent(interval=0.1) # Brief blocking call for current usage
                    ram_info = psutil.virtual_memory()
                    ram_usage_percent = ram_info.percent
                    logger.info(f"[{step_log_identifier}] Still running '{script_name_in_src}' (PID: {process.pid}). No new log output. System CPU: {cpu_usage}%, System RAM: {ram_usage_percent}% used.")
                except NameError:
                    logger.warning(f"[{step_log_identifier}] psutil not available. Cannot report CPU/RAM. Script '{script_name_in_src}' (PID: {process.pid}) is still running.")
                except Exception as psutil_e:
                    logger.warning(f"[{step_log_identifier}] Could not retrieve system stats: {psutil_e}. Script '{script_name_in_src}' (PID: {process.pid}) is still running.")
    
        # Subprocess finished, get final return code
        return_code = process.wait() # process resources are cleaned up and get final code.

        # Capture any final log lines written between the last check and process termination
        try:
            if os.path.exists(log_file_path):
                with open(log_file_path, 'r') as final_log_read_handle:
                    all_lines = final_log_read_handle.readlines()
                current_line_count = len(all_lines)
                if current_line_count > last_log_line_count:
                    logger.info(f"[{step_log_identifier}] Capturing final log lines from '{script_name_in_src}':")
                    for i in range(last_log_line_count, current_line_count):
                        logger.info(f"[{step_log_identifier} LOG]: {all_lines[i].strip()}")
        except Exception as final_log_read_e:
            logger.warning(f"[{step_log_identifier}] Error reading final log lines from {log_file_path}: {final_log_read_e}")
        
        # Process results based on return code
        if return_code != 0:
            logger.error(f"FATAL: [{step_log_identifier}] Script '{script_name_in_src}' failed with exit code {return_code}.")
            logger.error(f"FATAL: [{step_log_identifier}] Full log available at: {log_file_path}")
            if os.path.exists(log_file_path):
                try:
                    with open(log_file_path, 'r') as opened_log_file:
                        log_lines = opened_log_file.readlines()
                    number_of_log_lines = len(log_lines)
                    if number_of_log_lines > 0:
                        sys.stderr.write(f"\nINFO: [{step_log_identifier}] Displaying content from log file '{os.path.basename(log_file_path)}':\n")
                        sys.stderr.write("-" * 80 + "\n")
                        if number_of_log_lines > 100:
                            sys.stderr.write(f"--- First 50 lines of {os.path.basename(log_file_path)} ---\n")
                            for line_content_data in log_lines[:50]: sys.stderr.write(line_content_data)
                            sys.stderr.write("...\n[Log content truncated - full log has " + str(number_of_log_lines) + " lines]\n...\n")
                            sys.stderr.write(f"--- Last 50 lines of {os.path.basename(log_file_path)} ---\n")
                            for line_content_data in log_lines[-50:]: sys.stderr.write(line_content_data)
                        else:
                            sys.stderr.write(f"--- Full content of {os.path.basename(log_file_path)} (Total {number_of_log_lines} lines) ---\n")
                            for line_content_data in log_lines: sys.stderr.write(line_content_data)
                        sys.stderr.write("-" * 80 + "\n")
                        sys.stderr.write(f"--- End of displayed log content from {os.path.basename(log_file_path)} ---\n\n")
                    else:
                        logger.warning(f"WARNING: [{step_log_identifier}] Log file '{log_file_path}' was found but is empty.")
                except IOError as log_display_err:
                    logger.error(f"ERROR: [{step_log_identifier}] Could not read log file '{log_file_path}' for display after script failure: {log_display_err}")
            else:
                logger.warning(f"WARNING: [{step_log_identifier}] Log file '{log_file_path}' not found after script failure.")
            sys.exit(return_code)
        
        logger.info(f"INFO: [{step_log_identifier}] Script '{script_name_in_src}' completed successfully with exit code {return_code}.")
        return log_file_path

    except FileNotFoundError:
        logger.error(f"FATAL: [{step_log_identifier}] Python executable '{config.python_executable}' not found. Cannot run script '{script_name_in_src}'.")
        sys.exit(1)
    except OSError as e:
        logger.error(f"FATAL: [{step_log_identifier}] OS error while trying to run script '{script_name_in_src}': {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"FATAL: [{step_log_identifier}] An unexpected error occurred while managing subprocess for '{script_name_in_src}': {e}")
        logger.error(f"FATAL: [{step_log_identifier}] Check log for any partial output: {log_file_path}")
        sys.exit(1)
    finally:
        # the subprocess is terminated if the orchestrator exits unexpectedly
        if process and process.poll() is None:
            logger.warning(f"[{step_log_identifier}] Orchestrator is exiting; attempting to terminate subprocess '{script_name_in_src}' (PID: {process.pid}).")
            try:
                process.terminate() # Send SIGTERM
                process.wait(timeout=10) # Wait up to 10 seconds for graceful termination
                logger.info(f"[{step_log_identifier}] Subprocess '{script_name_in_src}' (PID: {process.pid}) terminated.")
            except subprocess.TimeoutExpired:
                logger.warning(f"[{step_log_identifier}] Subprocess '{script_name_in_src}' (PID: {process.pid}) did not terminate gracefully after 10s. Sending SIGKILL.")
                process.kill() # Force kill
                process.wait() # Wait for kill to complete
                logger.info(f"[{step_log_identifier}] Subprocess '{script_name_in_src}' (PID: {process.pid}) killed.")
            except Exception as term_err:
                logger.error(f"[{step_log_identifier}] Error during subprocess termination for '{script_name_in_src}' (PID: {process.pid if process else 'unknown'}): {term_err}")

def run_pipeline(config: PipelineConfig):
    """
    Main pipeline orchestration logic. Executes PRS analysis steps sequentially.
    """
    logger.info(f"Starting Python-based Polygenic Score (PGS) Pipeline. Run ID: {config.run_timestamp}")

    # Create a local working directory for this specific pipeline run.
    # This directory will store logs and any temporary files created by main.py.
    local_run_work_dir = os.path.join(os.getcwd(), f"pgs_pipeline_work_{config.run_timestamp}")
    try:
        os.makedirs(local_run_work_dir, exist_ok=True)
    except OSError as e:
        logger.error(f"FATAL: Could not create local working directory {local_run_work_dir}: {e}")
        sys.exit(1)
        
    local_script_log_dir = os.path.join(local_run_work_dir, "script_execution_logs")
    try:
        os.makedirs(local_script_log_dir, exist_ok=True)
    except OSError as e:
        logger.error(f"FATAL: Could not create local script log directory {local_script_log_dir}: {e}")
        sys.exit(1)

    logger.info(f"Local working directory for this run: {local_run_work_dir}")
    logger.info(f"Logs for individual script executions will be stored in: {local_script_log_dir}")

    # --- Step 1: Fetch Phenotype Cases ---
    logger.info("\n" + "="*30 + " Step 1: Fetching Phenotype Cases " + "="*30)
    # Define the GCS path where the phenotype CSV will be written by fetch_phenotypes.py
    gcs_phenotype_cases_csv_path = f"{config.gcs_reusable_intermediates_dir}/phenotype_data/{config.phenotype_target_name.replace(' ','_')}_cases.csv"
    
    fetch_pheno_args = [
        "--phenotype_name", config.phenotype_target_name,
        "--phenotype_concept_ids", ",".join(map(str, config.phenotype_concept_ids)),
        "--workspace_cdr", config.aou_workspace_cdr,
        "--output_phenotype_csv_gcs_path", gcs_phenotype_cases_csv_path,
        "--google_billing_project", config.google_billing_project
    ]
    execute_script(
        script_name_in_src="fetch_phenotypes.py", 
        command_args=fetch_pheno_args, 
        config=config, 
        log_dir_local=local_script_log_dir, 
        step_log_identifier="01_fetch_phenotypes",
        cwd_for_script=local_run_work_dir # Script writes no CWD-relative files needed by main.py
    )
    logger.info(f"Phenotype cases CSV is expected to be at GCS path: {gcs_phenotype_cases_csv_path}")

    # --- Step 2: Prepare Base VDS ---
    logger.info("\n" + "="*30 + " Step 2: Preparing Base VDS " + "="*30)
    # Define GCS paths for the outputs of prepare_base_vds.py
    gcs_base_vds_path = f"{config.gcs_reusable_intermediates_dir}/base_cohort_vds/base_cohort_wgs_ehr_unrelated.vds"
    gcs_wgs_ehr_ids_path = f"{config.gcs_reusable_intermediates_dir}/cohort_definitions/people_with_WGS_EHR_ids_for_vds_filtering.csv"

    prepare_vds_args = [
        "--workspace_cdr", config.aou_workspace_cdr,
        "--run_timestamp", config.run_timestamp, # For Hail log naming within prepare_base_vds.py
        "--gcs_temp_dir", config.gcs_reusable_intermediates_dir, # Base for stable VDS/ID list checkpoints
        "--gcs_hail_temp_dir", config.gcs_run_hail_temp_dir, # Run-specific Hail temporary directory
        "--wgs_vds_path", config.aou_wgs_vds_path,
        "--flagged_samples_gcs_path", config.aou_flagged_samples_gcs_path,
        "--base_cohort_vds_path_out", gcs_base_vds_path,
        "--wgs_ehr_ids_gcs_path_out", gcs_wgs_ehr_ids_path,
        "--target_phenotype_name", config.phenotype_target_name, # For logging during downsampling
        "--phenotype_cases_gcs_path_input", gcs_phenotype_cases_csv_path, # Output from Step 1
        "--n_cases_downsample", str(config.vds_prep_n_cases),
        "--n_controls_downsample", str(config.vds_prep_n_controls),
        "--downsampling_random_state", str(config.vds_prep_random_state),
        "--google_billing_project", config.google_billing_project,
        "--spark_configurations_json", config.spark_conf_json_str
    ]
    if config.vds_prep_enable_downsampling:
        prepare_vds_args.append("--enable_downsampling_for_vds")
    
    execute_script(
        script_name_in_src="prepare_base_vds.py", 
        command_args=prepare_vds_args, 
        config=config, 
        log_dir_local=local_script_log_dir, 
        step_log_identifier="02_prepare_base_vds",
        cwd_for_script=local_run_work_dir # Script writes no CWD-relative files needed by main.py
    )
    logger.info(f"Base VDS GCS path: {gcs_base_vds_path}")
    logger.info(f"WGS+EHR IDs CSV GCS path: {gcs_wgs_ehr_ids_path}")

    # --- Step 3: Process PRS Models (Iterate through models in models.csv) ---
    logger.info("\n" + "="*30 + " Step 3: Processing PRS Models " + "="*30)
    if not os.path.isfile(config.models_csv_path):
        logger.error(f"FATAL: Models CSV file '{config.models_csv_path}' not found. Cannot proceed with model processing.")
        sys.exit(1)
        
    try:
        models_df = pd.read_csv(config.models_csv_path)
    except pd.errors.EmptyDataError:
        logger.error(f"FATAL: Models CSV file '{config.models_csv_path}' is empty.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"FATAL: Failed to read Models CSV file '{config.models_csv_path}': {e}")
        sys.exit(1)

    required_model_cols = ['id', 'url', 'phenotype']
    if not all(col in models_df.columns for col in required_model_cols):
        logger.error(f"FATAL: Models CSV file '{config.models_csv_path}' must contain the columns: {', '.join(required_model_cols)}.")
        sys.exit(1)

    all_local_summary_files = [] # To store paths of locally generated summary text files

    for index, model_row in models_df.iterrows():
        model_id = str(model_row['id']).strip()
        model_url = str(model_row['url']).strip()
        # This is the phenotype label associated with the PRS model itself (from the catalog)
        model_catalog_phenotype_label = str(model_row['phenotype']).strip() 
        
        logger.info(f"\n--- Processing Model: {model_id} (Catalog Pheno: {model_catalog_phenotype_label}, URL: {model_url}) ---")

        # --- Step 3a: Execute process_prs_model.py ---
        step_3a_log_id = f"03a_process_prs_model_{model_id}"
        # Define GCS output paths for this specific model and run
        model_specific_run_output_gcs_dir = f"{config.gcs_run_outputs_dir}/scores/{model_id}"
        
        final_hail_table_gcs_path = f"{model_specific_run_output_gcs_dir}/hail_table/{model_id}_scores.ht"
        final_score_csv_gcs_path = f"{model_specific_run_output_gcs_dir}/score_csv/{model_id}_scores.csv"

        process_model_args = [
            "--prs_id", model_id,
            "--prs_url", model_url,
            "--base_cohort_vds_path", gcs_base_vds_path, # From Step 2
            "--gcs_temp_dir", config.gcs_reusable_intermediates_dir, # For weights/intervals checkpoints
            "--gcs_hail_temp_dir", config.gcs_run_hail_temp_dir, # Run-specific Hail temp
            "--run_timestamp", config.run_timestamp, # For Hail log naming within the script
            "--output_final_hail_table_gcs_path", final_hail_table_gcs_path,
            "--output_final_score_csv_gcs_path", final_score_csv_gcs_path,
            "--google_billing_project", config.google_billing_project,
            "--spark_configurations_json", config.spark_conf_json_str
        ]
        execute_script(
            script_name_in_src="process_prs_model.py", 
            command_args=process_model_args, 
            config=config, 
            log_dir_local=local_script_log_dir, 
            step_log_identifier=step_3a_log_id,
            cwd_for_script=local_run_work_dir # Script writes no CWD-relative files needed by main.py
        )
        gcs_score_csv_for_this_model = final_score_csv_gcs_path # Store for the analysis sub-step
        logger.info(f"Score CSV for model {model_id} written to GCS: {gcs_score_csv_for_this_model}")

        # --- Step 3b: Execute analyze_one_model_results.py ---
        step_3b_log_id = f"03b_analyze_results_{model_id}"
        # This script writes plots and standardized scores to GCS directly.
        # It produces a local summary text file in its Current Working Directory.
        local_summary_filename_for_this_model = f"{model_id}_analysis_summary.txt" 
        
        analyze_results_args = [
            "--prs_id", model_id,
            "--prs_phenotype_label", model_catalog_phenotype_label, # PRS model's original phenotype
            "--score_csv_gcs_path", gcs_score_csv_for_this_model, # From Step 3a
            "--wgs_ehr_ids_gcs_path", gcs_wgs_ehr_ids_path, # From Step 2
            "--phenotype_cases_csv_gcs_path", gcs_phenotype_cases_csv_path, # From Step 1
            "--phenotype_name", config.phenotype_target_name, # The pipeline's target phenotype for this analysis run
            "--gcs_base_output_dir_run", config.gcs_run_outputs_dir, # Base for plots, standardized scores for this run
            "--run_timestamp", config.run_timestamp, # For unique local plot dir naming if used by script
            "--google_billing_project", config.google_billing_project,
            "--output_summary_file_name", local_summary_filename_for_this_model # Name of local file script will create
        ]
        execute_script(
            script_name_in_src="analyze_prs_results.py", 
            command_args=analyze_results_args, 
            config=config, 
            log_dir_local=local_script_log_dir, 
            step_log_identifier=step_3b_log_id,
            cwd_for_script=local_run_work_dir # summary file is created in a known, predictable location
        )
        # Collect path to the locally created summary file
        full_local_summary_path_for_model = os.path.join(local_run_work_dir, local_summary_filename_for_this_model)
        if os.path.isfile(full_local_summary_path_for_model):
            all_local_summary_files.append(full_local_summary_path_for_model)
            logger.info(f"Analysis summary for model {model_id} created locally at: {full_local_summary_path_for_model}")
        else:
            logger.warning(f"Analysis summary file for model {model_id} was expected at {full_local_summary_path_for_model} but was not found after script execution.")
    
    # --- Step 4: Collect and Finalize All Model Summaries ---
    logger.info("\n" + "="*30 + " Step 4: Collecting Analysis Summaries " + "="*30)
    if all_local_summary_files:
        collected_summary_local_filename = f"all_models_analysis_summaries_{config.run_timestamp}.txt"
        collected_summary_local_path = os.path.join(local_run_work_dir, collected_summary_local_filename)
        
        try:
            with open(collected_summary_local_path, 'w') as outfile_handle:
                for i, individual_summary_path in enumerate(all_local_summary_files):
                    if i > 0: # Add a separator between summaries from different models
                        outfile_handle.write("\n\n" + "="*80 + "\n\n")
                    try:
                        with open(individual_summary_path, 'r') as infile_handle:
                            outfile_handle.write(infile_handle.read())
                    except FileNotFoundError:
                        logger.warning(f"Could not find individual summary file {individual_summary_path} during collection phase.")
                    except IOError as e:
                        logger.warning(f"IOError reading individual summary file {individual_summary_path}: {e}")
            logger.info(f"Collected all model summaries into local file: {collected_summary_local_path}")

            # Upload the collected summary file to GCS in the run's output directory
            collected_summary_gcs_path = f"{config.gcs_run_outputs_dir}/analysis_summary_collection/{collected_summary_local_filename}"
            logger.info(f"Attempting to upload collected summary to GCS: {collected_summary_gcs_path}")
            try:
                from src.utils import get_gcs_fs 
                gcs_fs_instance = get_gcs_fs(project_id_for_billing=config.google_billing_project)
                
                gcs_summary_parent_dir = os.path.dirname(collected_summary_gcs_path)
                if not gcs_fs_instance.exists(gcs_summary_parent_dir):
                     logger.info(f"Creating GCS directory for collected summary: {gcs_summary_parent_dir}")
                     gcs_fs_instance.mkdirs(gcs_summary_parent_dir, exist_ok=True)
                gcs_fs_instance.put(collected_summary_local_path, collected_summary_gcs_path)
                logger.info(f"Collected summary successfully uploaded to GCS: {collected_summary_gcs_path}")
            except ImportError: # If 'from src.utils import get_gcs_fs' fails
                logger.error("Failed to import 'get_gcs_fs' from 'src.utils'. Cannot upload collected summary.")
            except Exception as e: # Catch other GCS upload errors
                logger.error(f"Failed to upload collected summary file to GCS. Error: {e}")
        except IOError as e:
            logger.error(f"IOError writing collected summary file {collected_summary_local_path}: {e}")
    else:
        logger.warning("No model-specific summary files were generated or found; skipping collection and upload.")

    # --- Pipeline Completion Summary ---
    logger.info("\n" + "="*50)
    logger.info(f"Python-based PGS Pipeline run '{config.run_timestamp}' finished.")
    logger.info(f"  Main run outputs GCS location: {config.gcs_run_outputs_dir}")
    logger.info(f"  Reusable intermediates GCS location (VDS, etc.): {config.gcs_reusable_intermediates_dir}")
    logger.info(f"  Local working directory (logs, temp files): {local_run_work_dir}")
    logger.info("Pipeline execution complete.")


if __name__ == "__main__":
    # Setup command line argument parsing for main.py itself
    cli_parser = argparse.ArgumentParser(description="PGS Pipeline Orchestrator.")
    cli_parser.add_argument(
        "--config",
        type=str,
        default="config.yaml",  # Default filename
        help="Path to the pipeline configuration YAML file. Defaults to 'config.yaml' in the same directory as main.py."
    )
    pipeline_args = cli_parser.parse_args()

    if not os.path.isabs(pipeline_args.config):
        config_file_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), pipeline_args.config)
    else:
        config_file_path = pipeline_args.config
    
    try:
        logger.info("Initializing pipeline configuration...")
        # Instantiate the configuration object, which loads YAML and environment variables
        pipeline_config = PipelineConfig(config_file_path) # Use the resolved path
        run_pipeline(pipeline_config)
    except SystemExit:
        # This handles sys.exit() calls from within the pipeline, often due to logged FATAL errors.
        logger.info("Pipeline exited due to a previously logged fatal error.")
    except Exception as e:
        # Catch any other unexpected exceptions at the top level of the pipeline.
        logger.exception(f"An unexpected FATAL error occurred at the top level of the pipeline: {e}")
        sys.exit(1) # Exit with a non-zero code to indicate failure
