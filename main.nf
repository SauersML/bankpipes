#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------------------------------------------------------------
// P G S   A N A L Y S I S   N E X T F L O W   P I P E L I N E
// ----------------------------------------------------------------------------
// Description: Calculates Polygenic Scores (PGS) using Hail and analyzes results.
// ----------------------------------------------------------------------------

log.info """
============================================
  Pipeline Configuration Summary
============================================
  Run Name             : ${workflow.runName}
  Launch Directory     : ${workflow.launchDir}
  Work Directory       : ${workflow.workDir}
  Project Directory    : ${projectDir}
  Output Base (GCS)    : ${params.output_dir_base}
  Temp Base (GCS)      : ${params.gcs_temp_dir_base}
  Hail Temp (GCS)      : ${params.gcs_hail_temp_dir_base}
  --- Input Data ---
  WGS VDS Path         : ${params.wgs_vds_path}
  Flagged Samples Path : ${params.flagged_samples_gcs_path}
  Models Source        : ${params.models_csv}
  --- Phenotype ---
  Target Phenotype     : ${params.target_phenotype_name}
  Concept IDs          : ${params.phenotype_concept_ids}
  --- VDS Preparation ---
  Enable Downsampling  : ${params.enable_downsampling_for_vds_generation}
  Cases (Downsample)   : ${params.n_cases_downsample}
  Controls (Downsample): ${params.n_controls_downsample}
  Downsampling Seed    : ${params.downsampling_random_state}
  --- Execution ---
  Executor             : ${workflow.executor}
  Containerization     : Disabled (Bare Metal)
  Max Spark Forks      : ${params.max_spark_forks}
============================================
"""

workflow.onStart {
    log.info "Starting PGS Analysis Workflow: ${workflow.runName}"
}

// ----------------------------------------------------------------------------
// Channels
// ----------------------------------------------------------------------------

// Create a channel emitting model information parsed from the input CSV file
Channel
    .fromPath(params.models_csv)
    .ifEmpty { error "Models CSV file not found or empty: ${params.models_csv}" }
    .splitCsv(header:true, sep:',') // Assumes comma separated with header
    .map { row ->
        // Basic validation
        if (!row.id || !row.url) {
            error "Invalid row in models CSV: ${row}. Requires 'id' and 'url' columns."
        }
        // Create a map for each model
        [ id: row.id.trim(), url: row.url.trim() ]
    }
    .set { model_channel } // Output channel: model_channel = [ [id:'PGS...', url:'...'], ... ]

// Log the number of models found
model_channel.count().subscribe { count -> log.info "Found ${count} models to process from ${params.models_csv}" }


// ----------------------------------------------------------------------------
// Processes
// ----------------------------------------------------------------------------

// 1) Fetch phenotype cases from BigQuery (run once)
//    Uses BQ client, Pandas. Needs Python environment.
process fetch_phenotype_cases {
    tag "Fetch Pheno: ${params.target_phenotype_name}"
    label 'pyhail_env' // Apply beforeScript for Python packages

    publishDir "${workflow.launchDir}/results/01_fetch_pheno", mode:'copy', overwrite:true, pattern: '*'

    output:
    path "phenotype_cases_gcs_path.txt", emit: pheno_gcs_path_file // File containing the GCS path to the CSV

    script:
    // Define the expected GCS output path for the phenotype CSV
    def gcs_pheno_out_path = "${params.gcs_temp_dir_base}/phenotype_cases/${params.target_phenotype_name.replace(' ','_')}_cases.csv"
    """
    echo "INFO: Fetching phenotype data for '${params.target_phenotype_name}'..."
    python3 ${projectDir}/src/fetch_phenotypes.py \\
        --phenotype_name                "${params.target_phenotype_name}" \\
        --phenotype_concept_ids         "${params.phenotype_concept_ids.join(',')}" \\
        --workspace_cdr                 "${System.getenv('WORKSPACE_CDR')}" \\
        --output_phenotype_csv_gcs_path "${gcs_pheno_out_path}" \\
        --project_bucket                "${System.getenv('WORKSPACE_BUCKET')}"

    # Record the GCS path to the output file for downstream processes
    echo "INFO: Phenotype data saved to GCS: ${gcs_pheno_out_path}"
    echo "${gcs_pheno_out_path}" > phenotype_cases_gcs_path.txt
    echo "INFO: Phenotype GCS path recorded in phenotype_cases_gcs_path.txt"
    """
}

// 2) Prepare Base VDS (run once)
//    Uses Hail/Spark, Pandas, BQ client. Needs Python+Hail env and Spark resource limits.
//    Takes phenotype CSV path as input if downsampling.
process prepare_base_vds {
    tag "Prepare Base VDS"
    label 'pyhail_env' // Apply beforeScript for Python/Hail packages
    label 'spark_job'  // Apply resource limits (e.g., maxForks)

    publishDir "${workflow.launchDir}/results/02_prepare_vds", mode:'copy', overwrite:true, pattern: '*.txt' // Publish logs or path files

    input:
    path phenotype_cases_gcs_path_file from fetch_phenotype_cases.out.pheno_gcs_path_file // Path to the file containing the GCS path of phenotype cases CSV

    output:
    tuple path("base_vds_path.txt"), path("wgs_ehr_ids_gcs_path.txt"), emit: base_vds_files // Files containing GCS paths

    script:
    // Define expected GCS output paths for VDS and sample IDs
    def gcs_vds_out_path = "${params.gcs_temp_dir_base}/base_cohort_wgs_ehr_unrelated.vds" // Checkpoint VDS
    def gcs_ids_out_path = "${params.gcs_temp_dir_base}/people_with_WGS_EHR_ids.csv"      // Sample IDs used for VDS filtering

    // Read the GCS path of the phenotype CSV from the input file
    def pheno_input_gcs_path_cmd = "cat ${phenotype_cases_gcs_path_file}"

    """
    echo "INFO: Preparing Base VDS..."
    # Retrieve the GCS path for the phenotype CSV needed for downsampling (if enabled)
    PHENO_INPUT_GCS_PATH=\$(${pheno_input_gcs_path_cmd})
    echo "INFO: Using Phenotype Cases CSV from GCS path: \${PHENO_INPUT_GCS_PATH}"

    python3 ${projectDir}/src/prepare_base_vds.py \\
        --project_bucket                    "${System.getenv('WORKSPACE_BUCKET')}" \\
        --workspace_cdr                     "${System.getenv('WORKSPACE_CDR')}" \\
        --run_timestamp                     "${workflow.runName}" \\
        --gcs_temp_dir                      "${params.gcs_temp_dir_base}" \\
        --gcs_hail_temp_dir                 "${params.gcs_hail_temp_dir_base}" \\
        --wgs_vds_path                      "${params.wgs_vds_path}" \\
        --flagged_samples_gcs_path          "${params.flagged_samples_gcs_path}" \\
        --base_cohort_vds_path_out          "${gcs_vds_out_path}" \\
        --wgs_ehr_ids_gcs_path_out          "${gcs_ids_out_path}" \\
        --phenotype_cases_gcs_path_input    "\${PHENO_INPUT_GCS_PATH}" \\
        ${params.enable_downsampling_for_vds_generation ? '--enable_downsampling_for_vds' : ''} \\
        --n_cases_downsample                ${params.n_cases_downsample} \\
        --n_controls_downsample             ${params.n_controls_downsample} \\
        --downsampling_random_state         ${params.downsampling_random_state}

    # Record the GCS paths to the generated VDS and sample ID list
    echo "INFO: Base VDS written to GCS: ${gcs_vds_out_path}"
    echo "${gcs_vds_out_path}" > base_vds_path.txt
    echo "INFO: WGS+EHR IDs written to GCS: ${gcs_ids_out_path}"
    echo "${gcs_ids_out_path}" > wgs_ehr_ids_gcs_path.txt
    echo "INFO: Output GCS paths recorded."
    """
}

// 3) Process each PRS model in parallel
//    Uses Hail/Spark, Pandas. Needs Python+Hail env and Spark resource limits.
process process_prs_model {
    tag "Process PRS: ${model.id}"
    label 'pyhail_env' // Apply beforeScript for Python/Hail packages
    label 'spark_job'  // Apply resource limits (e.g., maxForks)

    publishDir "${workflow.launchDir}/results/03_process_model/${model.id}", mode:'copy', overwrite:true, pattern: '*.txt'

    input:
    tuple val(model), path(base_vds_path_file) // model = [id:'PGS...', url:'...'], base_vds_path_file contains GCS path to VDS

    output:
    tuple val(model.id), path("final_score_gcs_path.txt"), emit: score_files // File containing GCS path to the final score CSV

    script:
    // Read the GCS path to the base VDS from the input file
    def base_vds_path_cmd = "cat ${base_vds_path_file}"
    // Define GCS output paths for this model's final scores (Hail table and CSV) within the run-specific output directory
    def run_output_dir = "${params.output_dir_base}/${workflow.runName}"
    def final_hail_out_path = "${run_output_dir}/scores/${model.id}/hail_table" // Final Hail Table destination
    def final_csv_out_path  = "${run_output_dir}/scores/${model.id}/score/${model.id}_scores.csv" // Final score CSV destination

    """
    echo "INFO: Processing PRS model ${model.id}..."
    BASE_VDS_PATH=\$(${base_vds_path_cmd})
    echo "INFO: Using Base VDS from GCS path: \${BASE_VDS_PATH}"

    python3 ${projectDir}/src/process_prs_model.py \\
        --prs_id                            "${model.id}" \\
        --prs_url                           "${model.url}" \\
        --base_cohort_vds_path              "\${BASE_VDS_PATH}" \\
        --gcs_temp_dir                      "${params.gcs_temp_dir_base}" \\
        --gcs_hail_temp_dir                 "${params.gcs_hail_temp_dir_base}" \\
        --run_timestamp                     "${workflow.runName}" \\
        --output_final_hail_table_gcs_path  "${final_hail_out_path}" \\
        --output_final_score_csv_gcs_path   "${final_csv_out_path}"

    # Record the GCS path to the final score CSV for the analysis step
    echo "INFO: Final scores for ${model.id} saved to GCS: ${final_csv_out_path}"
    echo "${final_csv_out_path}" > final_score_gcs_path.txt
    echo "INFO: Final score CSV GCS path recorded."
    """
}

// 4) Analyze all results (run once)
//    Uses Pandas, Scikit-learn, Matplotlib etc. Needs Python env.
process analyze_all_results {
    tag "Analyze Results"
    label 'pyhail_env' // Apply beforeScript for Python packages

    publishDir "${workflow.launchDir}/results/04_analysis", mode:'copy', overwrite:true, pattern: '*' // Publish summary log and plots

    input:
    path score_gcs_path_files from process_prs_model.out.score_files.map { it[1] }.collect() // Collect all "final_score_gcs_path.txt" files
    path wgs_ehr_ids_gcs_path_file from prepare_base_vds.out.base_vds_files.map { it[1] }    // The "wgs_ehr_ids_gcs_path.txt" file
    path phenotype_cases_gcs_path_file from fetch_phenotype_cases.out.pheno_gcs_path_file  // The "phenotype_cases_gcs_path.txt" file

    output:
    path "analysis_summary.log", emit: summary_log // Final summary log file

    script:
    // Define the output directory on GCS where the analysis script will save plots etc.
    def run_output_dir = "${params.output_dir_base}/${workflow.runName}"

    """
    echo "INFO: Starting analysis of all PRS results..."

    # Create a manifest file listing the GCS paths of all score CSVs
    echo "INFO: Building score manifest file..."
    > score_manifest.txt
    for f in ${score_gcs_path_files.join(' ')}; do
      if [[ -f "\$f" ]]; then
        cat "\$f" >> score_manifest.txt
        echo >> score_manifest.txt # Add newline separator
      else
         echo "WARN: Score path file not found: \$f"
      fi
    done
    echo "INFO: Score manifest created:"
    cat score_manifest.txt

    # Read the GCS paths for WGS/EHR IDs and Phenotype Cases from their respective files
    FINAL_WGS_IDS_GCS_PATH=\$(cat ${wgs_ehr_ids_gcs_path_file})
    PHENO_CASES_GCS_PATH=\$(cat ${phenotype_cases_gcs_path_file})
    echo "INFO: Using WGS+EHR IDs from GCS: \${FINAL_WGS_IDS_GCS_PATH}"
    echo "INFO: Using Phenotype Cases from GCS: \${PHENO_CASES_GCS_PATH}"

    # Run the analysis script which loops through the manifest
    echo "INFO: Executing analysis script..."
    python3 ${projectDir}/src/analyze_prs_results.py \\
        --score_manifest               score_manifest.txt \\
        --wgs_ehr_ids_gcs_path         "\${FINAL_WGS_IDS_GCS_PATH}" \\
        --phenotype_cases_csv_gcs_path "\${PHENO_CASES_GCS_PATH}" \\
        --phenotype_name               "${params.target_phenotype_name}" \\
        --gcs_base_output_dir_run      "${run_output_dir}" \\
        --run_timestamp                "${workflow.runName}" \\
        --project_bucket               "${System.getenv('WORKSPACE_BUCKET')}" \\
        --output_summary_log           analysis_summary.log

    echo "INFO: Analysis script finished. Summary log created: analysis_summary.log"
    """
}

// ----------------------------------------------------------------------------
// Workflow Wiring and Execution
// ----------------------------------------------------------------------------

workflow {
    // Step 1: Fetch phenotype cases definitions and get the GCS path to the CSV
    pheno_file_ch = fetch_phenotype_cases()

    // Step 2: Prepare the base VDS, filtering and optionally downsampling.
    // Needs the phenotype CSV path for downsampling.
    // Outputs paths to the filtered VDS and the list of included sample IDs.
    base_ch = prepare_base_vds(pheno_file_ch)
            // base_ch emits: tuple path("base_vds_path.txt"), path("wgs_ehr_ids_gcs_path.txt")

    // Step 3: Process each model from the model_channel in parallel.
    // Combines each model with the path to the base VDS.
    // Outputs the path to the final score CSV for each model.
    prs_ch = process_prs_model(
                model_channel.combine(base_ch.map { it[0] }) // Combine models with base_vds_path.txt
             )
            // prs_ch emits: tuple val(model.id), path("final_score_gcs_path.txt")

    // Step 4: Run the final analysis step once all models are processed.
    // Collects all score CSV paths, the sample ID list, and the phenotype case list path.
    analyze_all_results (
        prs_ch.map { it[1] }.collect(),  // Collect all "final_score_gcs_path.txt" files
        base_ch.map { it[1] },          // Pass the "wgs_ehr_ids_gcs_path.txt" file path
        pheno_file_ch                   // Pass the "phenotype_cases_gcs_path.txt" file path
    )
}

workflow.onComplete {
    log.info """
    ============================================
      PGS Analysis Workflow Completed
    ============================================
      Run Name    : ${workflow.runName}
      Status      : ${workflow.success ? 'Success' : 'Failed'}
      Exit Status : ${workflow.exitStatus}
      Duration    : ${workflow.duration}
      Work Dir    : ${workflow.workDir}
      Launch Dir  : ${workflow.launchDir}
      Output Dir  : ${params.output_dir_base}/${workflow.runName}
    ============================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    ============================================
      PGS Analysis Workflow Error!
    ============================================
      Run Name    : ${workflow.runName}
      Status      : Failed
      Exit Status : ${workflow.exitStatus}
      Error report: ${workflow.errorReport}
      Trace file  : ${workflow.launchDir}/results/reports/trace.txt
      Check work directory for task logs: ${workflow.workDir}
    ============================================
    """.stripIndent()
}
