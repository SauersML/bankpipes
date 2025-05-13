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
  Conda Env Mgmt       : Enabled (Use Mamba: ${params.use_mamba})
  Conda Cache          : ${params.conda_cache_dir}
============================================
""" // Updated summary

workflow.onStart {
    log.info "Starting PGS Analysis Workflow: ${workflow.runName}"
}

// ----------------------------------------------------------------------------
// Channels
// ----------------------------------------------------------------------------

Channel
    .fromPath(params.models_csv)
    .ifEmpty { error "Models CSV file not found or empty: ${params.models_csv}" }
    .splitCsv(header:true, sep:',')
    .map { row ->
        if (!row.id || !row.url || !row.phenotype) { // Added phenotype for plotting label
            error "Invalid row in models CSV: ${row}. Requires 'id', 'url', and 'phenotype' columns."
        }
        // Ensure the map structure is consistent for downstream processes
        [ id: row.id.trim(), url: row.url.trim(), phenotype_label: row.phenotype.trim() ]
    }
    .set { model_info_ch } // Channel emits: map [ id:'PGS...', url:'...', phenotype_label:'Ischemic Stroke' ]

model_info_ch.count().subscribe { count -> log.info "Found ${count} models to process from ${params.models_csv}" }


// ----------------------------------------------------------------------------
// Processes
// ----------------------------------------------------------------------------

process fetch_phenotype_cases {
    tag "Fetch Pheno: ${params.target_phenotype_name}"
    label 'python_bq_env'

    publishDir "${params.output_dir_base}/${workflow.runName}/logs/01_fetch_pheno", mode:'copy', overwrite:true, pattern: '*.{log,txt}'

    output:
    path "phenotype_cases_gcs_path.txt", emit: pheno_cases_gcs_path_file // Single file path output

    script:
    def gcs_pheno_out_path = "${params.gcs_temp_dir_base}/phenotype_data/${params.target_phenotype_name.replace(' ','_')}_cases.csv"
    """
    echo "INFO: Fetching phenotype data for '${params.target_phenotype_name}'..."
    # Ensure python is resolved from the Conda environment
    python ${projectDir}/src/fetch_phenotypes.py \\
        --phenotype_name                "${params.target_phenotype_name}" \\
        --phenotype_concept_ids         "${params.phenotype_concept_ids.join(',')}" \\
        --workspace_cdr                 "${System.getenv('WORKSPACE_CDR')}" \\
        --output_phenotype_csv_gcs_path "${gcs_pheno_out_path}" \\
        --project_bucket                "${System.getenv('WORKSPACE_BUCKET')}" \\
        > fetch_phenotypes.log 2>&1

    echo "${gcs_pheno_out_path}" > phenotype_cases_gcs_path.txt
    echo "INFO: Phenotype GCS path recorded in phenotype_cases_gcs_path.txt"
    """
}

process prepare_base_vds {
    tag "Prepare Base VDS"
    label 'pyhail_env'
    label 'spark_job'

    publishDir "${params.output_dir_base}/${workflow.runName}/logs/02_prepare_vds", mode:'copy', overwrite:true, pattern: '*.{log,txt}'

    input:
    path phenotype_cases_gcs_path_file // From fetch_phenotype_cases

    output:
    // *** CHANGED: Emit two separate, named file paths ***
    path "base_vds_gcs_path.txt",  emit: base_vds_path
    path "wgs_ehr_ids_gcs_path.txt", emit: wgs_ids_path

    script:
    def gcs_vds_out_path = "${params.gcs_temp_dir_base}/base_cohort_wgs_ehr_unrelated.vds"
    def gcs_ids_out_path = "${params.gcs_temp_dir_base}/cohort_definitions/people_with_WGS_EHR_ids_for_vds_filtering.csv"

    """
    echo "INFO: Preparing Base VDS..."
    # Read the GCS path from the input file into a shell variable
    pheno_input_gcs_path=\$(cat ${phenotype_cases_gcs_path_file})
    echo "INFO: Using Phenotype Cases CSV from GCS path: \${pheno_input_gcs_path}"

    # Ensure python is resolved from the Conda environment
    python ${projectDir}/src/prepare_base_vds.py \\
        --project_bucket                    "${System.getenv('WORKSPACE_BUCKET')}" \\
        --workspace_cdr                     "${System.getenv('WORKSPACE_CDR')}" \\
        --run_timestamp                     "${workflow.runName}" \\
        --gcs_temp_dir                      "${params.gcs_temp_dir_base}" \\
        --gcs_hail_temp_dir                 "${params.gcs_hail_temp_dir_base}" \\
        --wgs_vds_path                      "${params.wgs_vds_path}" \\
        --flagged_samples_gcs_path          "${params.flagged_samples_gcs_path}" \\
        --base_cohort_vds_path_out          "${gcs_vds_out_path}" \\
        --wgs_ehr_ids_gcs_path_out          "${gcs_ids_out_path}" \\
        --target_phenotype_name             "${params.target_phenotype_name}" \\
        --phenotype_cases_gcs_path_input    "\${pheno_input_gcs_path}" \\
        ${params.enable_downsampling_for_vds_generation ? '--enable_downsampling_for_vds' : ''} \\
        --n_cases_downsample                ${params.n_cases_downsample} \\
        --n_controls_downsample             ${params.n_controls_downsample} \\
        --downsampling_random_state         ${params.downsampling_random_state} \\
        > prepare_base_vds.log 2>&1

    # Write the output paths to the respective output files
    echo "${gcs_vds_out_path}" > base_vds_gcs_path.txt
    echo "${gcs_ids_out_path}" > wgs_ehr_ids_gcs_path.txt
    echo "INFO: Output GCS paths recorded."
    """
}

process process_prs_model {
    tag "Process PRS: ${model.id}"
    label 'pyhail_env'
    label 'spark_job'

    publishDir "${params.output_dir_base}/${workflow.runName}/scores/${model.id}/logs_calculation", mode:'copy', overwrite:true, pattern: '*.{log,txt}'

    input:
    // Input matches the first .combine() output structure
    tuple val(model), path(base_vds_gcs_path_file) // model is the map [id, url, phenotype_label]

    output:
    // Pass the model map and the new score file path
    tuple val(model), path("final_score_gcs_path.txt"), emit: model_and_score_file

    script:
    // Read the GCS path from the input file into a shell variable
    def base_vds_path = "\$(cat ${base_vds_gcs_path_file})"
    def run_specific_output_dir = "${params.output_dir_base}/${workflow.runName}"
    def final_hail_out_path = "${run_specific_output_dir}/scores/${model.id}/hail_table/${model.id}_scores.ht"
    def final_csv_out_path  = "${run_specific_output_dir}/scores/${model.id}/score_csv/${model.id}_scores.csv"

    """
    echo "INFO: Processing PRS model ${model.id}..."
    echo "INFO: Using Base VDS from GCS path: ${base_vds_path}"

    # Ensure python is resolved from the Conda environment
    python ${projectDir}/src/process_prs_model.py \\
        --prs_id                            "${model.id}" \\
        --prs_url                           "${model.url}" \\
        --base_cohort_vds_path              "${base_vds_path}" \\
        --gcs_temp_dir                      "${params.gcs_temp_dir_base}" \\
        --gcs_hail_temp_dir                 "${params.gcs_hail_temp_dir_base}" \\
        --run_timestamp                     "${workflow.runName}" \\
        --output_final_hail_table_gcs_path  "${final_hail_out_path}" \\
        --output_final_score_csv_gcs_path   "${final_csv_out_path}" \\
        > process_prs_model_${model.id}.log 2>&1

    echo "${final_csv_out_path}" > final_score_gcs_path.txt
    echo "INFO: Final score CSV GCS path for ${model.id} recorded."
    """
}

process analyze_one_model_results {
    tag "Analyze: ${model.id} for ${params.target_phenotype_name}"
    label 'pyhail_env'

    publishDir "${params.output_dir_base}/${workflow.runName}/scores/${model.id}/analysis_results", mode:'copy', overwrite:true, pattern: '*.{log,png,txt}'

    input:
    // *** CHANGED: Input tuple structure reflects the final .combine() operation ***
    tuple val(model), path(score_gcs_path_file), path(wgs_ehr_ids_gcs_path_file), path(phenotype_cases_gcs_path_file)

    output:
    path "${model.id}_analysis_summary.txt", emit: summary_file

    script:
    // Read the GCS paths from the input files into shell variables
    def score_csv_path = "\$(cat ${score_gcs_path_file})"
    def wgs_ehr_ids_path = "\$(cat ${wgs_ehr_ids_gcs_path_file})"
    def pheno_cases_path = "\$(cat ${phenotype_cases_gcs_path_file})"
    def run_specific_output_dir = "${params.output_dir_base}/${workflow.runName}"

    """
    echo "INFO: Analyzing results for PRS model ${model.id}..."
    echo "INFO: Score CSV GCS path: ${score_csv_path}"
    echo "INFO: WGS+EHR IDs GCS path: ${wgs_ehr_ids_path}"
    echo "INFO: Phenotype Cases GCS path: ${pheno_cases_path}"

    # Ensure python is resolved from the Conda environment
    python ${projectDir}/src/analyze_prs_results.py \\
        --prs_id                       "${model.id}" \\
        --prs_phenotype_label          "${model.phenotype_label}" \\
        --score_csv_gcs_path           "${score_csv_path}" \\
        --wgs_ehr_ids_gcs_path         "${wgs_ehr_ids_path}" \\
        --phenotype_cases_csv_gcs_path "${pheno_cases_path}" \\
        --phenotype_name               "${params.target_phenotype_name}" \\
        --gcs_base_output_dir_run      "${run_specific_output_dir}" \\
        --run_timestamp                "${workflow.runName}" \\
        --project_bucket               "${System.getenv('WORKSPACE_BUCKET')}" \\
        --output_summary_file_name     "${model.id}_analysis_summary.txt" \\
        > analyze_results_${model.id}.log 2>&1

    echo "INFO: Analysis for ${model.id} finished."
    """
}


// ----------------------------------------------------------------------------
// Workflow Execution (Revised Structure for v22.04.5 Compatibility)
// ----------------------------------------------------------------------------

workflow {

    // Stage 1: Fetch phenotype data (single output)
    // Call the process and explicitly capture its single output channel
    fetch_phenotype_cases()
    pheno_path_ch = fetch_phenotype_cases.out.pheno_cases_gcs_path_file

    // Stage 2: Prepare VDS (takes single pheno path, produces two separate single paths)
    // Call the process, passing the channel from Stage 1
    prepare_base_vds(pheno_path_ch)
    // Capture the two *distinct* named output channels
    base_vds_path_ch = prepare_base_vds.out.base_vds_path
    wgs_ids_path_ch  = prepare_base_vds.out.wgs_ids_path

    // Stage 3: Combine model info (multi-item) with base VDS path (single item)
    // This prepares the input for process_prs_model
    process_prs_inputs_ch = model_info_ch.combine(base_vds_path_ch)

    // Stage 4: Run PRS model processing using the combined channel
    process_prs_model(process_prs_inputs_ch)
    // Capture the output channel containing model info and score path
    processed_scores_ch = process_prs_model.out.model_and_score_file

    // Stage 5: Prepare inputs for analyze_one_model_results
    // Combine the multi-item PRS results with the single WGS IDs path channel
    analysis_inputs_intermediate_ch = processed_scores_ch.combine(wgs_ids_path_ch)
    // Combine the intermediate result with the single Pheno path channel
    final_analysis_inputs_ch = analysis_inputs_intermediate_ch.combine(pheno_path_ch)
        // This final channel has the structure:
        // [ [id, url, label], final_score_path.txt, wgs_ids_path.txt, pheno_path.txt ]
        // which matches the input definition of analyze_one_model_results

    // Stage 6: Run analysis using the final combined channel
    analyze_one_model_results(final_analysis_inputs_ch)
    // Capture the output channel containing summary file paths
    analysis_summary_files_ch = analyze_one_model_results.out.summary_file

    // Stage 7: Collect the summary files
    analysis_summary_files_ch.collectFile(
        name: 'all_analysis_summaries.txt',
        storeDir: "${params.output_dir_base}/${workflow.runName}/analysis_summary_collection",
        overwrite: true // Optional: Added for convenience during testing
    )
}

// --- Workflow Handlers ---

workflow.onComplete {
    def summary_file = "${params.output_dir_base}/${workflow.runName}/analysis_summary_collection/all_analysis_summaries.txt"
    // Ensure params are defined or provide defaults if they might be missing in some context
    def outputDir = params.output_dir_base ?: 'undefined_output_dir'
    def tempDir = params.gcs_temp_dir_base ?: 'undefined_temp_dir'
    def condaCache = params.conda_cache_dir ?: 'undefined_conda_cache'

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
      Run Output  : ${outputDir}/${workflow.runName}
      Temp Files  : ${tempDir}
      Conda Cache : ${condaCache}
      Collected analysis summaries (if successful): ${summary_file} (path depends on execution context)
    ============================================
    """.stripIndent()
}

workflow.onError {
    // Ensure params are defined or provide defaults if they might be missing
    def outputDir = params.output_dir_base ?: 'undefined_output_dir'

    log.error """
    ============================================
      PGS Analysis Workflow Error!
    ============================================
      Run Name    : ${workflow.runName}
      Status      : Failed
      Exit Status : ${workflow.exitStatus}
      Error report: ${workflow.errorReport ?: 'N/A'}
      Trace file  : ${outputDir}/${workflow.runName}/reports/trace.txt (if generated and output_dir_base defined)
      Check work directory for task logs: ${workflow.workDir}
      Check Nextflow log (${workflow.logFile ?: '.nextflow.log'}) for more details.
    ============================================
    """.stripIndent()
}
