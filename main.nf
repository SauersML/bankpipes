#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------------------------------------------------------------
// P G S   A N A L Y S I S   N E X T F L O W   P I P E L I N E
// ----------------------------------------------------------------------------
// Description: Calculates Polygenic Scores (PGS) using Hail and analyzes results.
// Relies on a pre-configured Python environment as found on AoU Workbenches.
// ----------------------------------------------------------------------------

log.info """
============================================
  Pipeline Configuration Summary
============================================
  Run Name                 : ${workflow.runName}
  Launch Directory         : ${workflow.launchDir}
  Work Directory           : ${workflow.workDir}
  Project Directory        : ${projectDir}
  Python Executable        : ${params.python_executable}
  --- GCS Paths ---
  Output Base              : ${params.output_dir_base}/${workflow.runName}
  Reusable Temp Base       : ${params.gcs_temp_dir_base}
  Hail Run Temp Base       : ${params.gcs_hail_temp_dir_base}
  --- Input Data ---
  AoU WGS VDS Path         : ${params.aou_wgs_vds_path}
  AoU Flagged Samples Path : ${params.aou_flagged_samples_gcs_path}
  Models Source CSV        : ${params.models_csv}
  --- Phenotype ---
  Target Phenotype         : ${params.target_phenotype_name}
  Phenotype Concept IDs    : ${params.phenotype_concept_ids}
  --- VDS Preparation ---
  Enable Downsampling      : ${params.enable_downsampling_for_vds_generation}
  Cases (Downsample)       : ${params.n_cases_downsample}
  Controls (Downsample)    : ${params.n_controls_downsample}
  Downsampling Seed        : ${params.downsampling_random_state}
  --- Execution ---
  Executor                 : ${workflow.executor}
  Environment Management   : Relies on pre-configured system Python
============================================
"""

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
        if (!row.id || !row.url || !row.phenotype) {
            error "Invalid row in models CSV: ${row}. Requires 'id', 'url', and 'phenotype' columns."
        }
        [ id: row.id.trim(), url: row.url.trim(), phenotype_label: row.phenotype.trim() ]
    }
    .set { model_info_ch }

model_info_ch.count().subscribe { count -> log.info "Found ${count} models to process from ${params.models_csv}" }

// ----------------------------------------------------------------------------
// Processes
// ----------------------------------------------------------------------------

process fetch_phenotype_cases {
    tag "Fetch Pheno: ${params.target_phenotype_name}"

    publishDir "${params.output_dir_base}/${workflow.runName}/logs/01_fetch_pheno", mode:'copy', overwrite:true, pattern: '*.{log,txt}'

    output:
    path "phenotype_cases_gcs_path.txt", emit: pheno_cases_gcs_path_file

    script:
    // Define the GCS path where the phenotype CSV will be written by the python script
    // This path is within the reusable gcs_temp_dir_base
    def gcs_pheno_out_path = "${params.gcs_temp_dir_base}/phenotype_data/${params.target_phenotype_name.replace(' ','_')}_cases.csv"
    """
    echo "INFO: Fetching phenotype data for '${params.target_phenotype_name}'..."
    "${params.python_executable}" "${projectDir}/src/fetch_phenotypes.py" \\
        --phenotype_name                "${params.target_phenotype_name}" \\
        --phenotype_concept_ids         "${params.phenotype_concept_ids.join(',')}" \\
        --workspace_cdr                 "${params.aou_workspace_cdr}" \\
        --output_phenotype_csv_gcs_path "${gcs_pheno_out_path}" \\
        --google_billing_project        "${params.google_billing_project}" \\
        > fetch_phenotypes.log 2>&1

    echo "${gcs_pheno_out_path}" > phenotype_cases_gcs_path.txt
    echo "INFO: Phenotype GCS path recorded in phenotype_cases_gcs_path.txt: ${gcs_pheno_out_path}"
    """
}

process prepare_base_vds {
    tag "Prepare Base VDS"

    publishDir "${params.output_dir_base}/${workflow.runName}/logs/02_prepare_vds", mode:'copy', overwrite:true, pattern: '*.{log,txt}'

    input:
    path phenotype_cases_gcs_path_file // From fetch_phenotype_cases

    output:
    path "base_vds_gcs_path.txt",    emit: base_vds_path_file
    path "wgs_ehr_ids_gcs_path.txt", emit: wgs_ids_path_file

    script:
    // Define GCS paths for the outputs of this script
    def gcs_vds_out_path = "${params.gcs_temp_dir_base}/base_cohort_wgs_ehr_unrelated.vds" // Reusable checkpoint
    def gcs_ids_out_path = "${params.gcs_temp_dir_base}/cohort_definitions/people_with_WGS_EHR_ids_for_vds_filtering.csv" // Reusable checkpoint

    // Prepare Spark configurations to be passed as a JSON string
    def spark_confs_json_for_script = groovy.json.JsonOutput.toJson(params.spark_conf_list)

    """
    echo "INFO: Preparing Base VDS..."
    pheno_input_gcs_path=\$(cat "${phenotype_cases_gcs_path_file}")
    echo "INFO: Using Phenotype Cases CSV from GCS path: \${pheno_input_gcs_path}"

    "${params.python_executable}" "${projectDir}/src/prepare_base_vds.py" \\
        --workspace_cdr                 "${params.aou_workspace_cdr}" \\
        --run_timestamp                 "${workflow.runName}" \\
        --gcs_temp_dir                  "${params.gcs_temp_dir_base}" \\
        --gcs_hail_temp_dir             "${params.gcs_hail_temp_dir_base}" \\
        --wgs_vds_path                  "${params.aou_wgs_vds_path}" \\
        --flagged_samples_gcs_path      "${params.aou_flagged_samples_gcs_path}" \\
        --base_cohort_vds_path_out      "${gcs_vds_out_path}" \\
        --wgs_ehr_ids_gcs_path_out      "${gcs_ids_out_path}" \\
        --target_phenotype_name         "${params.target_phenotype_name}" \\
        --phenotype_cases_gcs_path_input "\${pheno_input_gcs_path}" \\
        ${params.enable_downsampling_for_vds_generation ? '--enable_downsampling_for_vds' : ''} \\
        --n_cases_downsample            ${params.n_cases_downsample} \\
        --n_controls_downsample         ${params.n_controls_downsample} \\
        --downsampling_random_state     ${params.downsampling_random_state} \\
        --google_billing_project        "${params.google_billing_project}" \\
        --spark_configurations_json     '${spark_confs_json_for_script}' \\
        > prepare_base_vds.log 2>&1

    echo "${gcs_vds_out_path}" > base_vds_gcs_path.txt
    echo "${gcs_ids_out_path}" > wgs_ehr_ids_gcs_path.txt
    echo "INFO: Output Base VDS GCS path recorded: ${gcs_vds_out_path}"
    echo "INFO: Output WGS+EHR IDs GCS path recorded: ${gcs_ids_out_path}"
    """
}

process process_prs_model {
    tag "Process PRS: ${model.id}"

    publishDir "${params.output_dir_base}/${workflow.runName}/scores/${model.id}/logs_calculation", mode:'copy', overwrite:true, pattern: '*.{log,txt}'

    input:
    tuple val(model), path(base_vds_gcs_path_file)

    output:
    tuple val(model), path("final_score_gcs_path.txt"), emit: model_and_score_file

    script:
    def base_vds_path = "\$(cat "${base_vds_gcs_path_file}")"
    // Outputs for this specific run
    def run_specific_output_dir = "${params.output_dir_base}/${workflow.runName}"
    def final_hail_out_path = "${run_specific_output_dir}/scores/${model.id}/hail_table/${model.id}_scores.ht"
    def final_csv_out_path  = "${run_specific_output_dir}/scores/${model.id}/score_csv/${model.id}_scores.csv"

    // Prepare Spark configurations
    def spark_confs_json_for_script = groovy.json.JsonOutput.toJson(params.spark_conf_list)

    """
    echo "INFO: Processing PRS model ${model.id}..."
    echo "INFO: Using Base VDS from GCS path: ${base_vds_path}"

    "${params.python_executable}" "${projectDir}/src/process_prs_model.py" \\
        --prs_id                            "${model.id}" \\
        --prs_url                           "${model.url}" \\
        --base_cohort_vds_path            "${base_vds_path}" \\
        --gcs_temp_dir                      "${params.gcs_temp_dir_base}" \\
        --gcs_hail_temp_dir                 "${params.gcs_hail_temp_dir_base}" \\
        --run_timestamp                     "${workflow.runName}" \\
        --output_final_hail_table_gcs_path  "${final_hail_out_path}" \\
        --output_final_score_csv_gcs_path   "${final_csv_out_path}" \\
        --google_billing_project            "${params.google_billing_project}" \\
        --spark_configurations_json         '${spark_confs_json_for_script}' \\
        > process_prs_model_${model.id}.log 2>&1

    echo "${final_csv_out_path}" > final_score_gcs_path.txt
    echo "INFO: Final score CSV GCS path for ${model.id} recorded: ${final_csv_out_path}"
    """
}

process analyze_one_model_results {
    tag "Analyze: ${model.id} for ${params.target_phenotype_name}"

    publishDir "${params.output_dir_base}/${workflow.runName}/scores/${model.id}/analysis_results", mode:'copy', overwrite:true, pattern: '*.{log,png,txt}'

    input:
    tuple val(model), path(score_gcs_path_file), path(wgs_ehr_ids_gcs_path_file), path(phenotype_cases_gcs_path_file)

    output:
    path "${model.id}_analysis_summary.txt", emit: summary_file // This is a local file in workDir

    script:
    def score_csv_path = "\$(cat "${score_gcs_path_file}")"
    def wgs_ehr_ids_path = "\$(cat "${wgs_ehr_ids_gcs_path_file}")"
    def pheno_cases_path = "\$(cat "${phenotype_cases_gcs_path_file}")"
    // GCS base path for this run, where plots etc. will be written by the python script
    def run_specific_gcs_output_dir = "${params.output_dir_base}/${workflow.runName}"

    """
    echo "INFO: Analyzing results for PRS model ${model.id}..."
    echo "INFO: Score CSV GCS path: ${score_csv_path}"
    echo "INFO: WGS+EHR IDs GCS path: ${wgs_ehr_ids_path}"
    echo "INFO: Phenotype Cases GCS path: ${pheno_cases_path}"

    "${params.python_executable}" "${projectDir}/src/analyze_prs_results.py" \\
        --prs_id                        "${model.id}" \\
        --prs_phenotype_label           "${model.phenotype_label}" \\
        --score_csv_gcs_path            "${score_csv_path}" \\
        --wgs_ehr_ids_gcs_path          "${wgs_ehr_ids_path}" \\
        --phenotype_cases_csv_gcs_path  "${pheno_cases_path}" \\
        --phenotype_name                "${params.target_phenotype_name}" \\
        --gcs_base_output_dir_run       "${run_specific_gcs_output_dir}" \\
        --run_timestamp                 "${workflow.runName}" \\
        --google_billing_project        "${params.google_billing_project}" \\
        --output_summary_file_name      "${model.id}_analysis_summary.txt" \\
        > analyze_results_${model.id}.log 2>&1

    echo "INFO: Analysis for ${model.id} finished. Summary file: ${model.id}_analysis_summary.txt"
    """
}

// ----------------------------------------------------------------------------
// Workflow Execution
// ----------------------------------------------------------------------------

workflow {
    // Stage 1: Fetch phenotype data (single output)
    fetch_phenotype_cases()
    pheno_path_ch = fetch_phenotype_cases.out.pheno_cases_gcs_path_file

    // Stage 2: Prepare VDS (takes single pheno path, produces two separate single paths)
    prepare_base_vds(pheno_path_ch)
    base_vds_path_ch = prepare_base_vds.out.base_vds_path_file
    wgs_ids_path_ch  = prepare_base_vds.out.wgs_ids_path_file

    // Stage 3: Combine model info (multi-item) with base VDS path (single item)
    process_prs_inputs_ch = model_info_ch.combine(base_vds_path_ch)

    // Stage 4: Run PRS model processing using the combined channel
    process_prs_model(process_prs_inputs_ch)
    processed_scores_ch = process_prs_model.out.model_and_score_file

    // Stage 5: Prepare inputs for analyze_one_model_results
    analysis_inputs_intermediate_ch = processed_scores_ch.combine(wgs_ids_path_ch)
    final_analysis_inputs_ch        = analysis_inputs_intermediate_ch.combine(pheno_path_ch)

    // Stage 6: Run analysis using the final combined channel
    analyze_one_model_results(final_analysis_inputs_ch)
    analysis_summary_files_ch = analyze_one_model_results.out.summary_file

    // Stage 7: Collect the summary files
    analysis_summary_files_ch.collectFile(
        name: 'all_models_analysis_summaries.txt', // Changed name for clarity
        storeDir: "${params.output_dir_base}/${workflow.runName}/analysis_summary_collection",
        overwrite: true
    )
}

// --- Workflow Handlers ---

workflow.onComplete {
    def summary_collection_file = "${params.output_dir_base}/${workflow.runName}/analysis_summary_collection/all_models_analysis_summaries.txt"
    // Using direct param access for paths, error checks in config should ensure they are set.
    log.info """
    ============================================
      PGS Analysis Workflow Completed
    ============================================
      Run Name         : ${workflow.runName}
      Status           : ${workflow.success ? 'Success' : 'Failed'}
      Exit Status      : ${workflow.exitStatus}
      Duration         : ${workflow.duration}
      Work Dir         : ${workflow.workDir}
      Launch Dir       : ${projectDir}
      Run Output Base  : ${params.output_dir_base}/${workflow.runName}
      Reusable Temp    : ${params.gcs_temp_dir_base}
      Collected analysis summaries (if successful): ${summary_collection_file}
    ============================================
    """.stripIndent()
}

workflow.onError {
    log.error """
    ============================================
      PGS Analysis Workflow Error!
    ============================================
      Run Name         : ${workflow.runName}
      Status           : Failed
      Exit Status      : ${workflow.exitStatus}
      Error report     : ${workflow.errorReport ?: 'N/A'}
      Trace file       : ${params.output_dir_base}/${workflow.runName}/reports/trace.txt
      Check work directory for task logs: ${workflow.workDir}
      Check Nextflow log (${workflow.logFile ?: '.nextflow.log'}) for more details.
    ============================================
    """.stripIndent()
}
