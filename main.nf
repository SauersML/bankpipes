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

Channel
    .fromPath(params.models_csv)
    .ifEmpty { error "Models CSV file not found or empty: ${params.models_csv}" }
    .splitCsv(header:true, sep:',')
    .map { row ->
        if (!row.id || !row.url || !row.phenotype) { // Added phenotype for plotting label
            error "Invalid row in models CSV: ${row}. Requires 'id', 'url', and 'phenotype' columns."
        }
        [ id: row.id.trim(), url: row.url.trim(), phenotype_label: row.phenotype.trim() ]
    }
    .set { model_info_ch } // Channel emits: [ id:'PGS...', url:'...', phenotype_label:'Ischemic Stroke' ]

model_info_ch.count().subscribe { count -> log.info "Found ${count} models to process from ${params.models_csv}" }


// ----------------------------------------------------------------------------
// Processes
// ----------------------------------------------------------------------------

process fetch_phenotype_cases {
    tag "Fetch Pheno: ${params.target_phenotype_name}"
    label 'python_bq_env' // Specific, lighter environment

    publishDir "${params.output_dir_base}/${workflow.runName}/logs/01_fetch_pheno", mode:'copy', overwrite:true, pattern: '*.{log,txt}'

    output:
    path "phenotype_cases_gcs_path.txt", emit: pheno_cases_gcs_path_file

    script:
    def gcs_pheno_out_path = "${params.gcs_temp_dir_base}/phenotype_data/${params.target_phenotype_name.replace(' ','_')}_cases.csv"
    """
    echo "INFO: Fetching phenotype data for '${params.target_phenotype_name}'..."
    python3 ${projectDir}/src/fetch_phenotypes.py \\
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
    tuple path("base_vds_gcs_path.txt"), path("wgs_ehr_ids_gcs_path.txt"), emit: base_vds_and_ids_files

    script:
    def gcs_vds_out_path = "${params.gcs_temp_dir_base}/base_cohort_wgs_ehr_unrelated.vds"
    def gcs_ids_out_path = "${params.gcs_temp_dir_base}/cohort_definitions/people_with_WGS_EHR_ids_for_vds_filtering.csv" // More specific name
    def pheno_input_gcs_path = "\$(cat ${phenotype_cases_gcs_path_file})"

    """
    echo "INFO: Preparing Base VDS..."
    echo "INFO: Using Phenotype Cases CSV from GCS path: ${pheno_input_gcs_path}"

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
        --target_phenotype_name             "${params.target_phenotype_name}" \\
        --phenotype_concept_ids             "${params.phenotype_concept_ids.join(',')}" \\
        ${params.enable_downsampling_for_vds_generation ? '--enable_downsampling_for_vds' : ''} \\
        --n_cases_downsample                ${params.n_cases_downsample} \\
        --n_controls_downsample             ${params.n_controls_downsample} \\
        --downsampling_random_state         ${params.downsampling_random_state} \\
        > prepare_base_vds.log 2>&1

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
    tuple val(model), path(base_vds_gcs_path_file) // model = [id, url, phenotype_label]

    output:
    tuple val(model), path("final_score_gcs_path.txt"), emit: model_and_score_file // Pass model info through for analysis

    script:
    def base_vds_path = "\$(cat ${base_vds_gcs_path_file})"
    def run_specific_output_dir = "${params.output_dir_base}/${workflow.runName}" // For final outputs of this run
    def final_hail_out_path = "${run_specific_output_dir}/scores/${model.id}/hail_table/${model.id}_scores.ht"
    def final_csv_out_path  = "${run_specific_output_dir}/scores/${model.id}/score_csv/${model.id}_scores.csv"

    """
    echo "INFO: Processing PRS model ${model.id}..."
    echo "INFO: Using Base VDS from GCS path: ${base_vds_path}"

    python3 ${projectDir}/src/process_prs_model.py \\
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
    label 'pyhail_env' // Needs pandas, sklearn, viz libs (same as Hail for simplicity now)

    publishDir "${params.output_dir_base}/${workflow.runName}/scores/${model.id}/analysis_results", mode:'copy', overwrite:true, pattern: '*.{log,png,txt}'

    input:
    tuple val(model), path(score_gcs_path_file), path(wgs_ehr_ids_gcs_path_file), path(phenotype_cases_gcs_path_file)
    // model = [id, url, phenotype_label]

    output:
    path "${model.id}_analysis_summary.txt", emit: summary_file // A summary output for this model

    script:
    def score_csv_path = "\$(cat ${score_gcs_path_file})"
    def wgs_ehr_ids_path = "\$(cat ${wgs_ehr_ids_gcs_path_file})"
    def pheno_cases_path = "\$(cat ${phenotype_cases_gcs_path_file})"
    def run_specific_output_dir = "${params.output_dir_base}/${workflow.runName}" // For plots etc.

    """
    echo "INFO: Analyzing results for PRS model ${model.id}..."
    echo "INFO: Score CSV GCS path: ${score_csv_path}"
    echo "INFO: WGS+EHR IDs GCS path: ${wgs_ehr_ids_path}"
    echo "INFO: Phenotype Cases GCS path: ${pheno_cases_path}"

    python3 ${projectDir}/src/analyze_prs_results.py \\
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
// Workflow Execution
// ----------------------------------------------------------------------------

workflow {
    pheno_path_ch = fetch_phenotype_cases() // Emits phenotype_cases_gcs_path.txt

    vds_and_ids_ch = prepare_base_vds(pheno_path_ch)
        // Emits: tuple path(base_vds_gcs_path.txt), path(wgs_ehr_ids_gcs_path.txt)

    // Prepare inputs for process_prs_model
    // It needs: tuple val(model_map), path(base_vds_gcs_path_file)
    process_prs_inputs_ch = model_info_ch.combine(vds_and_ids_ch.map { it[0] })
        // Emits: [ [id, url, phenotype_label], base_vds_gcs_path.txt ]

    processed_scores_ch = process_prs_model(process_prs_inputs_ch)
        // Emits: tuple val(model_map), path(final_score_gcs_path.txt)


    // Prepare inputs for analyze_one_model_results
    // It needs: tuple val(model_map), path(score_gcs_path_file), path(wgs_ehr_ids_gcs_path_file), path(phenotype_cases_gcs_path_file)
    // We have:
    //   1. processed_scores_ch: [ [id, url, phenotype_label], final_score_gcs_path.txt ]
    //   2. vds_and_ids_ch.map { it[1] }: wgs_ehr_ids_gcs_path.txt (value channel)
    //   3. pheno_path_ch: phenotype_cases_gcs_path.txt (value channel)

    // Combine processed_scores_ch with the wgs_ehr_ids_gcs_path.txt.
    // Since vds_and_ids_ch is a single emission, its map {it[1]} will also be.
    // We want to pair each model's score with this single WGS ID path.
    analysis_inputs_intermediate_ch = processed_scores_ch.combine(vds_and_ids_ch.map { it[1] })
        // Emits: [ [id, url, phenotype_label], final_score_gcs_path.txt, wgs_ehr_ids_gcs_path.txt ]

    // Now combine with the single phenotype_cases_gcs_path.txt
    final_analysis_inputs_ch = analysis_inputs_intermediate_ch.combine(pheno_path_ch)
        // Emits: [ [id, url, phenotype_label], final_score_gcs_path.txt, wgs_ehr_ids_gcs_path.txt, phenotype_cases_gcs_path.txt ]

    analysis_summary_files_ch = analyze_one_model_results(final_analysis_inputs_ch)

    analysis_summary_files_ch.collectFile(name:'all_analysis_summaries.txt', storeDir:"${params.output_dir_base}/${workflow.runName}/analysis_summary_collection")
}

workflow.onComplete {
    def summary_file = "${params.output_dir_base}/${workflow.runName}/analysis_summary_collection/all_analysis_summaries.txt"
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
      Run Output  : ${params.output_dir_base}/${workflow.runName}
      Temp Files  : ${params.gcs_temp_dir_base}
      Collected analysis summaries (if successful): ${summary_file} (local path may vary if not using local executor for final collection)
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
      Trace file  : ${params.output_dir_base}/${workflow.runName}/reports/trace.txt (if generated)
      Check work directory for task logs: ${workflow.workDir}
      Check Nextflow log for more details.
    ============================================
    """.stripIndent()
}
