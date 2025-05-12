#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// --- Log parameters ---
log.info """
         P G S  A N A L Y S I S  W O R K F L O W
         =======================================
         WGS VDS Path          : ${params.wgs_vds_path}
         Flagged Samples Path  : ${params.flagged_samples_gcs_path}
         Output Directory Base : ${params.output_dir_base}
         Temp Directory Base   : ${params.gcs_temp_dir_base}
         Phenotype Name        : ${params.target_phenotype_name}
         Models to run         : ${params.models_list.size()}
         ---
         Run Name              : ${workflow.runName}
         Run Output Dir        : ${workflow.launchDir}/results
         Run Work Dir          : ${workflow.workDir}
         """
         .stripIndent()

// --- Define Channels ---
// Channel for PRS models
model_channel = Channel.fromList(params.models_list)

// Path for the base VDS
BASE_COHORT_VDS_GCS_PATH = "${params.gcs_temp_dir_base}/base_cohort_wgs_ehr_unrelated.vds"

// --- Define Processes ---
process prepare_base_vds {
    publishDir "${workflow.launchDir}/results/base_vds_log", mode: 'copy', overwrite: true, pattern: "*.log"

    output:
    val base_vds_path

    script:
    """
    python3 ${projectDir}/scripts/prepare_base_vds.py \
        --base_vds_gcs_path "${BASE_COHORT_VDS_GCS_PATH}" \
        --wgs_vds_path "${params.wgs_vds_path}" \
        --flagged_samples_path "${params.flagged_samples_gcs_path}" \
        --phenotype_name "${params.target_phenotype_name}" \
        --phenotype_ids_str "${params.phenotype_concept_ids.join(',')}" \
        --enable_downsampling ${params.enable_downsampling_for_vds_generation} \
        --n_cases ${params.n_cases_downsample} \
        --n_controls ${params.n_controls_downsample} \
        --random_state ${params.downsampling_random_state} \
        --project_bucket "${System.getenv('WORKSPACE_BUCKET')}" \
        --gcs_temp_dir "${params.gcs_temp_dir_base}" \
        --hail_temp_dir "${params.gcs_hail_temp_dir_base}" \
        --output_wgs_ehr_ids_csv "${params.gcs_temp_dir_base}/people_with_WGS_EHR_ids.csv" \
        
    # Output the path for the next step
    echo "${BASE_COHORT_VDS_GCS_PATH}" > base_vds_path.txt
    """
    afterScript "echo ${BASE_COHORT_VDS_GCS_PATH} > base_vds_path_output.txt"
    output:
        path("base_vds_path_output.txt")
        path("${params.gcs_temp_dir_base}/people_with_WGS_EHR_ids.csv")

}


process process_prs_model {
    tag "${model.id}"

    input:
    tuple val(model), path(base_vds_path_file) // model is a map [id: 'PGS...', url: '...']

    output:
    tuple val(model.id), path("final_score_gcs_path.txt") // File containing the GCS path of the score CSV

    script:
    base_vds_gcs_path = base_vds_path_file.getText('UTF-8').trim()

    """
    python3 ${projectDir}/scripts/process_prs_model.py \
        --prs_id "${model.id}" \
        --prs_url "${model.url}" \
        --base_vds_path "${base_vds_gcs_path}" \
        --project_bucket "${System.getenv('WORKSPACE_BUCKET')}" \
        --gcs_temp_dir_base "${params.gcs_temp_dir_base}" \
        --gcs_output_dir_base "${params.output_dir_base}/${workflow.runName}" \
        --hail_temp_dir_base "${params.gcs_hail_temp_dir_base}" \
        --output_path_file "final_score_gcs_path.txt"

    """
}

process analyze_results {

    input:
    path all_score_files_info // this will be a list of files, each containing a GCS path
    path wgs_ehr_ids_csv_from_prep // from prepare_base_vds

    output:
    path "analysis_summary.log"
    script:
    """
    # Create a manifest of score CSV GCS paths
    > score_manifest.txt
    for f in $all_score_files_info; do
        cat \$f >> score_manifest.txt
        echo "" >> score_manifest.txt # new line
    done
    
    python3 ${projectDir}/scripts/analyze_results.py \
        --score_manifest "score_manifest.txt" \
        --phenotype_name "${params.target_phenotype_name}" \
        --phenotype_ids_str "${params.phenotype_concept_ids.join(',')}" \
        --wgs_ehr_ids_csv "${wgs_ehr_ids_csv_from_prep}" \
        --project_bucket "${System.getenv('WORKSPACE_BUCKET')}" \
        --gcs_output_dir_base "${params.output_dir_base}/${workflow.runName}" \
        --output_summary_log "analysis_summary.log"
    """
}


// --- Workflow Definition ---
workflow {
    // 1. Prepare Base VDS (runs once)
    //    Its output includes the path to the VDS and the WGS+EHR IDs CSV
    //    and potentially the phenotype_cases_df.pkl
    ch_base_data = prepare_base_vds()
    
    base_vds_path_output_file = ch_base_data.map { it[0] } // base_vds_path_output.txt
    wgs_ehr_ids_gcs_file = ch_base_data.map { it[1] }     // people_with_WGS_EHR_ids.csv

    // 2. Process each PRS model in parallel
    //    It takes the base_vds_path and the model info
    ch_prs_scores_info = process_prs_model(model_channel.combine(base_vds_path_output_file))
    // ch_prs_scores_info will emit tuples: [model_id, path_to_txt_file_with_gcs_score_path]

    // 3. Collect all score file paths and run analysis
    analyze_results(
        ch_prs_scores_info.map{ it[1] }.collect(), // Collect all the .txt files
        wgs_ehr_ids_gcs_file.first() // .first() because prepare_base_vds runs once
    )
}

workflow.onComplete {
    log.info ( workflow.success ? """
        Pipeline execution summary
        ---------------------------
        Completed at : ${workflow.complete}
        Duration     : ${workflow.duration}
        Success      : ${workflow.success}
        Work directory : ${workflow.workDir}
        Output results : ${workflow.launchDir}/results (and GCS paths logged by scripts)
        """ : """
        Failed: ${workflow.errorReport}
        Process: ${workflow.process}
        """ )
}
