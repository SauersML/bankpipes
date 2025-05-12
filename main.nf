#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------------------------------------------------------------------
// P G S   A N A L Y S I S   N E X T F L O W   P I P E L I N E
// ----------------------------------------------------------------------------

params.wgs_vds_path               = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/vds/hail.vds'
params.flagged_samples_gcs_path   = 'gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv'
params.phenotype_concept_ids       = [432510, 374919]
params.target_phenotype_name       = "Ischemic Stroke"

params.enable_downsampling_for_vds_generation = true
params.n_cases_downsample                     = 500
params.n_controls_downsample                  = 500
params.downsampling_random_state              = 2025

params.gcs_temp_dir_base            = "${System.getenv('WORKSPACE_BUCKET')}/prs_analysis_temp_nextflow"
params.gcs_hail_temp_dir_base       = "${System.getenv('WORKSPACE_BUCKET')}/hail_temp_nextflow"
params.output_dir_base              = "${System.getenv('WORKSPACE_BUCKET')}/prs_analysis_runs_nextflow"

params.models_list = [
    [ id: 'PGS002724', url: 'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS002724/ScoringFiles/Harmonized/PGS002724_hmPOS_GRCh38.txt.gz' ],
    [ id: 'PGS000039', url: 'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS000039/ScoringFiles/Harmonized/PGS000039_hmPOS_GRCh38.txt.gz' ]
]

workflow.onStart {
    log.info """
    ============================================
      P G S   A N A L Y S I S   WORKFLOW
    ============================================
      WGS VDS Path         : ${params.wgs_vds_path}
      Flagged Samples Path : ${params.flagged_samples_gcs_path}
      Phenotype            : ${params.target_phenotype_name}
      Models to run        : ${params.models_list.size()}
      Temp Dir Base        : ${params.gcs_temp_dir_base}
      Output Dir Base      : ${params.output_dir_base}
      Run Name             : ${workflow.runName}
      LaunchDir            : ${workflow.launchDir}/results
      WorkDir              : ${workflow.workDir}
    ============================================
    """.stripIndent()
}

// ----------------------------------------------------------------------------
// Channels
// ----------------------------------------------------------------------------

model_channel = Channel.fromList(params.models_list)

// ----------------------------------------------------------------------------
// Processes
// ----------------------------------------------------------------------------

// 1) Fetch phenotype cases (run once)
process fetch_phenotype_cases {
    tag "${params.target_phenotype_name}"

    publishDir "${workflow.launchDir}/results/fetch_pheno", mode:'copy', overwrite:true

    output:
    path "phenotype_cases_gcs_path.txt"

    script:
    // output path in GCS
    def gcs_out = "${params.gcs_temp_dir_base}/phenotype_cases/${params.target_phenotype_name.replace(' ','_')}_cases.csv"
    """
    python3 ${projectDir}/src/fetch_phenotypes.py \\
        --phenotype_name    "${params.target_phenotype_name}" \\
        --phenotype_concept_ids "${params.phenotype_concept_ids.join(',')}" \\
        --workspace_cdr     "${System.getenv('WORKSPACE_CDR')}" \\
        --output_phenotype_csv_gcs_path "${gcs_out}" \\
        --project_bucket    "${System.getenv('WORKSPACE_BUCKET')}"

    // record the GCS path for downstream
    echo "${gcs_out}" > phenotype_cases_gcs_path.txt
    """
}

// 2) Prepare Base VDS (run once)
process prepare_base_vds {
    publishDir "${workflow.launchDir}/results/base_vds", mode:'copy', overwrite:true

    output:
    path "base_vds_path.txt"
    path "wgs_ehr_ids_gcs_path.txt"

    script:
    def vds_out   = "${params.gcs_temp_dir_base}/base_cohort_wgs_ehr_unrelated.vds"
    def ids_out   = "${params.gcs_temp_dir_base}/people_with_WGS_EHR_ids.csv"
    """
    python3 ${projectDir}/src/prepare_base_vds.py \\
        --project_bucket             "${System.getenv('WORKSPACE_BUCKET')}" \\
        --workspace_cdr              "${System.getenv('WORKSPACE_CDR')}" \\
        --run_timestamp              "${workflow.runName}" \\
        --gcs_temp_dir               "${params.gcs_temp_dir_base}" \\
        --gcs_hail_temp_dir          "${params.gcs_hail_temp_dir_base}" \\
        --wgs_vds_path               "${params.wgs_vds_path}" \\
        --flagged_samples_gcs_path   "${params.flagged_samples_gcs_path}" \\
        --enable_downsampling_for_vds ${params.enable_downsampling_for_vds_generation} \\
        --n_cases_downsample         ${params.n_cases_downsample} \\
        --n_controls_downsample      ${params.n_controls_downsample} \\
        --downsampling_random_state  ${params.downsampling_random_state} \\
        --target_phenotype_name      "${params.target_phenotype_name}" \\
        --phenotype_concept_ids      "${params.phenotype_concept_ids.join(',')}" \\
        --base_cohort_vds_path_out   "${vds_out}" \\
        --wgs_ehr_ids_gcs_path_out   "${ids_out}"

    echo "${vds_out}" > base_vds_path.txt
    echo "${ids_out}" > wgs_ehr_ids_gcs_path.txt
    """
}

// 3) Process each PRS model in parallel
process process_prs_model {
    tag "${model.id}"

    input:
    tuple val(model), path("base_vds_path.txt")

    output:
    tuple val(model.id), path("final_score_gcs_path.txt")

    script:
    def vds_path = file("base_vds_path.txt").text.trim()
    def hail_tmp = params.gcs_hail_temp_dir_base
    def temp_dir = params.gcs_temp_dir_base
    def out_base = "${params.output_dir_base}/${workflow.runName}"

    def hail_out = "${out_base}/scores/${model.id}/hail_table"
    def csv_out  = "${out_base}/scores/${model.id}/score/${model.id}_scores.csv"

    """
    python3 ${projectDir}/src/process_prs_model.py \\
        --prs_id                            "${model.id}" \\
        --prs_url                           "${model.url}" \\
        --base_cohort_vds_path              "${vds_path}" \\
        --gcs_temp_dir                      "${temp_dir}" \\
        --gcs_hail_temp_dir                 "${hail_tmp}" \\
        --run_timestamp                     "${workflow.runName}" \\
        --output_final_hail_table_gcs_path  "${hail_out}" \\
        --output_final_score_csv_gcs_path   "${csv_out}"

    echo "${csv_out}" > final_score_gcs_path.txt
    """
}

// 4) Analyze all results (run once)
process analyze_all_results {
    tag "analysis"

    input:
    path score_txts   // collection of final_score_gcs_path.txt files
    path base_ids_txt // wgs_ehr_ids_gcs_path.txt from prepare_base_vds
    path pheno_txt    // phenotype_cases_gcs_path.txt from fetch_phenotype_cases

    publishDir "${workflow.launchDir}/results/analysis", mode:'copy', overwrite:true

    output:
    path "analysis_summary.log"

    script:
    """
    # build manifest
    > score_manifest.txt
    for f in ${score_txts.join(' ')}; do
      cat \$f >> score_manifest.txt
      echo >> score_manifest.txt
    done

    FINAL_WGS_IDS=\$(cat ${base_ids_txt})
    PHENO_CASES=\$(cat ${pheno_txt})

    python3 ${projectDir}/src/analyze_all_results.py \\
        --score_manifest               score_manifest.txt \\
        --wgs_ehr_ids_gcs_path         "\${FINAL_WGS_IDS}" \\
        --phenotype_cases_csv_gcs_path "\${PHENO_CASES}" \\
        --phenotype_name               "${params.target_phenotype_name}" \\
        --gcs_output_dir_base          "${params.output_dir_base}/${workflow.runName}" \\
        --run_timestamp                "${workflow.runName}" \\
        --project_bucket               "${System.getenv('WORKSPACE_BUCKET')}" \\
        --output_summary_log           analysis_summary.log
    """
}

// ----------------------------------------------------------------------------
// Workflow wiring
// ----------------------------------------------------------------------------

workflow {
    // fetch phenotype cases
    pheno_ch = fetch_phenotype_cases()

    // prepare base VDS + WGS/EHR IDs
    base_ch  = prepare_base_vds()

    // process each PRS model
    prs_ch   = process_prs_model( model_channel.combine(base_ch.map{ it[0] }) )
