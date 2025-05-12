import argparse
import datetime
import os
import sys
import time
import pandas as pd
import hail as hl
from common_utils import (
    init_hail, gcs_path_exists, hail_path_exists, delete_gcs_path,
    cache_result, get_gcs_fs, get_cache_dir
)

def parse_args():
    parser = argparse.ArgumentParser(description="Prepare Base Cohort VDS for PRS Analysis.")
    parser.add_argument("--project_bucket", required=True, help="GCS project bucket (gs://your-bucket).")
    parser.add_argument("--workspace_cdr", required=True, help="Workspace CDR (e.g., fc-aou-cdr-prod-ct.CYYYYQQRR).")
    parser.add_argument("--run_timestamp", required=True, help="Run timestamp (YYYYMMDD_HHMMSS).")
    parser.add_argument("--gcs_temp_dir", required=True, help="GCS temporary directory for stable checkpoints.")
    parser.add_argument("--gcs_hail_temp_dir", required=True, help="GCS temporary directory for Hail operations.")
    parser.add_argument("--wgs_vds_path", required=True, help="GCS path to the full WGS VDS.")
    parser.add_argument("--flagged_samples_gcs_path", required=True, help="GCS path to the flagged samples TSV.")
    parser.add_argument("--base_cohort_vds_path_out", required=True, help="GCS output path for the prepared base cohort VDS.")
    parser.add_argument("--wgs_ehr_ids_gcs_path_out", required=True, help="GCS output path for the WGS+EHR sample IDs CSV.")

    parser.add_argument("--enable_downsampling_for_vds", action='store_true', help="Enable downsampling for VDS generation.")
    parser.add_argument("--n_cases_downsample", type=int, default=500, help="Number of cases for downsampling.")
    parser.add_argument("--n_controls_downsample", type=int, default=500, help="Number of controls for downsampling.")
    parser.add_argument("--downsampling_random_state", type=int, default=2025, help="Random state for downsampling.")
    parser.add_argument("--target_phenotype_name", required=True, help="Target phenotype name for downsampling logic.")
    parser.add_argument("--phenotype_concept_ids", required=True, help="Comma-separated OMOP concept IDs for the phenotype.")
    return parser.parse_args()

# --- Phenotype Data Fetching (copied from original, adapted for script use) ---
@cache_result("phenotype_df") # Name prefix for caching
def get_phenotype_data_for_vds_prep(phenotype_name, concept_ids_str, cdr_env_var_value, prs_id=None): # prs_id not used here but cache_result might expect it
    """Fetches case status based on OMOP concept IDs."""
    # get local cache dir exists for this script execution
    get_cache_dir()
    print(f"Retrieving {phenotype_name} phenotype (Cases) from BigQuery...")

    if not cdr_env_var_value:
        print(f"FATAL ERROR: Workspace CDR value not provided.")
        sys.exit(1)

    try:
        concept_ids = [int(cid.strip()) for cid in concept_ids_str.split(',')]
    except ValueError:
        print(f"FATAL ERROR: Invalid phenotype_concept_ids format: {concept_ids_str}. Must be comma-separated integers.")
        sys.exit(1)

    if not concept_ids:
        print("FATAL ERROR: No phenotype concept IDs provided.")
        sys.exit(1)
    concept_id_string = ', '.join(map(str, concept_ids))

    phenotype_query = f"""
    SELECT DISTINCT co.person_id
    FROM `{cdr_env_var_value}.condition_occurrence` co
    WHERE co.condition_concept_id IN ({concept_id_string})
    """
    print("Executing query:\n", phenotype_query)
    try:
        # Assuming BigQuery client is configured in the environment Nextflow runs this.
        # For bare metal, 'gcloud auth application-default login' or service account needed.
        # The project for billing needs to be set up for pandas-gbq.
        # It might try to infer project from environment, or you can set it via:
        # pd.options.io.bigquery.project = 'your-billing-project-id'
        # Alternatively, use google.cloud.bigquery.Client() directly.
        cases_df = pd.read_gbq(phenotype_query, dialect='standard', progress_bar_type=None) # None for non-interactive
        n_cases = cases_df['person_id'].nunique()
        print(f"Query complete. Found {n_cases} unique persons meeting case criteria for {phenotype_name}.\n")

        if cases_df.empty:
            print(f"WARNING: No cases found for phenotype '{phenotype_name}'. Downsampling might be affected.")
            # Return empty df with correct columns to avoid downstream errors
            return pd.DataFrame(columns=['s', 'phenotype_status'])

        cases_df = cases_df.rename(columns={'person_id': 's'})
        cases_df['phenotype_status'] = 1
        return cases_df[['s', 'phenotype_status']].astype({'s': str, 'phenotype_status': int})
    except Exception as e:
        print(f"FATAL ERROR fetching phenotype data from BigQuery: {e}")
        print("get your environment is authenticated to Google Cloud and the BigQuery API is enabled, and pandas-gbq can find a billing project.")
        sys.exit(1)


# --- VDS Preparation Functions (adapted from original script) ---
def load_full_vds(path, fs_gcs): # Pass gcsfs object
    print(f"Loading full WGS VariantDataset from: {path}")
    try:
        # Hail uses its own GCS connector, should be configured with Spark.
        # The gcsfs object is for other GCS checks, not directly for Hail VDS read.
        vds = hl.vds.read_vds(path)
        n_samples = vds.variant_data.count_cols()
        if n_samples == 0:
            print(f"FATAL ERROR: Loaded VDS from {path} contains 0 samples.")
            sys.exit(1)
        print("VDS successfully loaded.")
        print(f"Initial VDS sample count: {n_samples}\n")
        return vds
    except Exception as e:
        if "requester_pays" in str(e).lower():
            print(f"FATAL ERROR: Failed to load VDS from {path} due to requester pays issue. Check Hail/Spark GCS connector config. Error: {e}")
        else:
            print(f"FATAL ERROR: Failed to load or verify VDS from {path}: {e}")
        sys.exit(1)

def load_excluded_samples_ht(path, flagged_samples_gcs_path_default, fs_gcs): # Pass gcsfs object
    print(f"Importing flagged (related or excluded) samples from: {path}")
    if not gcs_path_exists(path): # Uses common_utils gcs_path_exists which uses the fs_gcs
        error_message = f"ERROR: Flagged samples file not found or accessible at {path}. "
        if path == flagged_samples_gcs_path_default:
            error_message += "This is the default All of Us path. The file might have moved or permissions changed."
        else:
            error_message += "This is a custom path. Please check the path and permissions."
        print(error_message)
        print("CRITICAL INFO: Proceeding without relatedness/flagged sample filtering. This may impact analysis results.")
        return None
    try:
        ht = hl.import_table(path, key='sample_id', impute=True)
        count = ht.count()
        print(f"Flagged samples loaded. Count: {count}\n")
        if count == 0:
            print("WARNING: Flagged samples table is empty.")
        return ht
    except Exception as e:
        print(f"ERROR: Failed to load flagged samples from {path}: {e}")
        print("Proceeding without relatedness filtering due to load error.")
        return None

def filter_samples_by_exclusion_list(vds, excluded_ht):
    if excluded_ht is None:
        print("Skipping relatedness/flagged sample filtering as excluded table is not available.")
        return vds
    print("Filtering flagged samples out of the VDS...")
    try:
        n_before = vds.variant_data.count_cols()
        # VDS sample key is 's' (string). Excluded HT key is 'sample_id' (imputed, often int or str).
        # Key excluded_ht by 's' as string for robust joining.
        if 's' not in excluded_ht.key:
             excluded_ht = excluded_ht.annotate(s=hl.str(excluded_ht.sample_id)).key_by('s')

        cleaned_variant_data = vds.variant_data.filter_cols(
            hl.is_missing(excluded_ht[vds.variant_data.s])
        )
        cleaned_ref_data = vds.reference_data.filter_cols(
            hl.is_missing(excluded_ht[vds.reference_data.s])
        )
        cleaned_vds = hl.vds.VariantDataset(cleaned_ref_data, cleaned_variant_data)
        n_after = cleaned_vds.variant_data.count_cols()
        print(f"Samples filtered. Count before: {n_before}, Count after: {n_after}. Removed: {n_before - n_after}\n")
        if n_after == 0:
            print("FATAL ERROR: Filtering removed all samples.")
            sys.exit(1)
        return cleaned_vds
    except Exception as e:
        print(f"ERROR: Failed during sample filtering: {e}")
        print("Exiting due to failure during sample filtering.")
        sys.exit(1)

@cache_result("wgs_ehr_samples_df")
def get_wgs_ehr_samples_df(cdr_env_var_value, prs_id=None): # prs_id not used here
    """Fetches person_ids with both WGS and EHR data."""
    get_cache_dir()
    print("Retrieving WGS+EHR samples from BigQuery...")
    if not cdr_env_var_value:
        print(f"FATAL ERROR: Workspace CDR value not set.")
        sys.exit(1)

    wgs_ehr_query = f"""
    SELECT person_id
    FROM `{cdr_env_var_value}.person`
    WHERE person_id IN (
        SELECT DISTINCT person_id
        FROM `{cdr_env_var_value}.cb_search_person`
        WHERE has_ehr_data = 1
        AND has_whole_genome_variant = 1
    )
    """
    print("Executing BQ query for WGS+EHR samples.")
    try:
        df = pd.read_gbq(wgs_ehr_query, dialect='standard', progress_bar_type=None)
        n_found = df['person_id'].nunique()
        print(f"WGS+EHR query completed. Found {n_found} unique persons.\n")
        if df.empty:
            print("FATAL ERROR: No samples found with both WGS and EHR data from BigQuery.")
            sys.exit(1)
        return df
    except Exception as e:
        print(f"FATAL ERROR: Failed to query BigQuery for WGS+EHR samples: {e}")
        sys.exit(1)

def save_wgs_ehr_ids_to_gcs(df, gcs_path, fs_gcs): # Pass gcsfs object
    if df is None or df.empty:
        print("FATAL ERROR: Cannot save WGS+EHR IDs as DataFrame is None or empty.")
        sys.exit(1)
    print(f"Storing WGS+EHR sample IDs to GCS: {gcs_path}")
    try:
        with fs_gcs.open(gcs_path, 'w') as f:
            df.to_csv(f, index=False)
        print("WGS+EHR sample IDs saved.\n")
    except Exception as e:
        print(f"FATAL ERROR: Failed to save WGS+EHR IDs to {gcs_path}: {e}")
        sys.exit(1)

def filter_vds_to_wgs_ehr_list(vds, wgs_ehr_ids_gcs_path, fs_gcs): # Pass gcsfs object
    if vds is None:
        print("FATAL ERROR: VDS is None, cannot filter to WGS/EHR samples.")
        sys.exit(1)
    if wgs_ehr_ids_gcs_path is None:
        print("FATAL ERROR: WGS/EHR IDs GCS path is None, cannot filter VDS.")
        sys.exit(1)

    print(f"Importing WGS+EHR sample list from {wgs_ehr_ids_gcs_path} for VDS filtering...")
    try:
        if not gcs_path_exists(wgs_ehr_ids_gcs_path): # Uses common_utils
            print(f"FATAL ERROR: WGS EHR IDs file disappeared before import: {wgs_ehr_ids_gcs_path}")
            sys.exit(1)

        wgs_ehr_ht = hl.import_table(wgs_ehr_ids_gcs_path, delimiter=',', key='person_id', types={'person_id': hl.tstr})
        if wgs_ehr_ht.count() == 0:
            print("FATAL ERROR: Imported WGS+EHR HailTable is empty. Cannot filter VDS.")
            sys.exit(1)

        vds_variant_data = vds.variant_data
        vds_reference_data = vds.reference_data
        if vds_variant_data.s.dtype != hl.tstr:
            print("Warning: VDS sample key 's' in variant_data is not string. Attempting cast.")
            vds_variant_data = vds_variant_data.key_cols_by(s=hl.str(vds_variant_data.s))
        if vds_reference_data.s.dtype != hl.tstr:
            print("Warning: VDS sample key 's' in reference_data is not string. Attempting cast.")
            vds_reference_data = vds_reference_data.key_cols_by(s=hl.str(vds_reference_data.s))

        vds = hl.vds.VariantDataset(vds_reference_data, vds_variant_data) # Recreate VDS if keys were cast

        wgs_ehr_ht_keyed = wgs_ehr_ht.key_by(s=wgs_ehr_ht.person_id) # Key by 's'

        print("Filtering VDS variant_data to WGS+EHR list...")
        filtered_variant_data = vds.variant_data.semi_join_cols(wgs_ehr_ht_keyed)
        print("Filtering VDS reference_data to WGS+EHR list...")
        filtered_reference_data = vds.reference_data.semi_join_cols(wgs_ehr_ht_keyed)
        subset_vds = hl.vds.VariantDataset(filtered_reference_data, filtered_variant_data)

        n_after_filter = subset_vds.variant_data.count_cols()
        print(f"VDS filtered to WGS+EHR samples. Final count: {n_after_filter}\n")
        if n_after_filter == 0:
            print("FATAL ERROR: Filtering VDS to WGS+EHR samples resulted in 0 samples remaining.")
            sys.exit(1)
        return subset_vds
    except Exception as e:
        print(f"FATAL ERROR: Failed filtering VDS to WGS+EHR samples: {e}")
        sys.exit(1)


def main():
    args = parse_args()
    get_cache_dir() # For pandas caching if used by helper functions

    # Initialize GCS FileSystem
    fs = get_gcs_fs() # fs object from common_utils

    # Initialize Hail
    init_hail(args.gcs_hail_temp_dir, args.run_timestamp)

    base_cohort_vds = None

    if hail_path_exists(args.base_cohort_vds_path_out): # Uses common_utils
        print(f"[CHECKPOINT] Found potential Base Cohort VDS directory at: {args.base_cohort_vds_path_out}")
        try:
            print("Attempting to read VDS checkpoint...")
            base_cohort_vds = hl.vds.read_vds(args.base_cohort_vds_path_out)
            print("Performing sanity check on loaded VDS...")
            n_samples = base_cohort_vds.variant_data.count_cols()
            n_variants = base_cohort_vds.variant_data.count_rows()
            print(f"Sanity check: Found {n_samples} samples and {n_variants} variants.")
            if n_samples > 0 and n_variants > 0:
                print("[CHECKPOINT HIT] Successfully loaded and verified VDS checkpoint.\n")
            else:
                print("[CHECKPOINT CORRUPTED] Loaded VDS has 0 samples or 0 variants. Assuming invalid.")
                base_cohort_vds = None
                delete_gcs_path(args.base_cohort_vds_path_out, recursive=True) # Uses common_utils
        except Exception as e:
            error_message = str(e)
            print(f"[CHECKPOINT CORRUPTED] Failed to read VDS from {args.base_cohort_vds_path_out}. Error: {error_message}")
            if "corrupted" in error_message.lower() or "metadata.json.gz is missing" in error_message or "Not a VDS" in error_message:
                print("Error suggests checkpoint corruption. Deleting and regenerating...")
                base_cohort_vds = None
                delete_gcs_path(args.base_cohort_vds_path_out, recursive=True) # Uses common_utils
            else:
                print("Error was not identified as typical corruption. Exiting to investigate.")
                sys.exit(1)

    if base_cohort_vds is None:
        print(f"[CHECKPOINT MISS or CORRUPTED] Base Cohort VDS needs to be generated at {args.base_cohort_vds_path_out}.")
        print("--- Starting Common Data Preparation ---")

        full_vds = load_full_vds(args.wgs_vds_path, fs)
        excluded_ht = load_excluded_samples_ht(args.flagged_samples_gcs_path, args.flagged_samples_gcs_path, fs) # Pass default for comparison
        cleaned_vds = filter_samples_by_exclusion_list(full_vds, excluded_ht)
        del full_vds, excluded_ht

        all_wgs_ehr_individuals_df = get_wgs_ehr_samples_df(args.workspace_cdr)
        all_wgs_ehr_individuals_df = all_wgs_ehr_individuals_df.rename(columns={'person_id': 's'})
        all_wgs_ehr_individuals_df['s'] = all_wgs_ehr_individuals_df['s'].astype(str)

        final_ids_for_vds_df = None
        if args.enable_downsampling_for_vds:
            print(f"--- Downsampling WGS+EHR cohort to ~{args.n_cases_downsample} cases and ~{args.n_controls_downsample} controls for VDS generation ---")
            phenotype_df_for_sampling = get_phenotype_data_for_vds_prep(
                args.target_phenotype_name, args.phenotype_concept_ids, args.workspace_cdr
            ) # This is cached
            phenotype_df_for_sampling['s'] = phenotype_df_for_sampling['s'].astype(str)

            merged_cohort_phenotype_df = pd.merge(
                all_wgs_ehr_individuals_df[['s']],
                phenotype_df_for_sampling[['s', 'phenotype_status']],
                on='s',
                how='left'
            )
            merged_cohort_phenotype_df['phenotype_status'] = merged_cohort_phenotype_df['phenotype_status'].fillna(0).astype(int)
            print(f"Total WGS+EHR individuals considered for downsampling: {len(merged_cohort_phenotype_df)}")
            print(f"Phenotype distribution for downsampling:\n{merged_cohort_phenotype_df['phenotype_status'].value_counts(dropna=False)}")

            cases_df = merged_cohort_phenotype_df[merged_cohort_phenotype_df['phenotype_status'] == 1]
            controls_df = merged_cohort_phenotype_df[merged_cohort_phenotype_df['phenotype_status'] == 0]
            print(f"Available cases for downsampling: {len(cases_df)}")
            print(f"Available controls for downsampling: {len(controls_df)}")

            num_cases_to_sample = min(args.n_cases_downsample, len(cases_df))
            sampled_cases_df = cases_df.sample(n=num_cases_to_sample, random_state=args.downsampling_random_state, replace=False)
            print(f"Sampled {len(sampled_cases_df)} cases.")

            num_controls_to_sample = min(args.n_controls_downsample, len(controls_df))
            sampled_controls_df = controls_df.sample(n=num_controls_to_sample, random_state=args.downsampling_random_state, replace=False)
            print(f"Sampled {len(sampled_controls_df)} controls.")

            final_ids_for_vds_df = pd.concat([sampled_cases_df, sampled_controls_df], ignore_index=True)[['s']]
            final_ids_for_vds_df = final_ids_for_vds_df.drop_duplicates(subset=['s'])
            print(f"Total individuals in downsampled cohort for VDS generation: {len(final_ids_for_vds_df)}")
            if final_ids_for_vds_df.empty:
                print("FATAL ERROR: Downsampling resulted in an empty cohort.")
                sys.exit(1)
        else:
            print("--- Using full WGS+EHR cohort for VDS generation (downsampling disabled) ---")
            final_ids_for_vds_df = all_wgs_ehr_individuals_df[['s']]

        df_to_save_for_vds_filter = final_ids_for_vds_df.rename(columns={'s': 'person_id'})
        save_wgs_ehr_ids_to_gcs(df_to_save_for_vds_filter, args.wgs_ehr_ids_gcs_path_out, fs) # Uses common_utils
        del all_wgs_ehr_individuals_df, final_ids_for_vds_df, df_to_save_for_vds_filter
        if 'phenotype_df_for_sampling' in locals(): del phenotype_df_for_sampling
        if 'cases_df' in locals(): del cases_df
        if 'controls_df' in locals(): del controls_df
        if 'sampled_cases_df' in locals(): del sampled_cases_df
        if 'sampled_controls_df' in locals(): del sampled_controls_df
        if 'merged_cohort_phenotype_df' in locals(): del merged_cohort_phenotype_df

        base_cohort_vds = filter_vds_to_wgs_ehr_list(cleaned_vds, args.wgs_ehr_ids_gcs_path_out, fs) # Uses common_utils
        del cleaned_vds

        print(f"Writing filtered Base Cohort VDS checkpoint to: {args.base_cohort_vds_path_out}")
        try:
            base_cohort_vds.write(args.base_cohort_vds_path_out, overwrite=True)
            print("Base Cohort VDS checkpoint successfully written.")
            print("Verifying written checkpoint...")
            vds_check = hl.vds.read_vds(args.base_cohort_vds_path_out)
            vds_check_sample_count = vds_check.variant_data.count_cols() # Perform a light verification
            print(f"Checkpoint verified successfully. Sample count: {vds_check_sample_count}")
        except Exception as e:
            print(f"ERROR: Failed to write or verify Base Cohort VDS checkpoint: {e}")
            print("FATAL ERROR: Checkpoint write/verification failed.")
            delete_gcs_path(args.base_cohort_vds_path_out, recursive=True) # Uses common_utils
            sys.exit(1)
        print("--- Common Data Preparation Finished ---")

    if base_cohort_vds is None:
        print(f"FATAL ERROR: base_cohort_vds is None after attempting load/regeneration from {args.base_cohort_vds_path_out}. Cannot proceed.")
        sys.exit(1)
    else:
        final_sample_count = base_cohort_vds.variant_data.count_cols()
        if final_sample_count == 0:
            print(f"FATAL ERROR: base_cohort_vds at {args.base_cohort_vds_path_out} has 0 samples.")
            sys.exit(1)
        print(f"Base Cohort VDS at {args.base_cohort_vds_path_out} is ready with {final_sample_count} samples.\n")

    # The script's job is to produce args.base_cohort_vds_path_out and args.wgs_ehr_ids_gcs_path_out
    # Nextflow will know these paths.
    print(f"Script completed. Base VDS written to: {args.base_cohort_vds_path_out}")
    print(f"Script completed. WGS+EHR IDs written to: {args.wgs_ehr_ids_gcs_path_out}")

if __name__ == "__main__":
    main()
