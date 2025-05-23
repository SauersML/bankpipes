import argparse
import os
import sys
import pandas as pd
from utils import cache_result, make_cache_dir, gcs_path_exists, get_gcs_fs

@cache_result("phenotype_data_main") # Cache name prefix
def get_phenotype_data(phenotype_name, concept_ids_str, cdr_env_var_value, prs_id=None): # prs_id not used here, for cache_result consistency
    """Fetches case status based on OMOP concept IDs."""
    get_cache_dir() # Changed from make_cache_dir
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
        cases_df = pd.read_gbq(phenotype_query, dialect='standard', progress_bar_type=None)
        n_cases = cases_df['person_id'].nunique()
        print(f"Query complete. Found {n_cases} unique persons meeting case criteria for {phenotype_name}.\n")

        if cases_df.empty:
            print(f"WARNING: No cases found for phenotype '{phenotype_name}'. Analysis might be limited.")
            return pd.DataFrame(columns=['s', 'phenotype_status']) # Return empty with schema

        cases_df = cases_df.rename(columns={'person_id': 's'})
        cases_df['phenotype_status'] = 1
        return cases_df[['s', 'phenotype_status']].astype({'s': str, 'phenotype_status': int})
    except Exception as e:
        print(f"FATAL ERROR fetching phenotype data from BigQuery: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Fetch phenotype data from BigQuery.")
    parser.add_argument("--phenotype_name", required=True, help="Name of the target phenotype.")
    parser.add_argument("--phenotype_concept_ids", required=True, help="Comma-separated OMOP concept IDs.")
    parser.add_argument("--workspace_cdr", required=True, help="Workspace CDR string.")
    parser.add_argument("--output_phenotype_csv_gcs_path", required=True, help="GCS path to save the phenotype CSV.")
    parser.add_argument("--google_billing_project", required=True, help="Google Cloud Project ID for billing and GCS access.")

    args = parser.parse_args()
    get_cache_dir() # Ensure local cache directory is available for the @cache_result decorator; changed from make_cache_dir

    # Initialize GCS filesystem, explicitly providing the billing project
    fs = get_gcs_fs(project_id_for_billing=args.google_billing_project)

    phenotype_cases_df = get_phenotype_data(
        phenotype_name=args.phenotype_name,
        concept_ids_str=args.phenotype_concept_ids,
        cdr_env_var_value=args.workspace_cdr
    )

    if phenotype_cases_df is None: # Should not happen due to sys.exit in function
        print(f"FATAL: Phenotype data for {args.phenotype_name} could not be fetched.")
        sys.exit(1)

    print(f"--- Phenotype Definition for {args.phenotype_name} Finished ---")
    print(f"Number of cases found: {len(phenotype_cases_df[phenotype_cases_df['phenotype_status']==1])}")

    # Save the fetched phenotype data to the specified GCS path
    print(f"Saving phenotype data to GCS: {args.output_phenotype_csv_gcs_path}")
    try:
        # make parent GCS directory exist
        output_dir_gcs = os.path.dirname(args.output_phenotype_csv_gcs_path)
        if not gcs_path_exists(output_dir_gcs): # common_utils
            print(f"Creating GCS directory: {output_dir_gcs}")
            fs.mkdirs(output_dir_gcs, exist_ok=True)

        with fs.open(args.output_phenotype_csv_gcs_path, 'w') as f:
            phenotype_cases_df.to_csv(f, index=False, na_rep='NA')
        print(f"Phenotype data successfully saved to {args.output_phenotype_csv_gcs_path}")
    except Exception as e:
        print(f"FATAL ERROR: Failed to save phenotype data to {args.output_phenotype_csv_gcs_path}: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
