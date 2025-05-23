import os
from PheTK import PheWAS

# --- Configuration ---

# 1. Input File Paths
# This CSV file should contain:
#   - 'person_id'
#   - A column for independent variable (e.g., SV genotype).
#   - Columns for covariates (e.g., PCs, age).
#   - A column for sex at birth.
COHORT_DATA_CSV = "cohort_with_sv_and_pcs.csv"

# This CSV file should contain pre-calculated phecode counts for cohort:
#   - 'person_id'
#   - 'phecode' (the phecode string)
#   - 'count' (number of occurrences of the phecode for the person_id)
PHECODE_COUNTS_CSV = "precomputed_phecode_counts.csv"

# 2. Output File Name
PHEWAS_RESULTS_OUTPUT_CSV = "custom_phewas_run_results.csv"

# 3. PheWAS Core Parameters
PHECODE_VERSION_TO_USE = "X"  # Common options: "X" or "1.2"

# Column name in COHORT_DATA_CSV that holds main variable of interest (e.g., SV genotype)
# This variable should be numeric (e.g., 0 for control, 1 for case/has_SV).
INDEPENDENT_VARIABLE_COLUMN_NAME = "SV_STATUS"

# Column name in COHORT_DATA_CSV that holds sex at birth information
SEX_AT_BIRTH_COLUMN_NAME = "SEX_CODE" # Example: "SEX_CODE"
# Specify how sex is coded in the SEX_AT_BIRTH_COLUMN_NAME column
# True if males are coded as 1 and females as 0.
# False if males are coded as 0 and females as 1.
MALE_CODED_AS_ONE = True

# List of column names from COHORT_DATA_CSV to be used as covariates.
# These should include custom PCs and any other desired covariates like age.
# All covariate columns must contain numeric data.
COVARIATE_COLUMN_NAMES = [
    "PC1", "PC2", "PC3", "PC4", "PC5",
    "PC6", "PC7", "PC8", "PC9", "PC10",
    "AGE_AT_LAST_EVENT"
]

# 4. PheWAS Thresholds and Options
MINIMUM_CASES_PER_PHECODE = 20     # Minimum number of cases for a phecode to be tested.
MINIMUM_PHECODE_COUNT_FOR_CASE = 2 # Minimum occurrences of a phecode for an individual to be defined as a case.

# For Phecode v1.2, True to use standard phecode exclusion criteria. Ignored for Phecode v"X".
USE_PHECODE_EXCLUSIONS = False

# 5. Computational Settings
# Number of threads for parallel processing.
# Uses 2/3 of available CPUs by default in PheTK, adjust if needed.
try:
    NUM_THREADS_TO_USE = max(1, round(os.cpu_count() * 2 / 3) -1 ) # Default-ish from PheTK
    if NUM_THREADS_TO_USE == 0: NUM_THREADS_TO_USE = 1
except TypeError: # os.cpu_count() can return None
    NUM_THREADS_TO_USE = 4 # Fallback to a sensible number
    print(f"Warning: Could not determine CPU count, defaulting to {NUM_THREADS_TO_USE} threads.")


# --- Main PheWAS Execution ---

print(f"Starting PheWAS analysis with PheTK.")
print(f"Cohort data: {COHORT_DATA_CSV}")
print(f"Phecode counts: {PHECODE_COUNTS_CSV}")
print(f"Independent variable: {INDEPENDENT_VARIABLE_COLUMN_NAME}")
print(f"Covariates: {', '.join(COVARIATE_COLUMN_NAMES)}")

try:
    # Instantiate the PheWAS class
    phewas_runner = PheWAS(
        phecode_version=PHECODE_VERSION_TO_USE,
        phecode_count_csv_path=PHECODE_COUNTS_CSV,
        cohort_csv_path=COHORT_DATA_CSV,
        sex_at_birth_col=SEX_AT_BIRTH_COLUMN_NAME,
        covariate_cols=COVARIATE_COLUMN_NAMES,
        independent_variable_of_interest=INDEPENDENT_VARIABLE_COLUMN_NAME,
        male_as_one=MALE_CODED_AS_ONE,
        icd_version="US", # Primarily for loading phecode descriptions/metadata.
                          # Less critical if phecode_counts are already perfect.
        min_cases=MINIMUM_CASES_PER_PHECODE,
        min_phecode_count=MINIMUM_PHECODE_COUNT_FOR_CASE,
        use_exclusion=USE_PHECODE_EXCLUSIONS if PHECODE_VERSION_TO_USE == "1.2" else False,
        output_file_name=PHEWAS_RESULTS_OUTPUT_CSV,
        verbose=True # Provides more detailed console output during the run.
    )

    # Execute the PheWAS analysis
    phewas_runner.run(
        parallelization="multithreading", # This is currently the main supported method in PheTK source
        n_threads=NUM_THREADS_TO_USE
    )

    print(f"\nPheWAS analysis successfully completed.")
    print(f"Results saved to: {phewas_runner.output_file_name}") # PheTK might modify the output_file_name (e.g. add timestamp)
    print(f"Number of phecodes considered from input: {len(phewas_runner.phecode_list)}")
    print(f"Number of phecodes tested: {phewas_runner.tested_count}")
    if phewas_runner.tested_count > 0:
        print(f"Suggested Bonferroni correction threshold (-log10 scale): {phewas_runner.bonferroni:.4f}")
        print(f"Number of phecodes exceeding Bonferroni threshold: {phewas_runner.above_bonferroni_count}")
    else:
        print("No phecodes were tested. Please check input data and parameters (e.g., min_cases).")

except FileNotFoundError as e:
    print(f"\nERROR: An essential input file was not found.")
    print(f"Details: {e}")
    print(f"Please ensure '{COHORT_DATA_CSV}' and '{PHECODE_COUNTS_CSV}' exist at the specified paths.")
except KeyError as e:
    print(f"\nERROR: A specified column name was not found in CSV files.")
    print(f"Column details: {e}")
    print(f"Please check that '{INDEPENDENT_VARIABLE_COLUMN_NAME}', '{SEX_AT_BIRTH_COLUMN_NAME}', and all listed covariate names in '{COVARIATE_COLUMN_NAMES}' exist in '{COHORT_DATA_CSV}'.")
except Exception as e:
    print(f"\nAn unexpected error occurred during the PheWAS analysis:")
    print(f"Details: {e}")
    print("Please review input data formats and parameters.")

print("\nScript finished.")
