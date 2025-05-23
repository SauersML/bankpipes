"""
process_prs_model.py

Nextflow task script to:
  1. Download and prepare PRS weight table (checking for existing intermediates)
  2. Extract genomic intervals (checking for existing intermediates)
  3. Load base VDS and filter by those intervals
  4. Load prepared weights and annotate VDS
  5. Compute per-sample PRS scores using correct dosage calculation
  6. Save final Hail Table & CSV to GCS

Uses helper functions from utils.py.
"""

import argparse
import os
import sys
import pandas as pd
import hail as hl
from utils import (
    init_hail,
    get_gcs_fs,
    gcs_path_exists,
    delete_gcs_path,
    download_file
)

# --- Core PRS Calculation Logic ---

def prepare_weight_table(harmonized_path: str, prs_id: str) -> pd.DataFrame | None:
    """
    Reads, cleans (ensuring integer positions), and prepares the PRS weight table DataFrame.
    Returns the prepared DataFrame or None on failure.
    """
    if not os.path.exists(harmonized_path):
        print(f"[{prs_id}] ERROR: Missing harmonized file for preparation at {harmonized_path}")
        return None

    print(f"[{prs_id}] Reading harmonized PRS weight file: {harmonized_path}")
    try:
        score_df = pd.read_csv(harmonized_path, sep='\t', comment='#', low_memory=False)
        print(f"[{prs_id}] PRS weight data loaded. Initial variant count: {len(score_df)}")

        # Define comprehensive column mapping
        col_map = {
            'hm_chr': 'chr', 'chr': 'chr', 'chromosome': 'chr',
            'hm_pos': 'bp', 'pos': 'bp', 'base_pair_location': 'bp', 'position': 'bp', 'BP': 'bp',
            'effect_allele': 'effect_allele', 'allele': 'effect_allele', 'A1': 'effect_allele', 'hm_effect_allele': 'effect_allele',
            'other_allele': 'other_allele', 'A2': 'other_allele', 'hm_other_allele': 'other_allele', 'reference_allele': 'other_allele',
            'effect_weight': 'weight', 'beta': 'weight', 'OR': 'weight', 'effect': 'weight', 'hm_beta': 'weight'
        }
        rename_dict = {}
        found_cols = []
        needed_cols = ['chr', 'bp', 'effect_allele', 'weight'] # Essential columns
        # Optional: add 'other_allele' if needed for validation later, but not core calculation
        optional_cols = ['other_allele']
        all_needed_cols = needed_cols + optional_cols

        for potential_name, target_name in col_map.items():
             if potential_name in score_df.columns and target_name in all_needed_cols:
                 if target_name not in found_cols:
                     rename_dict[potential_name] = target_name
                     found_cols.append(target_name)

        print(f"[{prs_id}] Renaming columns: {rename_dict}")
        score_df = score_df.rename(columns=rename_dict)

        missing_essential = [col for col in needed_cols if col not in score_df.columns]
        if missing_essential:
            print(f"[{prs_id}] ERROR: Missing essential columns after renaming: {missing_essential}. Available columns: {score_df.columns.tolist()}")
            return None

        cols_to_keep = [col for col in all_needed_cols if col in score_df.columns]
        print(f"[{prs_id}] Selecting and cleaning required columns: {cols_to_keep}")
        score_df = score_df[cols_to_keep].copy() # Use copy to avoid SettingWithCopyWarning

        # string types for alleles and chr before cleaning numerics
        score_df['chr'] = score_df['chr'].astype(str)
        score_df['effect_allele'] = score_df['effect_allele'].astype(str)
        if 'other_allele' in score_df.columns:
            score_df['other_allele'] = score_df['other_allele'].astype(str)

        # Robust numeric conversion and cleaning
        score_df['bp'] = pd.to_numeric(score_df['bp'], errors='coerce')
        score_df['weight'] = pd.to_numeric(score_df['weight'], errors='coerce')

        # Drop rows where essential columns are now NaN or invalid strings
        initial_len = len(score_df)
        essential_for_dropna = ['bp', 'weight', 'effect_allele', 'chr']
        score_df = score_df.dropna(subset=essential_for_dropna)
        # Additional check for empty strings which might pass dropna
        score_df = score_df[score_df['chr'].str.strip() != '']
        score_df = score_df[score_df['effect_allele'].str.strip() != '']
        # Optional: Basic allele check (can be expanded)
        # score_df = score_df[score_df['effect_allele'].str.match('^[ACGT]+$', case=False, na=False)]

        dropped_count = initial_len - len(score_df)
        if dropped_count > 0:
            print(f"[{prs_id}] Warning: Dropped {dropped_count} rows during cleaning (missing/invalid pos, weight, allele, chr).")

        if score_df.empty:
            print(f"[{prs_id}] ERROR: DataFrame became empty after cleaning. Check input file format/content.")
            return None

        # Convert bp to integer *after* cleaning
        try:
            score_df['bp'] = score_df['bp'].astype(int)
        except ValueError as e:
             print(f"[{prs_id}] ERROR: Failed to convert 'bp' column to integer after cleaning. This should not happen. Error: {e}")
             return None

        # Add Hail-compatible columns
        print(f"[{prs_id}] Adding Hail-compatible columns: contig, position")
        # chr prefixes are handled
        score_df['contig'] = score_df['chr'].apply(lambda c: c if c.startswith('chr') else f'chr{c}')
        # only valid contigs remain (e.g., chr1-22, chrX, chrY, chrM) - basic check example
        valid_contigs = {f'chr{i}' for i in range(1, 23)} | {'chrX', 'chrY', 'chrM'}
        original_count = len(score_df)
        score_df = score_df[score_df['contig'].isin(valid_contigs)]
        if len(score_df) < original_count:
            print(f"[{prs_id}] WARNING: Dropped {original_count - len(score_df)} rows with non-standard contigs.")

        score_df['position'] = score_df['bp'] # Now guaranteed integer

        final_count = len(score_df)
        print(f"[{prs_id}] PRS weight table prepared. Final variant count: {final_count}.\n")
        if final_count == 0:
            print(f"[{prs_id}] ERROR: Prepared weight table has 0 variants after processing.")
            return None

        # Return only necessary columns for Hail import + interval extraction
        cols_to_return = ['contig', 'position', 'effect_allele', 'weight', 'bp']
        if 'other_allele' in score_df.columns:
            cols_to_return.append('other_allele')
        return score_df[cols_to_return]

    except Exception as e:
        print(f"[{prs_id}] ERROR preparing weight table: {e}")
        if 'score_df' in locals() and isinstance(score_df, pd.DataFrame):
             print(f"Columns at error: {score_df.columns.tolist()}")
        return None


def save_to_gcs(df: pd.DataFrame, prs_id: str, local_name: str, gcs_path: str, fs) -> bool:
    """Save a DataFrame locally then upload to GCS."""
    # Use a temporary local file name, potentially adding PID for parallelism safety if needed, though Nextflow tasks are isolated.
    tmp_local_path = f"{prs_id}_{local_name}_temp.csv"
    try:
        df.to_csv(tmp_local_path, index=False, na_rep='NA') # Explicit NA representation
        print(f"[{prs_id}] Saved {local_name} locally to {tmp_local_path}, uploading to {gcs_path}")

        parent_gcs_dir = os.path.dirname(gcs_path)
        if not gcs_path_exists(parent_gcs_dir): # Use helper from utils
            print(f"[{prs_id}] Creating GCS directory: {parent_gcs_dir}")
            fs.mkdirs(parent_gcs_dir, exist_ok=True)

        fs.put(tmp_local_path, gcs_path)
        print(f"[{prs_id}] Successfully uploaded {tmp_local_path} to {gcs_path}")
        os.remove(tmp_local_path) # Clean up local file after successful upload
        return True

    except Exception as e:
        print(f"[{prs_id}] ERROR saving or uploading {local_name} to {gcs_path}: {e}")
        # Clean up local file if it exists
        if os.path.exists(tmp_local_path):
            try: os.remove(tmp_local_path)
            except OSError: pass
        # Attempt to delete potentially incomplete GCS file
        delete_gcs_path(gcs_path, recursive=False) # Use helper from utils
        return False


def extract_intervals(df: pd.DataFrame, prs_id: str, gcs_path: str, fs) -> bool:
    """Extract contig, position, end TSV for Hail interval filtering and upload to GCS."""
    if df is None or df.empty:
        print(f"[{prs_id}] ERROR: Cannot extract intervals from empty or None DataFrame.")
        return False
    if not all(col in df.columns for col in ['contig', 'position']):
         print(f"[{prs_id}] ERROR: DataFrame missing 'contig' or 'position' for interval extraction.")
         return False

    print(f"[{prs_id}] Extracting variant intervals...")
    try:
        intervals_df = df[['contig', 'position']].copy()
        intervals_df['end'] = intervals_df['position'] # Create end column, inclusive point interval
        intervals_df = intervals_df[['contig', 'position', 'end']].drop_duplicates()

        if intervals_df.empty:
            print(f"[{prs_id}] ERROR: No unique intervals generated.")
            return False

        print(f"[{prs_id}] Generated {len(intervals_df)} unique intervals.")
        # Save locally temporarily
        tmp_local_tsv = f"{prs_id}_intervals_temp.tsv"
        intervals_df.to_csv(tmp_local_tsv, sep='\t', index=False, header=False)

        print(f"[{prs_id}] Saved intervals locally to {tmp_local_tsv}, uploading to {gcs_path}")
        parent_gcs_dir = os.path.dirname(gcs_path)
        if not gcs_path_exists(parent_gcs_dir): # Use helper from utils
            print(f"[{prs_id}] Creating GCS directory: {parent_gcs_dir}")
            fs.mkdirs(parent_gcs_dir, exist_ok=True)

        fs.put(tmp_local_tsv, gcs_path)
        print(f"[{prs_id}] Successfully uploaded intervals to {gcs_path}")
        os.remove(tmp_local_tsv) # Clean up local file
        return True

    except Exception as e:
        print(f"[{prs_id}] ERROR extracting or uploading intervals to {gcs_path}: {e}")
        # Clean up local file if it exists
        if 'tmp_local_tsv' in locals() and os.path.exists(tmp_local_tsv):
             try: os.remove(tmp_local_tsv)
             except OSError: pass
        # Attempt to delete potentially incomplete GCS file
        delete_gcs_path(gcs_path, recursive=False) # Use helper from utils
        return False


def import_weights_as_ht(gcs_csv_path: str, prs_id: str) -> hl.Table | None:
    """Import the prepared weight CSV from GCS into a Hail Table keyed by locus."""
    if not gcs_path_exists(gcs_csv_path): # Use helper from utils
        print(f"[{prs_id}] ERROR: Weight CSV file not found on GCS at {gcs_csv_path}")
        return None

    print(f"[{prs_id}] Importing PRS weights from GCS for Hail annotation: {gcs_csv_path}")
    try:
        # Import, specifying types crucial for locus creation and calculation
        prs_ht = hl.import_table(
            gcs_csv_path,
            delimiter=',',
            types={'contig': hl.tstr, 'position': hl.tint32, 'weight': hl.tfloat64, 'bp': hl.tint32},
            missing='NA',
            comment='#', # Should not be needed if prepared correctly, but safe.
            impute=True # Impute other columns like effect_allele
        )

        # --- Verification after import ---
        if prs_ht.count() == 0:
             print(f"[{prs_id}] ERROR: Imported PRS HailTable from {gcs_csv_path} has 0 rows.")
             return None
        required_fields = {'contig': hl.tstr, 'position': hl.tint32, 'weight': hl.tfloat64, 'effect_allele': hl.tstr}
        missing_or_wrong_type = []
        for field, expected_type in required_fields.items():
             if field not in prs_ht.row:
                  missing_or_wrong_type.append(f"'{field}' missing")
             elif prs_ht[field].dtype != expected_type:
                  # Allow float64->float32 implicitly if needed, but warn? Strict check here.
                  missing_or_wrong_type.append(f"'{field}' wrong type ({prs_ht[field].dtype}, expected {expected_type})")
        if missing_or_wrong_type:
             print(f"[{prs_id}] ERROR: Imported PRS HT has missing or incorrect types: {', '.join(missing_or_wrong_type)}. Schema: {prs_ht.row.dtype}")
             prs_ht.show(5)
             return None
        # --- End Verification ---

        # Create locus and key by it
        prs_ht = prs_ht.filter(hl.is_defined(prs_ht.contig) & hl.is_defined(prs_ht.position))
        prs_ht = prs_ht.annotate(locus=hl.locus(prs_ht.contig, prs_ht.position, reference_genome='GRCh38'))
        prs_ht = prs_ht.filter(hl.is_defined(prs_ht.locus)) # Filter out rows where locus creation failed

        n_rows_after_locus = prs_ht.count()
        if n_rows_after_locus == 0:
            print(f"[{prs_id}] ERROR: No rows remaining after locus creation/filtering.")
            return None

        prs_ht = prs_ht.key_by('locus')
        # Keep only necessary fields for annotation
        fields_to_keep = ['weight', 'effect_allele']
        if 'other_allele' in prs_ht.row: # Keep other allele if present for potential later checks
             fields_to_keep.append('other_allele')
        prs_ht = prs_ht.select(*fields_to_keep)

        print(f"[{prs_id}] PRS HT ready for annotation ({n_rows_after_locus} variants). Schema: {prs_ht.row.dtype}\n")
        return prs_ht

    except Exception as e:
        print(f"[{prs_id}] ERROR loading prepared PRS table '{gcs_csv_path}' into Hail: {e}")
        return None


def filter_vds_by_intervals(vds: hl.vds.VariantDataset, gcs_intervals_path: str, prs_id: str) -> hl.vds.VariantDataset | None:
    """Filters VDS to loci within the specified intervals TSV file on GCS."""
    if vds is None:
        print(f"[{prs_id}] ERROR: Input VDS is None, cannot filter.")
        return None
    if not gcs_path_exists(gcs_intervals_path): # Use helper from utils
        print(f"[{prs_id}] ERROR: Interval file missing at {gcs_intervals_path}, cannot filter VDS.")
        return None

    print(f"[{prs_id}] Filtering VDS by PRS intervals from: {gcs_intervals_path}")
    try:
        print(f"[{prs_id}] Importing locus intervals...")
        locus_intervals = hl.import_locus_intervals(
            gcs_intervals_path,
            reference_genome='GRCh38',
            skip_invalid_intervals=True # Important for robustness
        )
        interval_count = locus_intervals.count()
        if interval_count == 0:
            print(f"[{prs_id}] WARNING: No valid intervals loaded from {gcs_intervals_path}. Filtering will likely remove all variants.")
            # Proceeding will result in empty MT, which is handled later, but good to warn.
        else:
            print(f"[{prs_id}] Loaded {interval_count} intervals for filtering.")

        print(f"[{prs_id}] Applying interval filter to VDS variant data...")
        # Filtering VDS directly is more efficient than converting to MT first if only using variant_data
        filtered_variant_data = hl.filter_intervals(vds.variant_data, locus_intervals)
        filtered_ref_data = hl.filter_intervals(vds.reference_data, locus_intervals) # Filter ref data too
        vds_filtered = hl.vds.VariantDataset(filtered_ref_data, filtered_variant_data)

        n_variants_after = vds_filtered.variant_data.count_rows()
        print(f"[{prs_id}] VDS filtered to PRS intervals. Variants remaining in variant_data: {n_variants_after}\n")
        if n_variants_after == 0:
            print(f"[{prs_id}] WARNING: No variants remained in variant_data after interval filtering.")
            # This is a valid outcome if there's no overlap. Score computation will handle this.

        return vds_filtered
    except Exception as e:
        print(f"[{prs_id}] ERROR filtering VDS by intervals using {gcs_intervals_path}: {e}")
        return None

# --- CORRECT DOSAGE CALCULATION FUNCTION (from original notebook) ---
def calculate_effect_allele_dosage(mt_row):
    """
    Calculates effect allele dosage from Hail MatrixTable row.
    Assumes mt_row contains .alleles, .GT, and .prs.effect_allele fields.
    """
    # required fields are present (add checks if needed upstream)
    eff_allele = hl.str(mt_row.prs.effect_allele) # Get effect allele from annotated prs info
    ref_allele = mt_row.alleles[0]

    # Find the index of the effect allele if it's an alternate allele
    alt_indices = hl.range(hl.len(mt_row.alleles) - 1).map(
        lambda i: hl.struct(idx=i + 1, allele=mt_row.alleles[i+1])
    )
    effect_allele_alt_index_struct = hl.find(lambda x: x.allele == eff_allele, alt_indices)
    # Extract the index itself, will be missing if not found
    effect_allele_alt_index = hl.if_else(hl.is_defined(effect_allele_alt_index_struct),
                                          effect_allele_alt_index_struct.idx,
                                          hl.missing(hl.tint32))

    is_effect_ref = (ref_allele == eff_allele)
    is_effect_alt = hl.is_defined(effect_allele_alt_index)

    # Calculate dosage based on genotype and effect allele identity
    dosage = hl.cond(
        hl.is_missing(mt_row.GT),
        hl.missing(hl.tint32), # Return missing dosage if GT is missing
        hl.cond(
            mt_row.GT.is_hom_ref(),
            hl.if_else(is_effect_ref, 2, 0), # Dosage is 2 if effect allele is Ref, 0 otherwise
            hl.cond(
                mt_row.GT.is_het(),
                # Dosage for heterozygous calls
                (hl.if_else(is_effect_ref, 1, 0) + # Count 1 if effect allele is the reference
                 hl.if_else(is_effect_alt & mt_row.GT.contains(effect_allele_alt_index), 1, 0)), # Count 1 if effect allele is an ALT present in GT
                hl.cond(
                    mt_row.GT.is_hom_var(),
                    # Dosage for homozygous variant calls (only count 2 if hom var *for the effect allele*)
                    hl.if_else(is_effect_alt & mt_row.GT.ploidy == 2 & (mt_row.GT[0] == effect_allele_alt_index) & (mt_row.GT[1] == effect_allele_alt_index), 2, 0),
                    0 # Fallback for other non-diploid or unexpected GT states (should ideally be missing or handled)
                )
            )
        )
    )
    return dosage
# --- END CORRECT DOSAGE CALCULATION FUNCTION ---


def annotate_and_score(vds_filtered: hl.vds.VariantDataset, prs_ht: hl.Table, prs_id: str) -> hl.Table | None:
    """Annotates filtered VDS/MT with PRS weights and computes scores using correct dosage."""
    if vds_filtered is None or prs_ht is None:
        print(f"[{prs_id}] Skipping annotation/computation: Filtered VDS or PRS HT is None.")
        return None

    print(f"[{prs_id}] Annotating VDS with PRS info and computing scores...")
    try:
        # PRS HT is not empty before proceeding
        if prs_ht.count() == 0:
            print(f"[{prs_id}] ERROR: Input PRS HailTable for annotation is empty. Cannot annotate.")
            return None

        # Convert the variant_data part of the VDS to a MatrixTable for annotation/scoring
        mt = vds_filtered.variant_data
        n_variants_before_annot = mt.count_rows()
        if n_variants_before_annot == 0:
            print(f"[{prs_id}] WARNING: Filtered VDS variant_data has 0 rows before annotation. No scores can be computed.")
            # Return an empty table with the expected schema? Or None? Returning None is simpler.
            return None

        # Annotate rows with PRS info (effect allele, weight)
        print(f"[{prs_id}] Annotating {n_variants_before_annot} variants with PRS info...")
        mt = mt.annotate_rows(prs=prs_ht[mt.locus])

        # Filter to variants found in the PRS table (where annotation was successful)
        mt = mt.filter_rows(hl.is_defined(mt.prs) & hl.is_defined(mt.prs.weight) & hl.is_defined(mt.prs.effect_allele))

        rows_after_annot = mt.count_rows()
        print(f"[{prs_id}] Variants remaining after annotation (found in PRS table with needed info): {rows_after_annot}")
        if rows_after_annot == 0:
            print(f"[{prs_id}] WARNING: No variants overlapped between VDS and PRS table with valid info. No scores computed.")
            return None # Cannot compute scores

        # Calculate dosage using the correct function
        print(f"[{prs_id}] Calculating effect allele dosage...")
        mt = mt.unfilter_entries() # GT is available for all samples at these variants
        mt = mt.annotate_entries(effect_allele_count=calculate_effect_allele_dosage(mt)) # Pass the row context

        # Calculate per-variant contribution, handling missing dosage or weight
        # Treat missing contribution as zero for the sum aggregation.
        print(f"[{prs_id}] Calculating variant contributions...")
        mt = mt.annotate_entries(
            variant_contribution = hl.coalesce(mt.effect_allele_count * mt.prs.weight, 0.0)
        )

        # Aggregate scores per sample
        print(f"[{prs_id}] Aggregating scores per sample...")
        # Use select_cols for aggregation, more idiomatic than group_cols_by + aggregate
        prs_scores_ht = mt.select_cols(
            total_score=hl.agg.sum(mt.variant_contribution),
            # Count variants where calculation was possible (dosage defined)
            variant_count=hl.agg.count_where(hl.is_defined(mt.effect_allele_count))
        )

        # Normalize score (optional but good practice)
        prs_scores_ht = prs_scores_ht.annotate(
             normalized_score=hl.if_else(
                 prs_scores_ht.variant_count > 0,
                 prs_scores_ht.total_score / hl.float64(prs_scores_ht.variant_count),
                 hl.missing(hl.tfloat64) # Return missing if no variants contributed
             )
        )

        # Final checks and logging
        n_samples_scored = prs_scores_ht.count()
        print(f"[{prs_id}] PRS scores computed for {n_samples_scored} samples.")
        if n_samples_scored == 0:
            print(f"[{prs_id}] WARNING: Score computation resulted in 0 samples (might happen if filtered VDS had samples but none had overlapping PRS variants with genotypes).")
            # Return None as no useful scores were generated
            return None

        print(f"[{prs_id}] Result schema:")
        prs_scores_ht.describe()
        print(f"[{prs_id}] Example scores (first 5 samples):")
        prs_scores_ht.show(5)
        print("\n")
        return prs_scores_ht

    except Exception as e:
        print(f"[{prs_id}] ERROR during annotation or score computation: {e}")
        # Consider logging more details, e.g., mt.describe() at point of failure
        return None


def save_ht_and_csv(ht_scores: hl.Table | None, hail_table_gcs_path: str, score_csv_gcs_path: str, prs_id: str, fs) -> bool:
    """Write the resulting Hail Table and export its CSV representation to GCS."""
    if ht_scores is None:
        print(f"[{prs_id}] ERROR: Cannot save scores, input HailTable is None.")
        return False # Indicate failure

    # Check if the table is effectively empty before writing
    if ht_scores.count() == 0:
         print(f"[{prs_id}] WARNING: Input HailTable for saving scores is empty. Skipping write/export.")
         # Consider if this should be an error or just a skip. Let's treat as failure for pipeline task.
         return False # Indicate failure

    print(f"[{prs_id}] Preparing to save scores...")
    try:
        # parent directories exist on GCS
        parent_ht_dir = os.path.dirname(hail_table_gcs_path)
        parent_csv_dir = os.path.dirname(score_csv_gcs_path)
        if not gcs_path_exists(parent_ht_dir): fs.mkdirs(parent_ht_dir, exist_ok=True)
        if not gcs_path_exists(parent_csv_dir): fs.mkdirs(parent_csv_dir, exist_ok=True)

        # --- Save Hail Table ---
        print(f"[{prs_id}] Writing PRS results to Hail Table: {hail_table_gcs_path}")
        ht_scores.write(hail_table_gcs_path, overwrite=True)
        print(f"[{prs_id}] Hail Table write complete.")
        # Optional: Add verification read? Can be slow. hl.read_table(hail_table_gcs_path).count()

        # --- Export to CSV ---
        print(f"[{prs_id}] Exporting PRS scores to CSV: {score_csv_gcs_path}")
        # Select and rename columns for a clean CSV output
        ht_export = ht_scores.select(
            person_id=ht_scores.s, # Assuming the sample key is 's'
            prs_total_score=ht_scores.total_score,
            prs_variant_count=ht_scores.variant_count,
            prs_normalized_score=ht_scores.normalized_score
        )
        ht_export.export(score_csv_gcs_path, header=True, delimiter=',', missing='NA')
        print(f"[{prs_id}] PRS scores exported successfully to {score_csv_gcs_path}.\n")
        return True # Indicate success

    except Exception as e:
        print(f"[{prs_id}] ERROR during saving/exporting scores (HT: {hail_table_gcs_path}, CSV: {score_csv_gcs_path}): {e}")
        # Attempt cleanup of potentially failed writes
        delete_gcs_path(hail_table_gcs_path, recursive=True)
        delete_gcs_path(score_csv_gcs_path, recursive=False)
        return False # Indicate failure


# --- Main Execution ---

def main():
    parser = argparse.ArgumentParser(description="Process a single PRS model: download, prepare, filter, score, save.")
    parser.add_argument('--prs_id',                      required=True, help="PGS ID (e.g., PGS000123).")
    parser.add_argument('--prs_url',                     required=True, help="URL to the harmonized scoring file (txt.gz).")
    parser.add_argument('--base_cohort_vds_path',        required=True, help="GCS path to the prepared base cohort VDS.")
    parser.add_argument('--gcs_temp_dir',                required=True, help="GCS base directory for intermediate checkpoints (weights, intervals).")
    parser.add_argument('--gcs_hail_temp_dir',          required=True, help="GCS temporary directory for Hail operations.")
    parser.add_argument('--run_timestamp',                required=True, help="Pipeline run timestamp for logging.")
    parser.add_argument('--output_final_hail_table_gcs_path', required=True, help="GCS output path for the final scores Hail Table.")
    parser.add_argument('--output_final_score_csv_gcs_path',    required=True, help="GCS output path for the final scores CSV.")
    parser.add_argument('--google_billing_project',       required=True, help="Google Cloud Project ID for billing and GCS access.")
    parser.add_argument('--spark_configurations_json',  required=True, help="JSON string of Spark configurations for Hail initialization.")
    parser.add_argument(
        "--hail_cluster_mode", 
        choices=["local", "dataproc_yarn"], 
        default="local", 
        help="Hail execution mode: 'local' for local Spark, 'dataproc_yarn' for running on a Dataproc YARN cluster."
    )
    args = parser.parse_args()

    prs_id = args.prs_id

    print(f"--- Starting processing for PRS Model: {prs_id} ---")

    # Initialize Hail and GCS
    # Explicitly pass billing project to GCSFileSystem initialization
    fs = get_gcs_fs(project_id_for_billing=args.google_billing_project)
    init_hail(
        gcs_hail_temp_dir=args.gcs_hail_temp_dir,
        log_suffix=f"{args.run_timestamp}_{prs_id}", # Using run_timestamp and prs_id for specific log name
        spark_configurations_json_str=args.spark_configurations_json,
        cluster_mode=args.hail_cluster_mode # Pass the cluster mode to Hail initialization
    )
    # Set a default number of partitions for Hail operations.
    dynamic_partitions = max(200, hl.spark_context().defaultParallelism * 4)
    hl.utils.default_n_partitions(dynamic_partitions)

    # --- Define Intermediate Paths ---
    # These paths point to potentially reusable files across runs in the stable GCS_TEMP_DIR
    intermediate_weights_gcs_path = f"{args.gcs_temp_dir}/prs_weights/{prs_id}_prepared_weight_table.csv"
    intermediate_intervals_gcs_path = f"{args.gcs_temp_dir}/prs_intervals/{prs_id}_interval.tsv"

    # --- Step 1: Download and Prepare Weights/Intervals (with Checkpointing) ---
    weight_table_exists = gcs_path_exists(intermediate_weights_gcs_path)
    interval_file_exists = gcs_path_exists(intermediate_intervals_gcs_path)

    if weight_table_exists and interval_file_exists:
        print(f"[{prs_id}] [CHECKPOINT HIT] Found prepared weight table: {intermediate_weights_gcs_path}")
        print(f"[{prs_id}] [CHECKPOINT HIT] Found interval file: {intermediate_intervals_gcs_path}")
    else:
        print(f"[{prs_id}] [CHECKPOINT MISS] One or both intermediate files missing. Generating...")
        # Define local temporary file names
        local_harmonized_gz = f"{prs_id}_harmonized.txt.gz"
        local_harmonized_txt = f"{prs_id}_harmonized.txt"

        # Download
        print(f"[{prs_id}] Downloading PRS file...")
        local_downloaded_path = download_file(args.prs_url, local_harmonized_gz, local_harmonized_txt, prs_id) # From utils
        if local_downloaded_path is None:
            print(f"[{prs_id}] FATAL: Download failed.")
            sys.exit(1)

        score_df_prepared = None # Initialize before try block
        try:
            # Prepare DataFrame
            print(f"[{prs_id}] Preparing weight table DataFrame...")
            score_df_prepared = prepare_weight_table(local_downloaded_path, prs_id)
            if score_df_prepared is None:
                print(f"[{prs_id}] FATAL: Weight table preparation failed.")
                sys.exit(1) # Exit if preparation fails

        finally:
            # cleanup of the downloaded file regardless of preparation success/failure
            if local_downloaded_path and os.path.exists(local_downloaded_path):
                print(f"[{prs_id}] Cleaning up local downloaded file: {local_downloaded_path}")
                try:
                    os.remove(local_downloaded_path)
                except OSError as e:
                    print(f"[{prs_id}] WARNING: Could not remove local file {local_downloaded_path}: {e}")

        # Save/Upload Weights (if missing)
        if not weight_table_exists:
            print(f"[{prs_id}] Uploading prepared weights to checkpoint: {intermediate_weights_gcs_path}")
            if not save_to_gcs(score_df_prepared, prs_id, "weights", intermediate_weights_gcs_path, fs):
                print(f"[{prs_id}] FATAL: Failed to save/upload weights checkpoint.")
                sys.exit(1)
        else:
             print(f"[{prs_id}] Weight table already exists, skipping upload.")

        # Extract/Upload Intervals (if missing)
        if not interval_file_exists:
            print(f"[{prs_id}] Extracting and uploading intervals to checkpoint: {intermediate_intervals_gcs_path}")
            if not extract_intervals(score_df_prepared, prs_id, intermediate_intervals_gcs_path, fs):
                print(f"[{prs_id}] FATAL: Failed to extract/upload intervals checkpoint.")
                sys.exit(1)
        else:
             print(f"[{prs_id}] Interval file already exists, skipping upload.")

        # Clear DataFrame from memory if it was loaded
        del score_df_prepared
        print(f"[{prs_id}] Intermediate file preparation complete.")

    # --- Step 2: Load Base VDS ---
    print(f"[{prs_id}] Loading base cohort VDS: {args.base_cohort_vds_path}")
    try:
        base_vds = hl.vds.read_vds(args.base_cohort_vds_path)
        print(f"[{prs_id}] Base VDS loaded successfully.")
    except Exception as e:
        print(f"[{prs_id}] FATAL: Failed to load base VDS from {args.base_cohort_vds_path}: {e}")
        sys.exit(1)

    # --- Step 3: Filter VDS by Intervals ---
    vds_filtered = filter_vds_by_intervals(base_vds, intermediate_intervals_gcs_path, prs_id)
    if vds_filtered is None:
        print(f"[{prs_id}] FATAL: VDS filtering failed.")
        # Clean up base_vds from memory before exiting? Optional.
        del base_vds
        sys.exit(1)
    del base_vds # Free memory

    # --- Step 4: Load Prepared Weights into Hail ---
    prs_ht_for_annotation = import_weights_as_ht(intermediate_weights_gcs_path, prs_id)
    if prs_ht_for_annotation is None:
        print(f"[{prs_id}] FATAL: Failed to load weights into Hail.")
        del vds_filtered # Clean up memory
        sys.exit(1)

    # --- Step 5: Annotate and Compute Scores ---
    prs_scores_ht_result = annotate_and_score(vds_filtered, prs_ht_for_annotation, prs_id)
    # Clean up inputs immediately after use
    del vds_filtered
    del prs_ht_for_annotation

    if prs_scores_ht_result is None:
        print(f"[{prs_id}] FATAL: Score computation failed or produced no results.")
        # No scores table to clean up.
        sys.exit(1) # Treat no scores as a failure for this task

    # --- Step 6: Save Final Scores (Hail Table and CSV) ---
    print(f"[{prs_id}] Saving final scores...")
    save_success = save_ht_and_csv(
        prs_scores_ht_result,
        args.output_final_hail_table_gcs_path,
        args.output_final_score_csv_gcs_path,
        prs_id,
        fs
    )
    del prs_scores_ht_result # Clean up memory

    if not save_success:
        print(f"[{prs_id}] FATAL: Failed to save final scores.")
        sys.exit(1)

    print(f"--- Successfully completed processing for PRS Model: {prs_id} ---")
    print(f"Final Scores CSV written to: {args.output_final_score_csv_gcs_path}")
    print(f"Final Scores Hail Table written to: {args.output_final_hail_table_gcs_path}")


if __name__ == "__main__":
    main()
