import argparse
import datetime
import os
import sys
import shutil
import pandas as pd
import hail as hl
from common_utils import (
    init_hail, gcs_path_exists, hail_path_exists, delete_gcs_path,
    get_gcs_fs, download_file
)

# --- PRS Calculation Functions ---

def preview_harmonized_file(harmonized_file, prs_id):
    if not harmonized_file or not os.path.exists(harmonized_file):
        print(f"[{prs_id}] Cannot preview file, {harmonized_file} does not exist or path is invalid.")
        return
    print(f"[{prs_id}] Previewing the first 10 lines of {harmonized_file}:")
    try:
        with open(harmonized_file, 'r') as f:
            for i, line in enumerate(f):
                if i >= 10: break
                print(line.strip())
    except Exception as e:
        print(f"[{prs_id}] Error previewing file: {e}")
    print(f"[{prs_id}] Preview complete.\n")

def prepare_prs_weight_table_df(harmonized_file, prs_id):
    if not harmonized_file or not os.path.exists(harmonized_file):
        print(f"[{prs_id}] Skipping weight table preparation: Input file missing '{harmonized_file}'.")
        return None
    print(f"[{prs_id}] Reading harmonized PRS weight file: {harmonized_file}")
    try:
        score_df = pd.read_csv(harmonized_file, sep='\t', comment='#', low_memory=False)
        print(f"[{prs_id}] PRS weight data loaded. Initial variant count: {len(score_df)}")

        col_map = {
            'hm_chr': 'chr', 'chr': 'chr', 'chromosome': 'chr',
            'hm_pos': 'bp', 'pos': 'bp', 'base_pair_location': 'bp', 'position': 'bp', 'BP': 'bp',
            'effect_allele': 'effect_allele', 'allele': 'effect_allele', 'A1': 'effect_allele', 'hm_effect_allele': 'effect_allele',
            'other_allele': 'other_allele', 'A2': 'other_allele', 'hm_other_allele': 'other_allele', 'reference_allele': 'other_allele',
            'effect_weight': 'weight', 'beta': 'weight', 'OR': 'weight', 'effect': 'weight', 'hm_beta': 'weight'
        }
        rename_dict = {}
        found_cols_targets = [] # Track which target columns we've found a source for
        needed_cols_targets = ['chr', 'bp', 'effect_allele', 'weight'] # Target names

        # Greedily find mappings for needed columns
        for target_name in needed_cols_targets:
            for potential_source_name, mapped_target_name in col_map.items():
                if mapped_target_name == target_name and potential_source_name in score_df.columns:
                    if target_name not in found_cols_targets: # Only map once per target
                         rename_dict[potential_source_name] = target_name
                         found_cols_targets.append(target_name)
                         break # Found mapping for this target_name

        print(f"[{prs_id}] Renaming columns: {rename_dict}")
        score_df = score_df.rename(columns=rename_dict)

        missing_essential = [col for col in needed_cols_targets if col not in score_df.columns]
        if missing_essential:
            raise RuntimeError(f"Missing essential columns after renaming: {missing_essential}. Available columns: {score_df.columns.tolist()}")

        print(f"[{prs_id}] Selecting and cleaning required columns: {needed_cols_targets}")
        score_df = score_df[needed_cols_targets].copy()
        score_df['chr'] = score_df['chr'].astype(str)
        score_df['effect_allele'] = score_df['effect_allele'].astype(str)

        score_df['bp'] = pd.to_numeric(score_df['bp'], errors='coerce')
        score_df['weight'] = pd.to_numeric(score_df['weight'], errors='coerce')

        initial_len = len(score_df)
        essential_for_dropna = ['bp', 'weight', 'effect_allele', 'chr']
        score_df = score_df.dropna(subset=essential_for_dropna)
        dropped_count = initial_len - len(score_df)
        if dropped_count > 0:
            print(f"[{prs_id}] Warning: Dropped {dropped_count} rows during cleaning (missing/invalid pos, weight, allele, chr).")

        if score_df.empty:
            print(f"[{prs_id}] ERROR: DataFrame became empty after cleaning. Check input file format/content.")
            return None

        score_df['bp'] = score_df['bp'].astype(int)

        print(f"[{prs_id}] Adding Hail-compatible columns: contig, position (integer)")
        score_df = score_df[score_df['chr'].notna() & (score_df['chr'] != '')]
        score_df['contig'] = score_df['chr'].apply(lambda c: str(c) if str(c).startswith('chr') else f'chr{str(c)}')
        score_df['position'] = score_df['bp']

        final_count = len(score_df)
        print(f"[{prs_id}] PRS weight table prepared. Final variant count: {final_count}.\n")
        if final_count == 0:
            print(f"[{prs_id}] ERROR: Prepared weight table has 0 variants after processing.")
            return None
        return score_df
    except Exception as e:
        print(f"[{prs_id}] ERROR preparing weight table: {e}")
        if 'score_df' in locals() and isinstance(score_df, pd.DataFrame): print(f"Columns at error: {score_df.columns.tolist()}")
        return None


def save_and_upload_dataframe_to_gcs(df, prs_id, local_temp_filename_base, target_gcs_path, fs_gcs):
    if df is None or df.empty:
        print(f"[{prs_id}] Skipping save/upload for {local_temp_filename_base}: DataFrame is None or empty.")
        return False

    # Use a unique local temporary file name, possibly in a task-specific directory
    # For simplicity, using current dir, assuming Nextflow task dirs are isolated.
    local_temp_file_path = f"{local_temp_filename_base}_{prs_id}_temp.csv"

    print(f"[{prs_id}] Storing processed {local_temp_filename_base} locally as: {local_temp_file_path}")
    try:
        df.to_csv(local_temp_file_path, index=False, na_rep='NA')
    except Exception as e:
        print(f"[{prs_id}] ERROR saving {local_temp_filename_base} locally to {local_temp_file_path}: {e}")
        return False

    print(f"[{prs_id}] Copying {local_temp_filename_base} to GCS: {target_gcs_path}")
    try:
        target_dir = os.path.dirname(target_gcs_path)
        if not gcs_path_exists(target_dir): # common_utils
            print(f"[{prs_id}] Creating GCS directory: {target_dir}")
            fs_gcs.mkdirs(target_dir, exist_ok=True)

        print(f"[{prs_id}] Uploading {local_temp_file_path} to {target_gcs_path} using gcsfs.")
        fs_gcs.put(local_temp_file_path, target_gcs_path)
        print(f"[{prs_id}] Upload to GCS completed: {target_gcs_path}")
        os.remove(local_temp_file_path)
        print(f"[{prs_id}] Removed local temporary file: {local_temp_file_path}")
        return True
    except Exception as e:
        print(f"[{prs_id}] ERROR uploading {local_temp_file_path} to GCS path {target_gcs_path}: {e}")
        if os.path.exists(local_temp_file_path):
            try: os.remove(local_temp_file_path)
            except OSError: pass
        delete_gcs_path(target_gcs_path, recursive=False) # common_utils
        return False

def extract_save_upload_intervals_df(score_df, prs_id, target_gcs_path, fs_gcs):
    if score_df is None or score_df.empty:
        print(f"[{prs_id}] Skipping interval extraction: score_df is None or empty.")
        return False
    print(f"[{prs_id}] Extracting variant intervals for filtering VDS...")
    required_cols = ['contig', 'position']
    if not all(col in score_df.columns for col in required_cols):
        print(f"[{prs_id}] ERROR: Missing required columns ({required_cols}) in score_df for interval extraction.")
        return False
    if not pd.api.types.is_integer_dtype(score_df['position']):
        print(f"[{prs_id}] ERROR: 'position' column is not integer type ({score_df['position'].dtype}) before interval extraction.")
        return False

    try:
        intervals_df = score_df[required_cols].copy()
        intervals_df['end'] = intervals_df['position']
        intervals_df = intervals_df[intervals_df['contig'].notna() & (intervals_df['contig'] != '')]
        # Hail import_locus_intervals expects: contig<tab>start<tab>end
        # start and end are 1-based, inclusive. For single points, start=end=position.
        intervals_df = intervals_df[['contig', 'position', 'end']].drop_duplicates()
        # Ensure 'position' and 'end' are integers
        intervals_df['position'] = intervals_df['position'].astype(int)
        intervals_df['end'] = intervals_df['end'].astype(int)

    except Exception as e:
        print(f"[{prs_id}] ERROR during DataFrame manipulation for intervals: {e}")
        return False

    if intervals_df.empty:
        print(f"[{prs_id}] ERROR: No intervals generated after processing and deduplication.")
        return False

    # Save locally then upload
    local_temp_interval_file = f"{prs_id}_interval_temp.tsv" # In current Nextflow work_dir
    print(f"[{prs_id}] Saving {len(intervals_df)} unique intervals locally to: {local_temp_interval_file}")
    try:
        # Write as tab-separated, no header, for hl.import_locus_intervals
        intervals_df.to_csv(local_temp_interval_file, header=False, index=False, sep="\t")
    except Exception as e:
        print(f"[{prs_id}] ERROR saving intervals locally to {local_temp_interval_file}: {e}")
        return False

    print(f"[{prs_id}] Uploading interval list to GCS: {target_gcs_path}")
    try:
        target_dir = os.path.dirname(target_gcs_path)
        if not gcs_path_exists(target_dir): # common_utils
            print(f"[{prs_id}] Creating GCS directory: {target_dir}")
            fs_gcs.mkdirs(target_dir, exist_ok=True)

        fs_gcs.put(local_temp_interval_file, target_gcs_path)
        print(f"[{prs_id}] Intervals uploaded successfully to {target_gcs_path}")
        os.remove(local_temp_interval_file)
        print(f"[{prs_id}] Removed local temporary interval file: {local_temp_interval_file}")
        return True
    except Exception as e:
        print(f"[{prs_id}] ERROR uploading intervals {local_temp_interval_file} to GCS path {target_gcs_path}: {e}")
        if os.path.exists(local_temp_interval_file): os.remove(local_temp_interval_file)
        delete_gcs_path(target_gcs_path, recursive=False) # common_utils
        return False


def load_prepared_prs_table_for_hail_annotation(full_weight_gcs_path, prs_id, fs_gcs):
    print(f"[{prs_id}] Importing PRS weights from GCS for Hail annotation: {full_weight_gcs_path}")
    if not gcs_path_exists(full_weight_gcs_path): # common_utils
        print(f"[{prs_id}] ERROR: Weight file not found on GCS at {full_weight_gcs_path}")
        return None
    try:
        with fs_gcs.open(full_weight_gcs_path, 'r') as f: # Use gcsfs to peek
            header_line = f.readline()
            if not header_line: raise ValueError("PRS weight file is empty.")
            header = header_line.strip().split(',')
        print(f"[{prs_id}] File header: {header}")
        essential_cols = ["weight", "contig", "position", "effect_allele", "bp"]
        missing = [col for col in essential_cols if col not in header]
        if missing:
            raise ValueError(f"PRS weight table at {full_weight_gcs_path} is missing required columns: {missing}. Header: {header}")

        prs_ht = hl.import_table(full_weight_gcs_path,
                                 delimiter=',',
                                 types={'weight': hl.tfloat64, 'position': hl.tint32, 'bp': hl.tint32},
                                 missing='NA', comment='#', impute=False) # Set impute=False if types are strict

        if prs_ht.position.dtype != hl.tint32 or prs_ht.bp.dtype != hl.tint32:
            raise TypeError(f"[{prs_id}] Hail imported position/bp with wrong type: pos={prs_ht.position.dtype}, bp={prs_ht.bp.dtype}. Expected int32.")
        if prs_ht.weight.dtype != hl.tfloat64:
            print(f"[{prs_id}] WARNING: Hail imported weight column with type {prs_ht.weight.dtype}, expected float64. Attempting cast.")
            prs_ht = prs_ht.annotate(weight=hl.float64(prs_ht.weight))
        if prs_ht.contig.dtype != hl.tstr: prs_ht = prs_ht.annotate(contig=hl.str(prs_ht.contig))
        if prs_ht.effect_allele.dtype != hl.tstr: prs_ht = prs_ht.annotate(effect_allele=hl.str(prs_ht.effect_allele))


        n_rows_imported = prs_ht.count()
        if n_rows_imported == 0:
            print(f"[{prs_id}] ERROR: Imported PRS HailTable has 0 rows.")
            return None

        prs_ht = prs_ht.filter(hl.is_defined(prs_ht.contig) & hl.is_defined(prs_ht.position) & (prs_ht.contig != ""))
        prs_ht = prs_ht.annotate(locus=hl.locus(prs_ht.contig, prs_ht.position, reference_genome='GRCh38'))
        prs_ht = prs_ht.filter(hl.is_defined(prs_ht.locus))

        n_rows_locus = prs_ht.count()
        if n_rows_locus == 0:
            print(f"[{prs_id}] ERROR: No rows remaining after locus creation/filtering.")
            return None
        elif n_rows_locus < n_rows_imported:
            print(f"[{prs_id}] WARNING: Dropped {n_rows_imported - n_rows_locus} rows during locus creation/filtering.")

        prs_ht = prs_ht.key_by('locus')
        prs_ht = prs_ht.select('weight', 'effect_allele')
        print(f"[{prs_id}] PRS HT ready for annotation ({prs_ht.count()} variants). Schema: {prs_ht.row.dtype}\n")
        return prs_ht
    except Exception as e:
        print(f"[{prs_id}] ERROR loading prepared PRS table into Hail: {e}")
        return None

def calculate_effect_allele_dosage_hail(mt_row):
    eff_allele = hl.str(mt_row.prs_variant_info.effect_allele)
    ref_allele = mt_row.alleles[0]
    # Find if effect allele is one of the alt alleles and get its index (1-based for GT)
    alt_indices = hl.range(hl.len(mt_row.alleles) - 1).map(
        lambda i: hl.struct(idx=i + 1, allele=mt_row.alleles[i+1])
    )
    effect_allele_alt_struct = hl.find(lambda x: x.allele == eff_allele, alt_indices)
    effect_allele_alt_idx = hl.if_else(hl.is_defined(effect_allele_alt_struct),
                                       effect_allele_alt_struct.idx,
                                       hl.missing(hl.tint32))

    is_effect_ref = (ref_allele == eff_allele)
    is_effect_alt_and_defined = hl.is_defined(effect_allele_alt_idx)

    dosage = hl.cond(
        hl.is_missing(mt_row.GT),
        hl.missing(hl.tint32),
        hl.cond(
            mt_row.GT.is_hom_ref(),
            hl.if_else(is_effect_ref, 2, 0),
            hl.cond(
                mt_row.GT.is_het(),
                (hl.if_else(is_effect_ref, 1, 0) +
                 hl.if_else(is_effect_alt_and_defined & mt_row.GT.contains(effect_allele_alt_idx), 1, 0)),
                hl.cond(
                    mt_row.GT.is_hom_var(),
                    hl.if_else(is_effect_alt_and_defined & (mt_row.GT[0] == effect_allele_alt_idx) & (mt_row.GT[1] == effect_allele_alt_idx), 2, 0),
                    0 # Fallback for other non-diploid GT or unexpected states (e.g. haploid)
                )
            )
        )
    )
    return dosage


def filter_vds_by_hail_intervals(vds, interval_gcs_path, prs_id):
    if vds is None:
        print(f"[{prs_id}] ERROR: VDS is None, cannot filter by intervals.")
        return None
    if interval_gcs_path is None:
        print(f"[{prs_id}] WARNING: Interval GCS path is None. Skipping VDS interval filtering.")
        return vds
    print(f"[{prs_id}] Filtering VDS by PRS intervals from: {interval_gcs_path}")
    if not gcs_path_exists(interval_gcs_path): # common_utils
        print(f"[{prs_id}] ERROR: Interval file not found at {interval_gcs_path}. Cannot filter.")
        return None
    try:
        print(f"[{prs_id}] Importing locus intervals from {interval_gcs_path}...")
        locus_intervals = hl.import_locus_intervals(
            interval_gcs_path, reference_genome='GRCh38', skip_invalid_intervals=True
        )
        interval_count = locus_intervals.count()
        if interval_count == 0:
            print(f"[{prs_id}] WARNING: No valid intervals loaded from {interval_gcs_path}. Filtering will likely remove all variants.")
        else:
            print(f"[{prs_id}] Loaded {interval_count} intervals for filtering.")

        print(f"[{prs_id}] Applying interval filter to VDS...")
        vds_filtered_variant_data = hl.filter_intervals(vds.variant_data, locus_intervals)
        vds_filtered_ref_data = hl.filter_intervals(vds.reference_data, locus_intervals)
        vds_filtered = hl.vds.VariantDataset(vds_filtered_ref_data, vds_filtered_variant_data)

        n_variants_after = vds_filtered.variant_data.count_rows()
        print(f"[{prs_id}] VDS filtered to PRS intervals. Variants remaining: {n_variants_after}\n")
        if n_variants_after == 0:
            print(f"[{prs_id}] WARNING: No variants remained after interval filtering.")
        return vds_filtered
    except Exception as e:
        print(f"[{prs_id}] ERROR filtering VDS by intervals: {e}")
        return None

def annotate_and_compute_prs_scores_hail(vds_filtered, prs_ht, prs_id):
    if vds_filtered is None or prs_ht is None:
        print(f"[{prs_id}] Skipping annotation/computation: Filtered VDS or PRS HT is None.")
        return None
    print(f"[{prs_id}] Annotating VDS with PRS info and computing scores...")
    try:
        if prs_ht.count() == 0:
            print(f"[{prs_id}] ERROR: Input PRS HailTable is empty. Cannot annotate.")
            return None

        mt = vds_filtered.variant_data # Operate on variant_data MatrixTable
        n_variants_before_annot = mt.count_rows()
        if n_variants_before_annot == 0:
            print(f"[{prs_id}] WARNING: VDS variant_data has 0 rows before annotation. No scores can be computed.")
            return None

        print(f"[{prs_id}] Annotating {n_variants_before_annot} variants with PRS info...")
        mt = mt.annotate_rows(prs_variant_info=prs_ht[mt.locus])
        mt = mt.filter_rows(hl.is_defined(mt.prs_variant_info)) # Keep only variants in PRS

        rows_after_annot_filter = mt.count_rows()
        print(f"[{prs_id}] Variants remaining after annotation and filtering to PRS variants: {rows_after_annot_filter}")
        if rows_after_annot_filter == 0:
            print(f"[{prs_id}] WARNING: No variants overlapped between VDS and PRS table. No scores computed.")
            return None

        print(f"[{prs_id}] Calculating effect allele dosage...")
        mt = mt.unfilter_entries() # Ensure all entries are considered for dosage
        mt = mt.annotate_entries(effect_allele_count=calculate_effect_allele_dosage_hail(mt))

        print(f"[{prs_id}] Calculating variant contributions...")
        mt = mt.annotate_entries(
            variant_contribution=hl.if_else(
                hl.is_defined(mt.effect_allele_count) & hl.is_defined(mt.prs_variant_info.weight),
                mt.effect_allele_count * mt.prs_variant_info.weight,
                0.0 # Treat missing contribution as zero for summation
            )
        )

        print(f"[{prs_id}] Aggregating scores per sample...")
        prs_scores_ht = mt.select_cols(
            total_score=hl.agg.sum(mt.variant_contribution),
            variant_count=hl.agg.count_where(
                hl.is_defined(mt.effect_allele_count) & hl.is_defined(mt.prs_variant_info.weight) & (mt.variant_contribution != 0) # Count variants truly contributing
            )
        )
        prs_scores_ht = prs_scores_ht.annotate(
            normalized_score=hl.cond(
                prs_scores_ht.variant_count > 0,
                prs_scores_ht.total_score / hl.float64(prs_scores_ht.variant_count),
                hl.missing(hl.tfloat64)
            )
        )

        n_samples_scored = prs_scores_ht.count()
        print(f"[{prs_id}] PRS scores computed for {n_samples_scored} samples.")
        if n_samples_scored == 0:
            print(f"[{prs_id}] WARNING: Score computation resulted in 0 samples with scores.")
            return None
        print(f"[{prs_id}] Result schema:"); prs_scores_ht.describe(); prs_scores_ht.show(5); print("\n")
        return prs_scores_ht
    except Exception as e:
        print(f"[{prs_id}] ERROR during annotation or score computation: {e}")
        return None

def save_final_prs_scores_hail(prs_scores_ht, prs_id, output_hail_table_gcs_path, output_score_csv_gcs_path, fs_gcs):
    if prs_scores_ht is None:
        print(f"[{prs_id}] Skipping saving scores: prs_scores_ht is None.")
        return False
    if prs_scores_ht.count() == 0:
        print(f"[{prs_id}] Skipping saving scores: Input HailTable is empty.")
        return False

    # Save Hail Table
    print(f"[{prs_id}] Writing PRS results to Hail Table: {output_hail_table_gcs_path}")
    try:
        prs_scores_ht.write(output_hail_table_gcs_path, overwrite=True)
        print(f"[{prs_id}] Hail Table write complete.")
        hl.read_table(output_hail_table_gcs_path).count() # Verify
        print(f"[{prs_id}] Written Hail Table verified.")
    except Exception as e:
        print(f"[{prs_id}] ERROR writing or verifying Hail Table to {output_hail_table_gcs_path}: {e}")
        delete_gcs_path(output_hail_table_gcs_path, recursive=True) # common_utils
        return False

    # Export to CSV
    print(f"[{prs_id}] Exporting PRS scores to CSV: {output_score_csv_gcs_path}")
    try:
        # Ensure parent GCS directory for CSV exists
        csv_dir_gcs = os.path.dirname(output_score_csv_gcs_path)
        if not gcs_path_exists(csv_dir_gcs): # common_utils
            print(f"Creating GCS directory for CSV: {csv_dir_gcs}")
            fs_gcs.mkdirs(csv_dir_gcs, exist_ok=True)

        # Read from the just-written Hail table for export
        final_ht_to_export = hl.read_table(output_hail_table_gcs_path)
        final_ht_to_export = final_ht_to_export.select(
            person_id=final_ht_to_export.s,
            prs_total_score=hl.coalesce(final_ht_to_export.total_score, hl.missing(hl.tfloat64)),
            prs_variant_count=hl.coalesce(final_ht_to_export.variant_count, hl.missing(hl.tint64)),
            prs_normalized_score=hl.coalesce(final_ht_to_export.normalized_score, hl.missing(hl.tfloat64))
        )
        final_ht_to_export.export(output_score_csv_gcs_path, header=True, delimiter=',', missing='NA')
        print(f"[{prs_id}] PRS scores exported successfully to {output_score_csv_gcs_path}.\n")
        return True
    except Exception as e:
        print(f"[{prs_id}] ERROR exporting scores to CSV {output_score_csv_gcs_path}: {e}")
        delete_gcs_path(output_score_csv_gcs_path, recursive=False) # common_utils
        return False


def main():
    parser = argparse.ArgumentParser(description="Process a single PRS model.")
    parser.add_argument("--prs_id", required=True, help="PGS ID of the model.")
    parser.add_argument("--prs_url", required=True, help="URL to the harmonized GRCh38 scoring file.")
    parser.add_argument("--base_cohort_vds_path", required=True, help="GCS path to the prepared base cohort VDS.")
    parser.add_argument("--gcs_temp_dir", required=True, help="GCS temp directory for intermediate files (weights, intervals).")
    parser.add_argument("--gcs_output_score_dir", required=True, help="GCS base directory for final run-specific scores and outputs.")
    parser.add_argument("--gcs_hail_temp_dir", required=True, help="GCS temporary directory for Hail operations.")
    parser.add_argument("--run_timestamp", required=True, help="Current run timestamp.")
    # Output paths that this script will generate
    parser.add_argument("--output_final_score_csv_gcs_path", required=True, help="GCS path where the final score CSV for this model will be written.")
    parser.add_argument("--output_final_hail_table_gcs_path", required=True, help="GCS path where the final Hail score table for this model will be written.")

    args = parser.parse_args()

    prs_id = args.prs_id
    print(f"\n------- Processing PRS Model: {prs_id} -------")

    # Initialize GCS FS and Hail
    fs = get_gcs_fs()
    init_hail(args.gcs_hail_temp_dir, args.run_timestamp)

    # Define paths for intermediate files (in GCS_TEMP_DIR for potential reuse across runs)
    temp_weights_gcs_path = f"{args.gcs_temp_dir}/prs_weights/{prs_id}_weight_table.csv"
    temp_intervals_gcs_path = f"{args.gcs_temp_dir}/prs_intervals/{prs_id}_interval.tsv"

    # Check if final output already exists (e.g. from a resumed Nextflow run)
    # The main.nf script should handle this top-level checkpointing logic more robustly.
    # Here, we proceed assuming this script is called when outputs are needed.

    # --- Preparation Steps if intermediates are missing ---
    score_df_prepared_for_intervals = None # To hold df if needed for both weights and intervals

    if not gcs_path_exists(temp_weights_gcs_path) or not gcs_path_exists(temp_intervals_gcs_path):
        print(f"[{prs_id}] One or more intermediate files missing. Generating...")
        local_harmonized_gz = f"{prs_id}_hmPOS_GRCh38.txt.gz" # In current work dir
        local_harmonized_txt = local_harmonized_gz.replace('.gz', '')

        downloaded_file_path = download_file(args.prs_url, local_harmonized_gz, local_harmonized_txt, prs_id) # common_utils
        if downloaded_file_path is None:
            print(f"[{prs_id}] FATAL: Skipping model due to download failure."); sys.exit(1)
        preview_harmonized_file(downloaded_file_path, prs_id)

        score_df_prepared = prepare_prs_weight_table_df(downloaded_file_path, prs_id)
        if os.path.exists(downloaded_file_path): os.remove(downloaded_file_path)

        if score_df_prepared is None:
            print(f"[{prs_id}] FATAL: Skipping model due to weight table preparation failure."); sys.exit(1)
        score_df_prepared_for_intervals = score_df_prepared.copy() # Keep a copy for interval generation

        if not gcs_path_exists(temp_weights_gcs_path):
            print(f"[{prs_id}] Uploading prepared weights to TEMP: {temp_weights_gcs_path}")
            if not save_and_upload_dataframe_to_gcs(score_df_prepared, prs_id, "weights", temp_weights_gcs_path, fs):
                print(f"[{prs_id}] FATAL: Failed to save/upload weights to {temp_weights_gcs_path}."); sys.exit(1)
        del score_df_prepared # Free memory

        if not gcs_path_exists(temp_intervals_gcs_path):
            print(f"[{prs_id}] Uploading intervals to TEMP: {temp_intervals_gcs_path}")
            if not extract_save_upload_intervals_df(score_df_prepared_for_intervals, prs_id, temp_intervals_gcs_path, fs):
                print(f"[{prs_id}] FATAL: Failed to extract/save/upload intervals to {temp_intervals_gcs_path}."); sys.exit(1)
        if score_df_prepared_for_intervals is not None: del score_df_prepared_for_intervals

    else:
        print(f"[{prs_id}] Found prepared weight table at: {temp_weights_gcs_path}")
        print(f"[{prs_id}] Found interval file at: {temp_intervals_gcs_path}")

    # --- Hail Steps ---
    print(f"[{prs_id}] Loading base VDS from: {args.base_cohort_vds_path}")
    base_vds = hl.vds.read_vds(args.base_cohort_vds_path) # This must exist

    print(f"[{prs_id}] Filtering base VDS using intervals from {temp_intervals_gcs_path}")
    vds_filtered_for_model = filter_vds_by_hail_intervals(base_vds, temp_intervals_gcs_path, prs_id)
    del base_vds
    if vds_filtered_for_model is None:
        print(f"[{prs_id}] FATAL: Skipping model due to VDS filtering failure."); sys.exit(1)

    print(f"[{prs_id}] Loading prepared weights for annotation from {temp_weights_gcs_path}")
    prs_ht_for_annotation = load_prepared_prs_table_for_hail_annotation(temp_weights_gcs_path, prs_id, fs)
    if prs_ht_for_annotation is None:
        print(f"[{prs_id}] FATAL: Skipping model due to failure loading weights into Hail."); sys.exit(1)

    print(f"[{prs_id}] Annotating VDS and computing scores...")
    prs_scores_ht_result = annotate_and_compute_prs_scores_hail(vds_filtered_for_model, prs_ht_for_annotation, prs_id)
    del vds_filtered_for_model, prs_ht_for_annotation
    if prs_scores_ht_result is None:
        print(f"[{prs_id}] FATAL: Skipping model due to score computation failure."); sys.exit(1)

    print(f"[{prs_id}] Saving final scores...")
    # Output paths are now passed as arguments
    success = save_final_prs_scores_hail(prs_scores_ht_result, prs_id,
                                         args.output_final_hail_table_gcs_path,
                                         args.output_final_score_csv_gcs_path, fs)
    del prs_scores_ht_result
    if not success:
        print(f"[{prs_id}] FATAL: Failed to save final scores."); sys.exit(1)

    print(f"[{prs_id}] Successfully calculated and saved final scores to {args.output_final_score_csv_gcs_path} and {args.output_final_hail_table_gcs_path}.")
    print(f"------- Finished processing {prs_id} -------")

if __name__ == "__main__":
    main()
