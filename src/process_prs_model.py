"""
process_prs_model.py

Nextflow task script to:
  1. Download and prepare PRS weight table
  2. Extract genomic intervals
  3. Filter a base Hail VDS by those intervals
  4. Annotate and compute per-sample PRS scores in Hail
  5. Save final Hail Table & CSV to GCS

Uses the following helper functions exist in common_utils:
  - init_hail(tmp_dir: str, run_timestamp: str)
  - get_gcs_fs() -> gcsfs.GCSFileSystem
  - gcs_path_exists(path: str) -> bool
  - delete_gcs_path(path: str, recursive: bool = True) -> bool
  - download_file(url: str, gz_path: str, out_path: str, prs_id: str) -> str | None
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

def prepare_weight_table(harmonized_path: str, prs_id: str) -> pd.DataFrame | None:
    """Read, clean, and prepare PRS weight table as a DataFrame."""
    if not os.path.exists(harmonized_path):
        print(f"[{prs_id}] ERROR: Missing harmonized file at {harmonized_path}")
        return None

    df = pd.read_csv(harmonized_path, sep='\t', comment='#', low_memory=False)
    print(f"[{prs_id}] Loaded {len(df)} rows from weight file")

    col_map = {
        'hm_chr': 'chr', 'chromosome': 'chr',
        'hm_pos': 'bp', 'position': 'bp', 'BP': 'bp',
        'effect_allele': 'effect_allele', 'A1': 'effect_allele',
        'other_allele': 'other_allele', 'A2': 'other_allele',
        'effect_weight': 'weight', 'beta': 'weight', 'OR': 'weight'
    }
    # Rename columns by priority
    rename = {}
    for src, tgt in col_map.items():
        if tgt not in rename.values() and src in df.columns:
            rename[src] = tgt
    df = df.rename(columns=rename)
    needed = ['chr', 'bp', 'effect_allele', 'weight']
    missing = [c for c in needed if c not in df.columns]
    if missing:
        print(f"[{prs_id}] ERROR: Missing columns {missing} after rename")
        return None

    df = df[needed].copy()
    df['bp']     = pd.to_numeric(df['bp'], errors='coerce')
    df['weight'] = pd.to_numeric(df['weight'], errors='coerce')
    df = df.dropna(subset=needed)
    if df.empty:
        print(f"[{prs_id}] ERROR: No valid rows after cleaning")
        return None

    df['bp'] = df['bp'].astype(int)
    df['contig']  = df['chr'].astype(str).apply(lambda c: c if c.startswith('chr') else f'chr{c}')
    df['position'] = df['bp']

    print(f"[{prs_id}] Prepared weight table with {len(df)} variants")
    return df[['contig','position','effect_allele','weight','bp']]

def save_to_gcs(df: pd.DataFrame, prs_id: str, local_name: str, gcs_path: str, fs) -> bool:
    """Save a DataFrame locally then upload to GCS."""
    tmp = f"{prs_id}_{local_name}.csv"
    df.to_csv(tmp, index=False)
    print(f"[{prs_id}] Saved local {tmp}, uploading to {gcs_path}")
    parent = os.path.dirname(gcs_path)
    if not gcs_path_exists(parent):
        fs.mkdirs(parent, exist_ok=True)
    try:
        fs.put(tmp, gcs_path)
    except Exception as e:
        print(f"[{prs_id}] ERROR uploading {tmp} to {gcs_path}: {e}")
        os.remove(tmp)
        delete_gcs_path(gcs_path, recursive=False)
        return False
    os.remove(tmp)
    return True

def extract_intervals(df: pd.DataFrame, prs_id: str, gcs_path: str, fs) -> bool:
    """Extract contig:start-end TSV and upload to GCS."""
    iv = df[['contig','position']].copy()
    iv['end'] = iv['position']
    iv = iv.drop_duplicates()
    if iv.empty:
        print(f"[{prs_id}] ERROR: No intervals to extract")
        return False

    tmp = f"{prs_id}_intervals.tsv"
    iv.to_csv(tmp, sep='\t', index=False, header=False)
    print(f"[{prs_id}] Saved local intervals TSV, uploading to {gcs_path}")
    parent = os.path.dirname(gcs_path)
    if not gcs_path_exists(parent):
        fs.mkdirs(parent, exist_ok=True)
    try:
        fs.put(tmp, gcs_path)
    except Exception as e:
        print(f"[{prs_id}] ERROR uploading intervals to {gcs_path}: {e}")
        os.remove(tmp)
        delete_gcs_path(gcs_path, recursive=False)
        return False
    os.remove(tmp)
    return True

def import_weights_as_ht(gcs_csv: str, prs_id: str) -> hl.Table | None:
    """Import the prepared CSV into a Hail Table keyed by locus."""
    if not gcs_path_exists(gcs_csv):
        print(f"[{prs_id}] ERROR: Weight CSV not found at {gcs_csv}")
        return None

    ht = hl.import_table(
        gcs_csv,
        delimiter=',',
        types={'contig': hl.tstr, 'position': hl.tint32, 'weight': hl.tfloat64},
        missing='NA',
        comment='#'
    )
    ht = ht.annotate(locus=hl.locus(ht.contig, ht.position, 'GRCh38')).key_by('locus')
    print(f"[{prs_id}] Imported HailTable with {ht.count()} variants")
    return ht.select('weight','effect_allele')

def filter_vds_by_intervals(vds: hl.VariantDataset, gcs_intervals: str, prs_id: str) -> hl.VariantDataset | None:
    """Filter the base VDS to only those intervals."""
    if not gcs_path_exists(gcs_intervals):
        print(f"[{prs_id}] ERROR: Interval file missing at {gcs_intervals}")
        return None
    intervals = hl.import_locus_intervals(gcs_intervals, 'GRCh38', skip_invalid_intervals=True)
    print(f"[{prs_id}] Loaded {intervals.count()} intervals")
    vds_f = vds.filter_intervals(intervals)
    print(f"[{prs_id}] VDS filtered: {vds_f.variant_data.count_rows()} variants remain")
    return vds_f

def annotate_and_score(vds: hl.VariantDataset, ht: hl.Table, prs_id: str) -> hl.Table | None:
    """Annotate VDS with PRS Table and compute per-sample scores."""
    mt = vds.variant_data.annotate_rows(prs=ht[vds.variant_data.locus])
    mt = mt.filter_rows(hl.is_defined(mt.prs))
    if mt.count_rows() == 0:
        print(f"[{prs_id}] WARNING: No overlapping variants; no scores computed")
        return None

    def dosage_fn(r):
        return hl.if_else(
            hl.is_missing(r.GT), 0,
            r.GT.n_alt_alleles()  # works if effect allele=ALT1; else more complex
        )

    mt = mt.annotate_entries(dosage=dosage_fn(mt))
    mt = mt.annotate_entries(contrib = mt.dosage * mt.prs.weight)
    ht_scores = mt.group_cols_by(mt.s).aggregate(
        total=hl.agg.sum(mt.contrib),
        count=hl.agg.count_where(mt.contrib != 0)
    )
    return ht_scores.annotate(
        normalized = hl.if_else(ht_scores.count>0, ht_scores.total/hl.float64(ht_scores.count), hl.missing(hl.tfloat64))
    )

def save_ht_and_csv(ht: hl.Table, hail_out: str, csv_out: str, prs_id: str, fs) -> None:
    """Write Hail Table and export CSV to GCS."""
    parent_ht = os.path.dirname(hail_out)
    parent_csv = os.path.dirname(csv_out)
    if not gcs_path_exists(parent_ht): fs.mkdirs(parent_ht, exist_ok=True)
    if not gcs_path_exists(parent_csv): fs.mkdirs(parent_csv, exist_ok=True)

    ht.write(hail_out, overwrite=True)
    print(f"[{prs_id}] Wrote HailTable to {hail_out}")

    ht_export = ht.select(
        person_id=ht.s,
        prs_total_score=ht.total,
        prs_variant_count=ht.count,
        prs_normalized_score=ht.normalized
    )
    ht_export.export(csv_out, header=True, delimiter=',', missing='NA')
    print(f"[{prs_id}] Exported scores CSV to {csv_out}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--prs_id',                      required=True)
    p.add_argument('--prs_url',                     required=True)
    p.add_argument('--base_cohort_vds_path',        required=True)
    p.add_argument('--gcs_temp_dir',                required=True)
    p.add_argument('--gcs_hail_temp_dir',           required=True)
    p.add_argument('--run_timestamp',               required=True)
    p.add_argument('--output_final_hail_table_gcs_path', required=True)
    p.add_argument('--output_final_score_csv_gcs_path',   required=True)
    args = p.parse_args()

    fs = get_gcs_fs()
    init_hail(args.gcs_hail_temp_dir, args.run_timestamp)

    # Intermediate GCS paths
    w = f"{args.gcs_temp_dir}/prs_weights/{args.prs_id}.csv"
    i = f"{args.gcs_temp_dir}/prs_intervals/{args.prs_id}.tsv"

    # Step 1: download & prepare
    needs_prep = not (gcs_path_exists(w) and gcs_path_exists(i))
    if needs_prep:
        gz = f"{args.prs_id}.gz"
        txt = f"{args.prs_id}.txt"
        local = download_file(args.prs_url, gz, txt, args.prs_id)
        df = prepare_weight_table(local, args.prs_id)
        os.remove(local)
        save_to_gcs(df, args.prs_id, "weights", w, fs)
        extract_intervals(df, args.prs_id, i, fs)

    # Step 2: load & filter VDS
    vds = hl.vds.read_vds(args.base_cohort_vds_path)
    vds_f = filter_vds_by_intervals(vds, i, args.prs_id)

    # Step 3: annotate & score
    ht_w = import_weights_as_ht(w, args.prs_id)
    ht_scores = annotate_and_score(vds_f, ht_w, args.prs_id)

    # Step 4: save outputs
    save_ht_and_csv(ht_scores,
                    args.output_final_hail_table_gcs_path,
                    args.output_final_score_csv_gcs_path,
                    args.prs_id,
                    fs)

if __name__ == "__main__":
    main()
