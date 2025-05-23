import argparse
import datetime
import os
import sys
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.linear_model import LogisticRegression
import scipy.stats as stats
from utils import gcs_path_exists, get_gcs_fs

def main():
    # Define expected arguments for analyzing a single model's results
    parser = argparse.ArgumentParser(description="Analyze and visualize PRS results for a single model.")
    parser.add_argument("--prs_id", required=True, help="PGS ID of the model (e.g., PGS000123).")
    parser.add_argument("--prs_phenotype_label", required=True, help="Phenotype label for this PGS model (e.g., 'Ischemic Stroke') used for plots.")
    parser.add_argument("--score_csv_gcs_path", required=True, help="GCS path to the calculated PRS scores CSV for this specific model.")
    parser.add_argument("--wgs_ehr_ids_gcs_path", required=True, help="GCS path to the CSV containing all WGS+EHR sample IDs ('s' column) used for the VDS.")
    parser.add_argument("--phenotype_cases_csv_gcs_path", required=True, help="GCS path to the CSV containing all phenotype case sample IDs ('s', 'phenotype_status' columns) for the target condition.")
    parser.add_argument("--phenotype_name", required=True, help="Name of the target phenotype condition being analyzed against (e.g., 'Ischemic Stroke').")
    parser.add_argument("--gcs_base_output_dir_run", required=True, help="GCS base output directory for this specific pipeline run (used for saving plots, standardized scores, etc.).")
    parser.add_argument("--run_timestamp", required=True, help="Pipeline run timestamp (YYYYMMDD_HHMMSS), used for unique local directories if needed.")
    parser.add_argument("--google_billing_project", required=True, help="Google Cloud Project ID for billing and GCS access.")
    parser.add_argument("--output_summary_file_name", required=True, help="Local file name to write the analysis summary statistics.")

    args = parser.parse_args()

    # Assign arguments to variables for clarity
    prs_id = args.prs_id
    model_phenotype_label = args.prs_phenotype_label # Label specific to the model (from models.csv)
    target_phenotype_name = args.phenotype_name # The condition being predicted
    phenotype_col_name = f"{target_phenotype_name.replace(' ', '_')}_status" # Standardized column name for phenotype status

    # Initialize GCS filesystem, explicitly providing the billing project
    fs = get_gcs_fs(project_id_for_billing=args.google_billing_project) # Initialize GCS filesystem from utils

    print(f"\n------- Analyzing Results for PRS Model: {prs_id} ({target_phenotype_name}) -------")
    analysis_start_time = datetime.datetime.now()

    # --- Load Full WGS/EHR Cohort List ---
    print(f"[{prs_id}] Loading WGS+EHR cohort list from: {args.wgs_ehr_ids_gcs_path}")
    if not gcs_path_exists(args.wgs_ehr_ids_gcs_path): # common_utils
        print(f"[{prs_id}] FATAL ERROR: WGS/EHR sample list not found at {args.wgs_ehr_ids_gcs_path}.")
        sys.exit(1)
    try:
        with fs.open(args.wgs_ehr_ids_gcs_path, 'r') as f:
            # Assuming the CSV from prepare_base_vds.py has 'person_id' column
            wgs_ehr_full_df = pd.read_csv(f).rename(columns={'person_id': 's'})
        if wgs_ehr_full_df.empty:
            print(f"[{prs_id}] FATAL ERROR: WGS/EHR sample list at {args.wgs_ehr_ids_gcs_path} is empty.")
            sys.exit(1)
        wgs_ehr_full_df['s'] = wgs_ehr_full_df['s'].astype(str)
        print(f"[{prs_id}] Loaded {len(wgs_ehr_full_df)} total individuals from WGS+EHR cohort list.")
    except Exception as e:
        print(f"[{prs_id}] FATAL ERROR loading WGS/EHR cohort list: {e}"); sys.exit(1)

    # --- Load Phenotype Cases Data ---
    print(f"[{prs_id}] Loading phenotype case data from: {args.phenotype_cases_csv_gcs_path}")
    if not gcs_path_exists(args.phenotype_cases_csv_gcs_path): # common_utils
        print(f"[{prs_id}] FATAL ERROR: Phenotype cases CSV not found at {args.phenotype_cases_csv_gcs_path}.")
        sys.exit(1)
    try:
        with fs.open(args.phenotype_cases_csv_gcs_path, 'r') as f:
            phenotype_cases_df = pd.read_csv(f) # Expects 's' and 'phenotype_status'
        if phenotype_cases_df.empty and 's' not in phenotype_cases_df.columns: # if empty, it might not have columns
             phenotype_cases_df = pd.DataFrame(columns=['s', 'phenotype_status']) # ensure schema for merge
        elif 's' not in phenotype_cases_df.columns:
            print(f"[{prs_id}] FATAL ERROR: Phenotype cases CSV {args.phenotype_cases_csv_gcs_path} must contain an 's' column.")
            sys.exit(1)

        phenotype_cases_df['s'] = phenotype_cases_df['s'].astype(str)
        # Ensure 'phenotype_status' exists, even if no cases were found (it would be empty then)
        if 'phenotype_status' not in phenotype_cases_df.columns:
            phenotype_cases_df['phenotype_status'] = 1 # Assume all in this file are cases if column is missing. Or handle as error.
                                                      # The fetch_phenotypes script should always produce this column.
        print(f"[{prs_id}] Loaded {len(phenotype_cases_df)} entries from phenotype case data.")
    except Exception as e:
        print(f"[{prs_id}] FATAL ERROR loading phenotype case data: {e}"); sys.exit(1)

    # --- Create Full Phenotype DataFrame (Cases + Controls from WGS/EHR cohort) ---
    phenotype_df_full = pd.merge(wgs_ehr_full_df[['s']], phenotype_cases_df[['s', 'phenotype_status']], on='s', how='left')
    phenotype_df_full['phenotype_status'] = phenotype_df_full['phenotype_status'].fillna(0).astype(int)
    phenotype_df_full = phenotype_df_full.rename(columns={'phenotype_status': phenotype_col_name})
    print(f"\n[{prs_id}] Phenotype distribution in analyzed cohort ({target_phenotype_name}):")
    print(phenotype_df_full[phenotype_col_name].value_counts())
    print(f"Total N in cohort for analysis: {len(phenotype_df_full)}")


    # --- Load PRS Scores for the current model ---
    print(f"[{prs_id}] Loading final PRS scores from: {args.score_csv_gcs_path}")
    if not gcs_path_exists(args.score_csv_gcs_path): # common_utils
        print(f"[{prs_id}] FATAL ERROR: Score CSV for model {prs_id} not found at {args.score_csv_gcs_path}.")
        sys.exit(1)
    try:
        with fs.open(args.score_csv_gcs_path, 'r') as f:
            prs_scores_df = pd.read_csv(f)
        print(f"[{prs_id}] Loaded PRS scores for {prs_id}. Shape: {prs_scores_df.shape}")
        if prs_scores_df.empty:
            print(f"[{prs_id}] WARNING: Loaded PRS score file for {prs_id} is empty. Skipping analysis for this model.")
            sys.exit(0) # Not a fatal error for the pipeline, just this branch.

        # Standardize sample ID column name to 's'
        if 's' not in prs_scores_df.columns:
            if 'person_id' in prs_scores_df.columns:
                prs_scores_df = prs_scores_df.rename(columns={'person_id': 's'})
            else:
                raise ValueError("Missing required sample ID column ('s' or 'person_id') in score file.")
        prs_scores_df['s'] = prs_scores_df['s'].astype(str)
    except Exception as e:
        print(f"[{prs_id}] ERROR loading score file {args.score_csv_gcs_path}: {e}"); sys.exit(1)


    # --- Merge Scores with Phenotype ---
    print(f"[{prs_id}] Merging PRS scores with {target_phenotype_name} phenotype data...")
    try:
        # Inner join: only keep samples present in both scores and the full phenotype definition
        merged_scores_phenotypes = pd.merge(prs_scores_df, phenotype_df_full, on='s', how='inner')
        print(f"[{prs_id}] Merged DataFrame shape: {merged_scores_phenotypes.shape}")
        if merged_scores_phenotypes.empty:
            print(f"[{prs_id}] WARNING: Merged dataframe is empty (no overlap between scores and phenotype cohort). Skipping analysis for {prs_id}.")
            sys.exit(0) # Not fatal for pipeline
    except Exception as e:
        print(f"[{prs_id}] ERROR during merge of scores and phenotype: {e}"); sys.exit(1)


    # --- PRS Score Processing & Standardization ---
    print(f"[{prs_id}] Processing PRS scores (Z-score, percentiles, categories)...")
    score_column_for_analysis = 'prs_total_score' # Or 'prs_normalized_score'
    if score_column_for_analysis not in merged_scores_phenotypes.columns:
        print(f"[{prs_id}] ERROR: Score column '{score_column_for_analysis}' not found. Available: {merged_scores_phenotypes.columns.tolist()}")
        sys.exit(1)

    merged_scores_phenotypes = merged_scores_phenotypes.replace([np.inf, -np.inf], np.nan)
    merged_scores_phenotypes[phenotype_col_name] = pd.to_numeric(merged_scores_phenotypes[phenotype_col_name], errors='coerce')
    rows_before_drop = len(merged_scores_phenotypes)
    merged_scores_phenotypes = merged_scores_phenotypes.dropna(subset=[score_column_for_analysis, phenotype_col_name])
    rows_after_drop = len(merged_scores_phenotypes)
    if rows_after_drop < rows_before_drop:
        print(f"[{prs_id}] WARNING: Dropped {rows_before_drop - rows_after_drop} rows with NaN/Inf scores or invalid phenotype.")
    if merged_scores_phenotypes.empty:
        print(f"[{prs_id}] ERROR: No valid data remaining after cleaning scores/phenotypes."); sys.exit(1)
    merged_scores_phenotypes[phenotype_col_name] = merged_scores_phenotypes[phenotype_col_name].astype(int)

    mean_score = merged_scores_phenotypes[score_column_for_analysis].mean()
    std_score = merged_scores_phenotypes[score_column_for_analysis].std()
    can_do_quantiles = True
    if std_score == 0 or pd.isna(std_score) or len(merged_scores_phenotypes[score_column_for_analysis].unique()) < 10: # Added check for unique values
        print(f"[{prs_id}] WARNING: Std dev of PRS is zero/NaN ({std_score}) or too few unique scores. Assigning default standardized values/skipping quantiles.")
        merged_scores_phenotypes['prs_zscore'] = 0.0
        merged_scores_phenotypes['prs_percentile'] = 50.0
        merged_scores_phenotypes['risk_category'] = 'Average'
        merged_scores_phenotypes['prs_decile'] = 'D5' # Placeholder
        can_do_quantiles = False
    else:
        merged_scores_phenotypes['prs_zscore'] = (merged_scores_phenotypes[score_column_for_analysis] - mean_score) / std_score
        merged_scores_phenotypes['prs_percentile'] = merged_scores_phenotypes['prs_zscore'].rank(pct=True) * 100
        try:
            merged_scores_phenotypes['risk_category'] = pd.qcut(merged_scores_phenotypes['prs_zscore'], q=5, labels=[f'Q{i}' for i in range(1,6)], duplicates='drop')
            merged_scores_phenotypes['prs_decile'] = pd.qcut(merged_scores_phenotypes['prs_zscore'], q=10, labels=[f'D{i}' for i in range(1,11)], duplicates='drop')
        except ValueError as e:
            print(f"[{prs_id}] WARNING: Could not create quantiles (likely too few unique scores/samples). Error: {e}.")
            merged_scores_phenotypes['risk_category'] = 'NA'
            merged_scores_phenotypes['prs_decile'] = 'NA'
            can_do_quantiles = False

    print(f"[{prs_id}] Descriptive stats for '{score_column_for_analysis}':\n{merged_scores_phenotypes[score_column_for_analysis].describe()}\n")
    if 'prs_zscore' in merged_scores_phenotypes:
        print(f"[{prs_id}] Descriptive stats for PRS Z-score:\n{merged_scores_phenotypes['prs_zscore'].describe()}\n")

    # --- Save Standardized Scores ---
    # Output path relative to this run's GCS output dir
    standardized_scores_gcs_path = f"{args.gcs_base_output_dir_run}/scores/{prs_id}/calculated_scores/score/{prs_id}_scores_standardized.csv"
    print(f"[{prs_id}] Saving standardized PRS scores for this run to: {standardized_scores_gcs_path}")
    try:
        std_scores_dir_gcs = os.path.dirname(standardized_scores_gcs_path)
        if not gcs_path_exists(std_scores_dir_gcs): # common_utils
            print(f"[{prs_id}] Creating GCS directory for standardized scores: {std_scores_dir_gcs}")
            fs.mkdirs(std_scores_dir_gcs, exist_ok=True)
        with fs.open(standardized_scores_gcs_path, 'w') as f:
            merged_scores_phenotypes.to_csv(f, index=False, na_rep='NA')
        print(f"[{prs_id}] Standardized scores saved.\n")
    except Exception as e:
        print(f"[{prs_id}] ERROR saving standardized scores: {e}\n") # Non-fatal for subsequent analysis with in-memory df


    # --- Visualization & Performance Evaluation ---
    sns.set_style("whitegrid"); sns.set_context("talk")
    palette = sns.color_palette("viridis", as_cmap=True)
    # Local plot directory within Nextflow's task work directory
    plots_dir_local = f"plots_local_{prs_id}_{args.run_timestamp}"
    os.makedirs(plots_dir_local, exist_ok=True)
    # GCS plot directory within this run's output
    plots_dir_gcs = f"{args.gcs_base_output_dir_run}/scores/{prs_id}/plots"

    # Plot PRS Distribution
    if 'prs_zscore' in merged_scores_phenotypes:
        print(f"[{prs_id}] Visualizing PRS distributions...")
        try:
            fig_dist, axs_dist = plt.subplots(1, 2, figsize=(18, 7)); plt.subplots_adjust(wspace=0.3)
            sns.histplot(data=merged_scores_phenotypes, x='prs_zscore', bins=50, ax=axs_dist[0], color=palette(0.2))
            axs_dist[0].set_title(f'{prs_id} PRS Z-score Distribution', fontsize=16, fontweight='bold')
            axs_dist[0].set_xlabel('PRS Z-score', fontsize=14); axs_dist[0].set_ylabel('Frequency', fontsize=14)
            sns.kdeplot(data=merged_scores_phenotypes, x='prs_zscore', fill=True, ax=axs_dist[1], color=palette(0.6))
            axs_dist[1].set_title(f'{prs_id} PRS Z-score Density', fontsize=16, fontweight='bold')
            axs_dist[1].set_xlabel('PRS Z-score', fontsize=14); axs_dist[1].set_ylabel('Density', fontsize=14)
            plt.suptitle(f"{target_phenotype_name} PRS ({prs_id}) Distribution", fontsize=20, fontweight='bold')
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            dist_plot_local_path = os.path.join(plots_dir_local, f"{prs_id}_distribution_plots.png")
            plt.savefig(dist_plot_local_path)
            print(f"[{prs_id}] Saved distribution plot locally: {dist_plot_local_path}")
            dist_plot_gcs_path = f"{plots_dir_gcs}/{prs_id}_distribution_plots.png"
            fs.put(dist_plot_local_path, dist_plot_gcs_path) # Upload
            print(f"[{prs_id}] Uploaded distribution plot to: {dist_plot_gcs_path}")
            plt.close(fig_dist)
        except Exception as e: print(f"[{prs_id}] ERROR generating distribution plot: {e}"); plt.close('all')
    else: print(f"[{prs_id}] Skipping distribution plots as 'prs_zscore' is not available.")

    # ROC-AUC
    auc_val = np.nan
    if 'prs_zscore' in merged_scores_phenotypes and merged_scores_phenotypes[phenotype_col_name].nunique() >= 2:
        print(f"[{prs_id}] Calculating ROC-AUC for {target_phenotype_name}...")
        try:
            y_true = merged_scores_phenotypes[phenotype_col_name]
            y_score = merged_scores_phenotypes['prs_zscore']
            auc_val = roc_auc_score(y_true, y_score)
            fpr_arr, tpr_arr, _ = roc_curve(y_true, y_score)
            plt.figure(figsize=(10, 8))
            sns.lineplot(x=fpr_arr, y=tpr_arr, color='blue', lw=2, label=f'ROC (AUC = {auc_val:.3f})')
            sns.lineplot(x=[0,1], y=[0,1], color='gray', linestyle='--', label='Random (AUC = 0.500)')
            plt.xlim([0, 1]); plt.ylim([0, 1.05])
            plt.xlabel('False Positive Rate (1 - Specificity)', fontsize=14); plt.ylabel('True Positive Rate (Sensitivity)', fontsize=14)
            plt.title(f'{target_phenotype_name} PRS ({prs_id}) ROC Curve', fontsize=18, fontweight='bold')
            plt.legend(loc="lower right"); plt.tight_layout()
            roc_plot_local_path = os.path.join(plots_dir_local, f"{prs_id}_roc_curve.png")
            plt.savefig(roc_plot_local_path)
            print(f"[{prs_id}] Saved ROC plot locally: {roc_plot_local_path}")
            roc_plot_gcs_path = f"{plots_dir_gcs}/{prs_id}_roc_curve.png"
            fs.put(roc_plot_local_path, roc_plot_gcs_path) # Upload
            print(f"[{prs_id}] Uploaded ROC plot to: {roc_plot_gcs_path}")
            plt.close()
            print(f"[{prs_id}] ROC-AUC = {auc_val:.4f}")
        except ValueError as e: print(f"[{prs_id}] ERROR ROC (ValueError): {e}. Check data diversity.")
        except Exception as e: print(f"[{prs_id}] ERROR generating ROC plot: {e}"); plt.close('all')
    else: print(f"[{prs_id}] Skipping ROC: Z-score not available or only one phenotype class.")

    # Prevalence by Quantile
    if can_do_quantiles and 'prs_decile' in merged_scores_phenotypes and merged_scores_phenotypes['prs_decile'].nunique() > 1 : # Check if deciles are varied
        print(f"\n[{prs_id}] Analyzing {target_phenotype_name} prevalence across PRS deciles...")
        try:
            decile_stats = merged_scores_phenotypes.groupby('prs_decile', observed=False).agg(
                count=(phenotype_col_name, 'size'),
                cases=(phenotype_col_name, 'sum')
            ).reset_index()
            decile_stats['prevalence_percent'] = np.where(
                decile_stats['count'] > 0, (decile_stats['cases'] / decile_stats['count']) * 100, 0
            )
            print(f"[{prs_id}] {target_phenotype_name} Prevalence by PRS Decile:\n{decile_stats[['prs_decile', 'count', 'cases', 'prevalence_percent']].round(2)}\n")
            plt.figure(figsize=(12, 7))
            bar_plot = sns.barplot(data=decile_stats, x='prs_decile', y='prevalence_percent', palette='viridis')
            plt.title(f'{target_phenotype_name} Prevalence by PRS Decile ({prs_id})', fontsize=18, fontweight='bold')
            plt.xlabel('PRS Decile (D1=Lowest Risk, D10=Highest Risk)', fontsize=14); plt.ylabel(f'{target_phenotype_name} Prevalence (%)', fontsize=14)
            for index, row in decile_stats.iterrows():
                if pd.notna(row['prevalence_percent']):
                    bar_plot.text(index, row['prevalence_percent'] + 0.2, f"{row['prevalence_percent']:.2f}%", color='black', ha="center", va='bottom', fontsize=10)
            plt.tight_layout()
            prev_plot_local_path = os.path.join(plots_dir_local, f"{prs_id}_prevalence_by_decile.png")
            plt.savefig(prev_plot_local_path)
            print(f"[{prs_id}] Saved prevalence plot locally: {prev_plot_local_path}")
            prev_plot_gcs_path = f"{plots_dir_gcs}/{prs_id}_prevalence_by_decile.png"
            fs.put(prev_plot_local_path, prev_plot_gcs_path) # Upload
            print(f"[{prs_id}] Uploaded prevalence plot to: {prev_plot_gcs_path}")
            plt.close()
        except Exception as e: print(f"[{prs_id}] Error plotting prevalence by decile: {e}"); plt.close('all')
    else: print(f"[{prs_id}] Skipping prevalence by decile (quantiles not calculated or not varied).")

    # Odds Ratio per SD
    odds_ratio_per_SD = np.nan
    if 'prs_zscore' in merged_scores_phenotypes and std_score != 0 and not pd.isna(std_score) and merged_scores_phenotypes[phenotype_col_name].nunique() >= 2:
        print(f"\n[{prs_id}] Calculating odds ratio per standard deviation (SD) increase in PRS...")
        try:
            X = merged_scores_phenotypes[['prs_zscore']]
            y = merged_scores_phenotypes[phenotype_col_name]
            log_model = LogisticRegression(solver='liblinear', C=1.0, class_weight='balanced', penalty='l2', random_state=42)
            log_model.fit(X, y)
            log_odds_per_SD = log_model.coef_[0][0]
            odds_ratio_per_SD = np.exp(log_odds_per_SD)
            print(f"[{prs_id}] Odds Ratio per 1 SD increase in PRS Z-score: {odds_ratio_per_SD:.3f}\n")
        except Exception as e: print(f"[{prs_id}] ERROR calculating OR per SD: {e}")
    else: print(f"[{prs_id}] Skipping OR per SD (Z-score invalid, SD zero, or one phenotype class).")

    # Top vs Bottom 5%
    odds_ratio_5_v_95 = np.nan
    if 'prs_zscore' in merged_scores_phenotypes and std_score != 0 and not pd.isna(std_score) and merged_scores_phenotypes[phenotype_col_name].nunique() >= 2:
        print(f"\n[{prs_id}] Comparing Top 5% vs Bottom 5% of PRS distribution...")
        try:
            zscores_clean = merged_scores_phenotypes['prs_zscore'].dropna()
            if len(zscores_clean) < 40: # Need enough for 5th and 95th percentiles to be meaningful & distinct
                print(f"[{prs_id}] WARNING: Too few samples ({len(zscores_clean)}) for reliable 5th/95th percentile comparison.")
            else:
                top_5_threshold = np.percentile(zscores_clean, 95)
                bottom_5_threshold = np.percentile(zscores_clean, 5)
                if top_5_threshold <= bottom_5_threshold: # Check if thresholds are distinct
                    print(f"[{prs_id}] WARNING: Top 5% ({top_5_threshold:.3f}) and Bottom 5% ({bottom_5_threshold:.3f}) thresholds are not distinct. Cannot compare.")
                else:
                    lowest_5 = merged_scores_phenotypes[merged_scores_phenotypes['prs_zscore'] <= bottom_5_threshold]
                    highest_5 = merged_scores_phenotypes[merged_scores_phenotypes['prs_zscore'] >= top_5_threshold]
                    if lowest_5.empty or highest_5.empty:
                        print(f"[{prs_id}] WARNING: Top or Bottom 5% group is empty. Cannot compare.")
                    else:
                        bottom_5_count = len(lowest_5); bottom_5_cases = lowest_5[phenotype_col_name].sum()
                        top_5_count = len(highest_5); top_5_cases = highest_5[phenotype_col_name].sum()
                        bottom_5_prev = (bottom_5_cases / bottom_5_count * 100) if bottom_5_count > 0 else 0
                        top_5_prev = (top_5_cases / top_5_count * 100) if top_5_count > 0 else 0
                        print(f"Bottom 5% (Z <= {bottom_5_threshold:.3f}): Prev={bottom_5_prev:.2f}%, N={bottom_5_count}, Cases={bottom_5_cases}")
                        print(f"Top 5% (Z >= {top_5_threshold:.3f}): Prev={top_5_prev:.2f}%, N={top_5_count}, Cases={top_5_cases}")
                        
                        # Contingency table: [[HighRisk_Cases, HighRisk_Controls], [LowRisk_Cases, LowRisk_Controls]]
                        table = np.array([
                            [top_5_cases, top_5_count - top_5_cases],
                            [bottom_5_cases, bottom_5_count - bottom_5_cases]
                        ])
                        if np.any(table < 0): print(f"[{prs_id}] ERROR: Invalid counts in contingency table: {table}.")
                        else:
                            try:
                                if np.any(table == 0): # Haldane-Anscombe correction for zeros
                                    print(f"[{prs_id}] Applying 0.5 correction for Fisher's exact test due to zero cell.")
                                    odds_ratio_5_v_95, p_value = stats.fisher_exact(table + 0.5)
                                else:
                                    odds_ratio_5_v_95, p_value = stats.fisher_exact(table)
                                print(f"Odds Ratio (Top 5% vs Bottom 5%): {odds_ratio_5_v_95:.2f} (p={p_value:.2g})\n")
                            except ValueError as fe_err: print(f"[{prs_id}] ERROR Fisher's exact: {fe_err}"); odds_ratio_5_v_95 = np.nan

                        plt.figure(figsize=(8, 6))
                        group_labels = ['Bottom 5%', 'Top 5%']; group_prevalences = [bottom_5_prev, top_5_prev]
                        bars = plt.bar(group_labels, group_prevalences, color=['#1f77b4', '#d62728'])
                        if bottom_5_count > 0: plt.text(bars[0].get_x() + bars[0].get_width()/2., bars[0].get_height(), f'N={bottom_5_count}\nCases={bottom_5_cases}', ha='center', va='bottom', fontsize=12)
                        if top_5_count > 0: plt.text(bars[1].get_x() + bars[1].get_width()/2., bars[1].get_height(), f'N={top_5_count}\nCases={top_5_cases}', ha='center', va='bottom', fontsize=12)
                        plt.title(f'{target_phenotype_name} Prev in Top vs Bottom 5% PRS ({prs_id})', fontsize=16, fontweight='bold', pad=20)
                        plt.ylabel(f'{target_phenotype_name} Prevalence (%)', fontsize=14); plt.ylim(0, max(group_prevalences or [1]) * 1.2 + 5)
                        plt.tight_layout()
                        comp_plot_local_path = os.path.join(plots_dir_local, f"{prs_id}_top_bottom_comparison.png")
                        plt.savefig(comp_plot_local_path)
                        print(f"[{prs_id}] Saved top/bottom comparison plot locally: {comp_plot_local_path}")
                        comp_plot_gcs_path = f"{plots_dir_gcs}/{prs_id}_top_bottom_comparison.png"
                        fs.put(comp_plot_local_path, comp_plot_gcs_path) # Upload
                        print(f"[{prs_id}] Uploaded top/bottom comparison plot to: {comp_plot_gcs_path}")
                        plt.close()
        except Exception as e: print(f"[{prs_id}] ERROR Top vs Bottom 5% comparison: {e}"); plt.close('all')
    else: print(f"[{prs_id}] Skipping Top vs Bottom 5% (Z-score invalid, SD zero, or one phenotype class).")

    # --- Summary ---
    print(f"\n--- {prs_id} ANALYSIS SUMMARY ({target_phenotype_name}) ---")
    print(f"Score Column Used: {score_column_for_analysis}")
    print(f"N Analyzed: {len(merged_scores_phenotypes)}")
    print(f"ROC-AUC: {auc_val:.4f}" if not pd.isna(auc_val) else "ROC-AUC: Not Calculated")
    print(f"Odds Ratio per 1 SD(Z-score): {odds_ratio_per_SD:.3f}" if not pd.isna(odds_ratio_per_SD) else "OR per SD: Not Calculated")
    print(f"Odds Ratio (Top 5% vs Bottom 5%): {odds_ratio_5_v_95:.2f}" if not pd.isna(odds_ratio_5_v_95) else "OR Top/Bottom 5%: Not Calculated")
    analysis_end_time = datetime.datetime.now()
    print(f"Analysis duration for {prs_id}: {analysis_end_time - analysis_start_time}")
    print(f"------- Finished analysis for {prs_id} -------")

    # Clean up local plot directory
    if os.path.exists(plots_dir_local):
        shutil.rmtree(plots_dir_local, ignore_errors=True)
        print(f"[{prs_id}] Cleaned up local plot directory: {plots_dir_local}")

    # Create an output file that Nextflow can track for this process completion for this model
    # The actual outputs (plots, standardized CSV) are written to GCS directly.
    # This output file name should be unique per prs_id.
    completion_marker_filename = f"{prs_id}_analysis_complete.txt"
    with open(completion_marker_filename, "w") as f:
        f.write(f"Analysis for {prs_id} completed at {datetime.datetime.now()}.\n")
        f.write(f"ROC-AUC: {auc_val:.4f}\n")
        f.write(f"OR per SD: {odds_ratio_per_SD:.3f}\n")
        f.write(f"OR Top/Bottom 5%: {odds_ratio_5_v_95:.2f}\n")
    print(f"[{prs_id}] Analysis completion marker created: {completion_marker_filename}")

if __name__ == "__main__":
    main()
