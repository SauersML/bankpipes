"""
main.py  –  stand-alone orchestrator for the AoU PRS pipeline.

*   Executes Steps:
        1) fetch_phenotypes.py
        2) prepare_base_vds.py
        3) For each PGS entry:
               process_prs_model.py   →   analyze_prs_results.py
        4) Collate summaries
*   Each model runs in its own executor worker; workers == PARALLEL_WORKERS
    env-var (else = #models).  No artificial CPU/RAM caps.
"""

import argparse
import concurrent.futures
import datetime as _dt
import json
import logging
import os
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd
import psutil
import yaml

# ────────────────────────────────────────────────────────────────────────────────
#  CONSTANTS
# ────────────────────────────────────────────────────────────────────────────────
MODELS_CSV_FILENAME = "models.csv"
PYTHON_EXECUTABLE = "/opt/conda/bin/python3"

FLAGGED_SAMPLES_FILE_SUFFIX = (
    "/wgs/short_read/snpindel/aux/relatedness/relatedness_flagged_samples.tsv"
)

GCS_REUSABLE_INTERMEDIATES_SUFFIX = "pgs_pipeline_python_orchestrator/intermediates"
GCS_HAIL_TEMP_RUN_SPECIFIC_SUFFIX = "pgs_pipeline_python_orchestrator/hail_temp"
GCS_RUN_OUTPUTS_SUFFIX = "pgs_pipeline_python_orchestrator/runs"

VDS_PREP_ENABLE_DOWNSAMPLING = True
VDS_PREP_N_CASES_DOWNSAMPLE = 500
VDS_PREP_N_CONTROLS_DOWNSAMPLE = 500
VDS_PREP_DOWNSAMPLING_RANDOM_STATE = 2025

LOG_FMT = "%(asctime)s [%(levelname)s] %(module)s – %(message)s"
logging.basicConfig(level=logging.INFO, format=LOG_FMT, handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger("main")


# ────────────────────────────────────────────────────────────────────────────────
#  UTILITIES
# ────────────────────────────────────────────────────────────────────────────────
def _run_subprocess(cmd: list[str], env: dict[str, str], log_path: Path, step: str) -> int:
    """Run *cmd* synchronously, write its output to *log_path*, and stream both status and new log lines to the console."""
    with log_path.open("w") as fh:
        proc = subprocess.Popen(cmd, env=env, stdout=fh, stderr=subprocess.STDOUT, text=True)

    log.info("%s  –  spawned PID %d  → %s", step, proc.pid, log_path)

    # open a second handle for tail-follow
    reader = log_path.open("r")
    reader.seek(0, os.SEEK_END)

    start_time = time.time()
    last_report = start_time

    while proc.poll() is None:
        time.sleep(5)

        # stream any new log content
        new_data = reader.read()
        if new_data:
            sys.stdout.write(new_data)
            sys.stdout.flush()

        now = time.time()
        if now - last_report >= 15:
            cpu = psutil.cpu_percent()
            ram = psutil.virtual_memory().percent
            disk = psutil.disk_usage("/").percent
            elapsed = now - start_time
            log.info(
                "%s  –  running %.0fs  CPU %.1f %%  RAM %.1f %%  Disk %.1f %%",
                step, elapsed, cpu, ram, disk
            )
            last_report = now

    # flush any remaining log output
    final_data = reader.read()
    if final_data:
        sys.stdout.write(final_data)
        sys.stdout.flush()
    reader.close()

    rc = proc.wait()
    if rc:
        try:
            err_text = log_path.read_text()
            sys.stderr.write(f"\n--- {step} STDERR ---\n{err_text}\n--- END STDERR ---\n")
        except Exception:
            log.error("Could not read %s for error output.", log_path)
    return rc

def _build_env(base_py: str, script_dir: Path) -> dict[str, str]:
    """Return env dict for child scripts."""
    env = os.environ.copy()
    env["PYSPARK_PYTHON"] = base_py
    env["PYTHONPATH"] = f"{script_dir/'src'}{os.pathsep}{env.get('PYTHONPATH','')}"
    return env


# ────────────────────────────────────────────────────────────────────────────────
#  CONFIG CLASS
# ────────────────────────────────────────────────────────────────────────────────
class Config:
    def __init__(self, yaml_path: Path):
        self.cfg = yaml.safe_load(Path(yaml_path).read_text())
        pheno = self.cfg.get("phenotype_definition") or {}
        self.pheno_name: str = pheno.get("target_name") or sys.exit("target_name missing")
        self.pheno_concept_ids: list[int] = pheno.get("concept_ids") or sys.exit("concept_ids missing")
        self.env = {k: os.getenv(k) for k in (
            "GOOGLE_PROJECT", "WORKSPACE_BUCKET", "WGS_VDS_PATH",
            "WORKSPACE_CDR", "CDR_STORAGE_PATH"
        )}
        if any(v is None for v in self.env.values()):
            missing = [k for k, v in self.env.items() if v is None]
            sys.exit(f"Missing required env vars: {', '.join(missing)}")
        self.bucket = self.env["WORKSPACE_BUCKET"].removeprefix("gs://")
        ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
        self.ts = ts
        self.script_dir = Path(__file__).resolve().parent
        # GCS paths
        self.gcs_intermediate = f"gs://{self.bucket}/{GCS_REUSABLE_INTERMEDIATES_SUFFIX}"
        self.gcs_hail_tmp = f"gs://{self.bucket}/{GCS_HAIL_TEMP_RUN_SPECIFIC_SUFFIX}/{ts}"
        self.gcs_outputs = f"gs://{self.bucket}/{GCS_RUN_OUTPUTS_SUFFIX}/{ts}"
        # Spark conf
        spark_conf = [
            "spark.hadoop.fs.gs.requester.pays.mode=AUTO",
            f"spark.hadoop.fs.gs.requester.pays.project.id={self.env['GOOGLE_PROJECT']}",
            f"spark.hadoop.fs.gs.project.id={self.env['GOOGLE_PROJECT']}",
            "spark.hadoop.fs.gs.block.size=134217728",
            # leave dynamicAllocation unset (not supported) → cluster alloc determines resources
            "spark.serializer=org.apache.spark.serializer.KryoSerializer",
            "spark.executor.defaultJavaOptions=-XX:+UseG1GC -XX:InitiatingHeapOccupancyPercent=35"
        ]
        extra = os.getenv("EXTRA_SPARK_CONF")
        if extra:
            try:
                spark_conf.extend(json.loads(extra))
            except Exception as e:
                log.warning("Could not parse EXTRA_SPARK_CONF (%s)", e)
        self.spark_conf_json = json.dumps(spark_conf)
        # Models CSV
        self.models_csv = self.script_dir / MODELS_CSV_FILENAME
        if not self.models_csv.exists():
            sys.exit(f"Missing {self.models_csv}")
        log.info("Config initialised (%s models listed)", sum(1 for _ in self.models_csv.open()) - 1)


# ────────────────────────────────────────────────────────────────────────────────
#  MAIN PIPELINE
# ────────────────────────────────────────────────────────────────────────────────
def main(cfg: Config) -> None:
    work_dir = Path.cwd() / f"pgs_pipeline_work_{cfg.ts}"
    log_dir = work_dir / "logs"
    work_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)
    env = _build_env(PYTHON_EXECUTABLE, cfg.script_dir)

    # ── Step 1 – fetch_phenotypes ──────────────────────────────────────────────
    pheno_csv_gcs = f"{cfg.gcs_intermediate}/phenotype/{cfg.pheno_name.replace(' ','_')}_cases.csv"
    rc = _run_subprocess(
        [
            PYTHON_EXECUTABLE, cfg.script_dir/"src"/"fetch_phenotypes.py",
            "--phenotype_name", cfg.pheno_name,
            "--phenotype_concept_ids", ",".join(map(str, cfg.pheno_concept_ids)),
            "--workspace_cdr", cfg.env["WORKSPACE_CDR"],
            "--output_phenotype_csv_gcs_path", pheno_csv_gcs,
            "--google_billing_project", cfg.env["GOOGLE_PROJECT"],
        ],
        env, log_dir/"01_fetch_phenotypes.log", "fetch_phenotypes"
    )
    if rc:
        sys.exit("fetch_phenotypes failed")

    # ── Step 2 – prepare_base_vds ──────────────────────────────────────────────
    base_vds_gcs = f"{cfg.gcs_intermediate}/base_vds/base_cohort.vds"
    wgs_ehr_ids_gcs = f"{cfg.gcs_intermediate}/base_vds/cohort_ids.csv"
    prep_args = [
        PYTHON_EXECUTABLE, cfg.script_dir/"src"/"prepare_base_vds.py",
        "--workspace_cdr", cfg.env["WORKSPACE_CDR"],
        "--run_timestamp", cfg.ts,
        "--gcs_temp_dir", cfg.gcs_intermediate,
        "--gcs_hail_temp_dir", cfg.gcs_hail_tmp,
        "--wgs_vds_path", cfg.env["WGS_VDS_PATH"],
        "--flagged_samples_gcs_path", f"{cfg.env['CDR_STORAGE_PATH'].rstrip('/')}{FLAGGED_SAMPLES_FILE_SUFFIX}",
        "--base_cohort_vds_path_out", base_vds_gcs,
        "--wgs_ehr_ids_gcs_path_out", wgs_ehr_ids_gcs,
        "--target_phenotype_name", cfg.pheno_name,
        "--phenotype_cases_gcs_path_input", pheno_csv_gcs,
        "--n_cases_downsample", str(VDS_PREP_N_CASES_DOWNSAMPLE),
        "--n_controls_downsample", str(VDS_PREP_N_CONTROLS_DOWNSAMPLE),
        "--downsampling_random_state", str(VDS_PREP_DOWNSAMPLING_RANDOM_STATE),
        "--google_billing_project", cfg.env["GOOGLE_PROJECT"],
        "--spark_configurations_json", cfg.spark_conf_json,
        # Instruct Hail to run on the Dataproc YARN cluster
        "--hail_cluster_mode", "dataproc_yarn"
    ]
    if VDS_PREP_ENABLE_DOWNSAMPLING:

    rc = _run_subprocess(prep_args, env, log_dir/"02_prepare_vds.log", "prepare_base_vds")
    if rc:
        sys.exit("prepare_base_vds failed")

    # ── Step 3 – per-model work in parallel ───────────────────────────────────
    models = pd.read_csv(cfg.models_csv)
    required_cols = {"id", "url", "phenotype"}
    if not required_cols.issubset(models.columns):
        sys.exit(f"{cfg.models_csv} must have columns {required_cols}")

    max_workers = int(os.getenv("PARALLEL_WORKERS", len(models)))
    log.info("Launching %d parallel model workers", max_workers)

    def _process_model(row):
        m_id, m_url, m_pheno = str(row["id"]), row["url"], row["phenotype"]
        step_tag = f"model_{m_id}"
        logs = {
            "proc": log_dir/f"{step_tag}_01_process.log",
            "anal": log_dir/f"{step_tag}_02_analyze.log",
        }
        model_out_dir = f"{cfg.gcs_outputs}/scores/{m_id}"
        csv_gcs = f"{model_out_dir}/score_csv/{m_id}_scores.csv"
        ht_gcs = f"{model_out_dir}/hail_table/{m_id}_scores.ht"

        # process_prs_model.py
        rc1 = _run_subprocess(
            [
                PYTHON_EXECUTABLE, cfg.script_dir/"src"/"process_prs_model.py",
                "--prs_id", m_id,
                "--prs_url", m_url,
                "--base_cohort_vds_path", base_vds_gcs,
                "--gcs_temp_dir", cfg.gcs_intermediate,
                "--gcs_hail_temp_dir", cfg.gcs_hail_tmp,
                "--run_timestamp", cfg.ts,
                "--output_final_hail_table_gcs_path", ht_gcs,
                "--output_final_score_csv_gcs_path", csv_gcs,
                "--google_billing_project", cfg.env["GOOGLE_PROJECT"],
                "--spark_configurations_json", cfg.spark_conf_json,
                # Instruct Hail to run on the Dataproc YARN cluster
                "--hail_cluster_mode", "dataproc_yarn"
            ],
            env, logs["proc"], f"{m_id}:process"
        )
        if rc1:
            log.error("%s failed (processing)", m_id)
            return None

        # analyze_prs_results.py
        rc2 = _run_subprocess(
            [
                PYTHON_EXECUTABLE, cfg.script_dir/"src"/"analyze_prs_results.py",
                "--prs_id", m_id,
                "--prs_phenotype_label", m_pheno,
                "--score_csv_gcs_path", csv_gcs,
                "--wgs_ehr_ids_gcs_path", wgs_ehr_ids_gcs,
                "--phenotype_cases_csv_gcs_path", pheno_csv_gcs,
                "--phenotype_name", cfg.pheno_name,
                "--gcs_base_output_dir_run", cfg.gcs_outputs,
                "--run_timestamp", cfg.ts,
                "--google_billing_project", cfg.env["GOOGLE_PROJECT"],
                "--output_summary_file_name", f"{m_id}_summary.txt"
            ],
            env, logs["anal"], f"{m_id}:analyze"
        )
        if rc2:
            log.error("%s failed (analysis)", m_id)
            return None
        return logs["anal"]

    summaries: list[Path] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as ex:
        fut_to_id = {ex.submit(_process_model, row): row["id"] for _, row in models.iterrows()}
        for fut in concurrent.futures.as_completed(fut_to_id):
            mid = fut_to_id[fut]
            try:
                result = fut.result()
                if result:
                    summaries.append(result)
                    log.info("%s complete", mid)
                else:
                    log.warning("%s produced no summary (see logs)", mid)
            except Exception as e:
                log.exception("Worker for %s raised exception: %s", mid, e)

    # ── Step 4 – concatenate summaries ────────────────────────────────────────
    if summaries:
        concat_path = work_dir/f"all_models_summary_{cfg.ts}.txt"
        with concat_path.open("w") as out:
            for i, p in enumerate(summaries):
                if i:
                    out.write("\n" + "="*120 + "\n\n")
                out.write(p.read_text())
        log.info("Summaries combined → %s", concat_path)
    else:
        log.warning("No summaries generated!")

    log.info("Pipeline finished.  Outputs root: %s", cfg.gcs_outputs)


# ────────────────────────────────────────────────────────────────────────────────
#  CLI
# ────────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AoU PRS orchestrator")
    parser.add_argument("--config", default="config.yaml", help="Path to YAML config")
    args = parser.parse_args()
    cfg = Config(args.config)
    main(cfg)
