#!/usr/bin/env bash
# =====================================================================
# train_structmummi.sh
# Runs the StructMuMMI training in one command:
#   1. training/execute_johnsen.sh       train the XGBoost models (repeat_1 ... repeat_N)
#   2. training/run_prioritize_snps.sh   meta-model + SNP prioritization + VEP annotation
# =====================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"   # .../structmummi
TRAIN_DIR="${SCRIPT_DIR}/training"

start_from=1

# Options are collected per script:
TRAIN_ARGS=()     # for execute_johnsen.sh
PRIO_ARGS=()      # for run_prioritize_snps.sh

usage() {
cat <<EOF
Usage: train_structmummi.sh [options]

Pipeline:
  --start_from N          1 = training + prioritization (default)
                          2 = prioritization only (models already trained)

Summary option at each steps:
(/!/ Description of the option are given below)

Options given to BOTH steps (must be the same for training and prioritization):
  --genotype_path, --outdir, --device, --test_size, --rd_state_test, 
  --disease_name 
 
Options given to step 1 only (to train the model):
  --n_repeats, --continue, --seed_mode, --seeds, --n_jobs, --iteration_evals, --nsplits_shap_rank,
  --train_build_size, --min_k, --max_k
 
Options given to step 2 only (to prioritize SNPs):
  --ref_hg_target, --n_jobs_prio, --gpu_interaction, --inter_batch_size, --consensus_repeats, 
  --top_nb_snps, --ntop_inter1, --ntop_inter2
 
EOF
    echo
    echo "------------- Step 1: training model ----------------------"
    bash "${TRAIN_DIR}/execute_johnsen.sh" --help
    echo
    echo "------------- Step 2: Prioritizing SNPs (run_prioritize_snps.sh) -------------"
    bash "${TRAIN_DIR}/run_prioritize_snps.sh" --help
}


while [[ $# -gt 0 ]]; do
  case "$1" in
    # --- both steps
    --genotype_path|--outdir|--device|--test_size|--rd_state_test|--disease_name)
        TRAIN_ARGS+=( "$1" "${2:-}" ); PRIO_ARGS+=( "$1" "${2:-}" ); shift 2 ;;
    # --- training only
    --n_repeats|--seed_mode|--seeds|--n_jobs|--iteration_evals|--train_build_size|--min_k|--max_k|--continue)
        TRAIN_ARGS+=( "$1" "${2:-}" ); shift 2 ;;
    # --- prioritization only
    --ref_hg_target|--n_jobs_prio|--gpu_interaction|--inter_batch_size|--consensus_repeats|--top_nb_snps|--ntop_inter1|--ntop_inter2)
        PRIO_ARGS+=( "$1" "${2:-}" ); shift 2 ;;
    # --- pipeline
    --start_from)
        start_from="${2:-}"; shift 2 ;;
    -h|--help)
        usage; exit 0 ;;
    *)
        echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

[[ "$start_from" == "1" || "$start_from" == "2" ]] || { echo "ERROR: --start_from must be 1 or 2" >&2; exit 1; }
for s in execute_johnsen.sh run_prioritize_snps.sh; do
    [[ -f "${TRAIN_DIR}/${s}" ]] || { echo "ERROR: ${TRAIN_DIR}/${s} not found" >&2; exit 1; }
done

# --- Keep a copy of everything printed on screen in a log file
log_file="train_structmummi.log"
exec > >(tee -a "$log_file") 2>&1

# --- Run one step, print its name and duration
n_steps=2
run_step() {
    local n="$1"; shift
    local name="$1"; shift
    if (( n < start_from )); then
        echo ">>> Step $n/$n_steps skipped ($name)"
        return 0
    fi
    echo
    echo "=================================================================="
    echo ">>> Step $n/$n_steps: $name   [$(date '+%Y-%m-%d %H:%M:%S')]"
    echo "=================================================================="
    local t0=$SECONDS
    if ! "$@"; then
        echo
        echo "!!! Step $n ($name) FAILED. Fix the problem, then resume with:" >&2
        echo "!!!   train_structmummi.sh --start_from $n [same options as before]" >&2
        exit 1
    fi
    echo ">>> Step $n done in $(( (SECONDS - t0) / 60 )) min"
}

echo "train_structmummi started $(date '+%Y-%m-%d %H:%M:%S') in $(pwd)"

run_step 1 execute_johnsen.sh     bash "${TRAIN_DIR}/execute_johnsen.sh"     "${TRAIN_ARGS[@]}"
run_step 2 run_prioritize_snps.sh bash "${TRAIN_DIR}/run_prioritize_snps.sh" "${PRIO_ARGS[@]}"

echo
echo "train_structmummi finished $(date '+%Y-%m-%d %H:%M:%S')"
echo "Log: $(pwd)/$log_file"
