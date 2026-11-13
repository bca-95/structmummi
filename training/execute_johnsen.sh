#!/usr/bin/env bash
# ----------------------------------------------
# Training the input datasets using the model based on Johnson et al.
# The input should be a genotype matrix
# If the preprocessing pipelines from StructMuMMI where performed, 
# simply run : execute_johnsen.sh

# Default input: the output of the preprocessing, i.e.
#   qc_imputed/outputs_eur_imputed/merge_recodeA.raw   (in the current folder)
# Any PLINK --recode A matrix (.raw) can be given with --genotype.
#
# Output: one folder per repeat, <outdir>/repeat_1 ... repeat_N
# A repeat that already has a final model is skipped (use --force to redo it)
# ----------------------------------------------
set -euo pipefail

# --- Find conda automatically
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found in PATH" >&2; exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

# conda's activation scripts use unset variables, which breaks "set -u"
set +u
conda deactivate
conda activate jupyter_env2
set -u

SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
PYTHON_SCRIPT="${SCRIPT_DIR}/xgboost_simple_johnsen.py"

export TMPDIR="${TMPDIR:-/tmp}"

PYTHON_BIN="python3"  

# --------------------------------------------------------
# config template
# --------------------------------------------------------
GENOTYPE_PATH="qc_imputed/outputs_eur_imputed/merge_recodeA.raw"
OUTDIR="ML_johnsen"
N_REPEATS=10
FIRST_REPEAT=1

#BASE_PARAM_OUT="/wibicomfs/STBS/carene/sanger/processed_data/merged_control_case/wtccc2_psoriasis/MAF_01/qc_imputed/outputs_eur_imputed/ML/model_johnsen/version3.0/"

TRAIN_BUILD_SIZE=0.4
TEST_SIZE=0.2
NSPLITS_SHAP_RANK=4
RD_STATE_TEST=42          # kept FIXED across all repeats -- same held-out test set every time
MAX_K=4500
MIN_K=1200
N_JOBS=46
DEVICE="cuda:0"
ITERATION_EVALS=200
NAME_IMAGE="Disease"
LOG_FILE="outlog_johnsen.txt"

# Define the seed to split the rank and build sets, assigned for each repeats
RD_STATE_BUILDS=(42 71 118 127 201 13 27 97 148 174)

# The seeds can be modified
SEED_MODE="default"
SEED_MIN=12
SEED_MAX=999
USER_SEEDS=""
CONTINUE="no"

usage() {
cat <<EOF
Usage: run_training.sh [options]
 
Input / output:
  --genotype PATH         Genotype matrix from plink --recode A (.raw)
                          (default: $GENOTYPE_PATH)
  --outdir DIR            Output folder (default: ./$OUTDIR)
  --name_image NAME             Label used in plot/file names (default: $NAME_IMAGE)
 
Repeats:
  --n_repeats N           Number of repeats. This creates directory repeat_1 ... repeat_N. (default: $N_REPEATS)
  --continue              Keep the repeats already done in --outdir and add N new ones
                          after the last finished repeat (e.g. repeat_11 ... repeat_15)
                          Without it, repeats start at repeat_1 and old results are overwritten 

Seeds (split of the samples into rank/build sets assigned for each repeat):
  --seed_mode MODE        default   : fixed list (${RD_STATE_BUILDS[*]}) for repeats 1-10,
                                      random seeds for repeats above 10
                          random    : random seeds between $SEED_MIN and $SEED_MAX
                          user_seed : your own seeds, given with --seeds
                          (default: $SEED_MODE)
  --seeds "S1,S2,..."     With --seed_mode user_seed: one seed per repeat, in order
  The seeds used are written to <outdir>/seeds_used.tsv
 
 
Model:
  --device DEV            cuda:0 (GPU) or cpu (default: $DEVICE)
  --n_jobs N              CPU threads (default: $N_JOBS)
  --iteration_evals N     Hyper-parameter tuning trials (default: $ITERATION_EVALS)
  --test_size F           Fraction of samples in the test set (default: $TEST_SIZE)
  --train_build_size F    Fraction of the remaining samples in the build set (default: $TRAIN_BUILD_SIZE)
  --min_k F               Min number of top features tried (default: $MIN_K)
  --max_k F               Max number of top features tried (default: $MAX_K)
  --nsplits_shap_rank N   When computing feature importances, the rank datasets is split into (default: $NSPLITS_SHAP_RANK)
  --rd_state_test F       Assigned seed for partitioning the test set (default: $RD_STATE_TEST)
EOF
}


while [[ $# -gt 0 ]]; do
  case "$1" in
    --genotype_path)          GENOTYPE_PATH="${2:-}";    shift 2 ;;
    --outdir)                 OUTDIR="${2:-}";           shift 2 ;;
    --name_image)             NAME_IMAGE="${2:-}";       shift 2 ;;
    --n_repeats)         N_REPEATS="${2:-}";        shift 2 ;;
    --first_repeat)      FIRST_REPEAT="${2:-}";     shift 2 ;;
    --continue)          CONTINUE="yes";            shift 1 ;;
    --seed_mode)         SEED_MODE="${2:-}";        shift 2 ;;
    --seeds)             USER_SEEDS="${2:-}";       shift 2 ;;
    --device)            DEVICE="${2:-}";           shift 2 ;;
    --n_jobs)            N_JOBS="${2:-}";           shift 2 ;;
    --iteration_evals)   ITERATION_EVALS="${2:-}";  shift 2 ;;
    --test_size)         TEST_SIZE="${2:-}";        shift 2 ;;
    --train_build_size)  TRAIN_BUILD_SIZE="${2:-}"; shift 2 ;;
    --min_k)             MIN_K="${2:-}";       shift 2 ;;
    --max_k)             MAX_K="${2:-}";       shift 2 ;;
    -h|--help)           usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done


# --------------------------------------------------------
# Checks
# --------------------------------------------------------
[[ -f "$PYTHON_SCRIPT" ]] || { echo "ERROR: $PYTHON_SCRIPT not found" >&2; exit 1; }
[[ -f "$GENOTYPE_PATH" ]] || { echo "ERROR: genotype matrix not found: $GENOTYPE_PATH (If not yet done, create a genotype matrix for instance as follow \
plink2 --bfile <plink_datasets_format> --recode A --out <genotype_diseasesX> )" >&2; exit 1; }
[[ "$N_REPEATS" =~ ^[1-9][0-9]*$ ]]    || { echo "ERROR: --n_repeats must be a positive number" >&2; exit 1; }


# --------------------------------------------------------

GENOTYPE_PATH="$(readlink -f "$GENOTYPE_PATH")"

mkdir -p "$OUTDIR"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# --- Log of the seeds used (one line per repeat run, appended across runs)
SEED_LOG="${OUTDIR}/seeds_used.tsv"

# --------------------------------------------------------
# Which repeats to run: START_REPEAT ... END_REPEAT
#   normal run   : repeat_1 ... repeat_N
#   --continue   : starts after the last finished repeat in OUTDIR
#                  (finished = its folder contains model_final_johnsen.json)
# --------------------------------------------------------
START_REPEAT=1
if [[ "$CONTINUE" == "yes" ]]; then
    last_done=0
    for d in "$OUTDIR"/repeat_*; do
        [[ -f "$d/model_final_johnsen.json" ]] || continue
        num="${d##*repeat_}"                      # ".../repeat_7" -> "7"
        if [[ "$num" =~ ^[0-9]+$ ]] && (( num > last_done )); then
            last_done=$num
        fi
    done
    START_REPEAT=$(( last_done + 1 ))
fi
END_REPEAT=$(( START_REPEAT + N_REPEATS - 1 )) 



# --------------------------------------------------------
# Seeds: build the list SEEDS_RUN, one seed per repeat that will run (in order)
# --------------------------------------------------------
# Seeds that random draws must avoid: the fixed list + every seed already used in OUTDIR
EXCLUDE_SEEDS=( "${RD_STATE_BUILDS[@]}" )
if [[ -f "$SEED_LOG" ]]; then
    mapfile -t used < <(tail -n +2 "$SEED_LOG" | cut -f2)
    EXCLUDE_SEEDS+=( "${used[@]}" )
fi
 
# draw N random seeds between SEED_MIN and SEED_MAX, all different,
# and not in RD_STATE_BUILDS 
# This draw n random seeds that are not present in the default seed RD_STATE_BUILDS
draw_random_seeds() {
    local n="$1"
    seq "$SEED_MIN" "$SEED_MAX" | grep -vxF -f <(printf '%s\n' "${EXCLUDE_SEEDS[@]}") | shuf -n "$n"
}

n_available=$(seq "$SEED_MIN" "$SEED_MAX" | grep -vxF -f <(printf '%s\n' "${EXCLUDE_SEEDS[@]}") | wc -l)

# Here we have 3 modes, default, random and user_seed
# case is like multiple if in a row
SEEDS_RUN=()
case "$SEED_MODE" in
 
    default)
        #-z means "is empty". If someone gave --seeds without choosing user_seed
        [[ -z "$USER_SEEDS" ]] || { echo "ERROR: --seeds is only used with --seed_mode user_seed" >&2; exit 1; }
 
        # repeat_1 ... repeat_10 use the fixed list, repeats above 10 get random seeds
        n_extra=0
        for (( r = START_REPEAT; r <= END_REPEAT; r++ )); do
            if (( r > ${#RD_STATE_BUILDS[@]} )); then
                n_extra=$(( n_extra + 1 ))
            fi
        done
        (( n_extra <= n_available )) || { echo "ERROR: too many repeats for the random seed range" >&2; exit 1; }
        extra=()
        if (( n_extra > 0 )); then
            mapfile -t extra < <(draw_random_seeds "$n_extra")
        fi
        i_extra=0
        for (( r = START_REPEAT; r <= END_REPEAT; r++ )); do
            if (( r <= ${#RD_STATE_BUILDS[@]} )); then
                SEEDS_RUN+=( "${RD_STATE_BUILDS[$((r - 1))]}" )
            else
                SEEDS_RUN+=( "${extra[$i_extra]}" )
                i_extra=$(( i_extra + 1 ))
            fi
        done ;;
 
    random)
        [[ -z "$USER_SEEDS" ]] || { echo "ERROR: --seeds is only used with --seed_mode user_seed" >&2; exit 1; }
        (( N_REPEATS <= n_available )) || { echo "ERROR: --n_repeats is larger than the number of possible seeds ($n_available)" >&2; exit 1; }
 
        # This draws all the seeds at random, one per repeat
        mapfile -t SEEDS_RUN < <(draw_random_seeds "$N_REPEATS") ;;
 
    user_seed)
        #-n means "is not empty"
        [[ -n "$USER_SEEDS" ]] || { echo "ERROR: --seed_mode user_seed needs --seeds \"S1,S2,...\"" >&2; exit 1; }
        IFS=',' read -r -a SEEDS_RUN <<< "${USER_SEEDS// /}"    # "11, 22,33" -> (11 22 33)
        for sd in "${SEEDS_RUN[@]}"; do
            [[ "$sd" =~ ^[0-9]+$ ]] || { echo "ERROR: --seeds must be whole numbers separated by commas (got '$sd')" >&2; exit 1; }
        done
        (( ${#SEEDS_RUN[@]} == N_REPEATS )) || { echo "ERROR: --seeds gives ${#SEEDS_RUN[@]} seeds but --n_repeats is $N_REPEATS: give one seed per repeat" >&2; exit 1; } ;;
    *)
        echo "ERROR: --seed_mode must be default, random or user_seed (got '$SEED_MODE')" >&2; exit 1 ;;
esac

 
# -------------------------------------------------------
echo "Genotype matrix : $GENOTYPE_PATH"
echo "Output folder   : $OUTDIR"
echo "Repeats         : $FIRST_REPEAT to $(( FIRST_REPEAT + N_REPEATS - 1 ))"
echo "Device          : $DEVICE"
echo "Seed mode       : $SEED_MODE"
echo "Seeds           : ${SEEDS_RUN[*]}"

# --- Log of the seeds used (one line per repeat run, appended across runs)
[[ -f "$SEED_LOG" ]] || printf "repeat_id\trd_state_build\tseed_mode\trd_state_test\tdate\n" > "$SEED_LOG"

# ------------------------------------------------------------
for (( r = START_REPEAT; r <= END_REPEAT; r++ )); do
#for (( r = 1; r <= N_REPEATS; r++ )); do
    repeat_id="repeat_${r}"
    rd_state_build="${SEEDS_RUN[$((r - 1))]}"

    outdir="${OUTDIR}/${repeat_id}"
    mkdir -p "$outdir"

    config_path="${outdir}/config_johnsen_${repeat_id}.json"
    cat > "$config_path" << EOF
{
"repeat_id": "${repeat_id}",
"genotype_path": "${GENOTYPE_PATH}",
"param_out": "${OUTDIR}/",
"train_build_size": ${TRAIN_BUILD_SIZE},
"test_size": ${TEST_SIZE},
"nsplits_shap_rank": ${NSPLITS_SHAP_RANK},
"rd_state_test": ${RD_STATE_TEST},
"rd_state_build": ${rd_state_build},
"n_jobs": ${N_JOBS},
"device": "${DEVICE}",
"iteration_evals": ${ITERATION_EVALS},
"name_image": "${NAME_IMAGE}",
"log_file": "${LOG_FILE}",
"max_k": ${MAX_K},
"min_k": ${MIN_K}
}
EOF

    echo "============================================================"
    echo " Running repead_id : ${repeat_id} --  with seed (rd_state_build=${rd_state_build})"
    echo " Config: ${config_path}"
    echo "============================================================"

    "$PYTHON_BIN" "$PYTHON_SCRIPT" --config "$config_path"

    echo " Finished ${repeat_id}"
done
