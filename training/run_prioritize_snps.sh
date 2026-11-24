#!/usr/bin/env bash
# -----------------------------------------------------
# This script execute the prioritization snps
# -----------------------------------------------------
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
REPO_DIR="$(dirname "$SCRIPT_DIR")"
PYTHON_SCRIPT="${SCRIPT_DIR}/run_prioritize_snps.py"

export TMPDIR="${TMPDIR:-/tmp}"
mkdir -p "$TMPDIR"
PYTHON_BIN="python3"  


# --------------------------------------------------------
# config template
# --------------------------------------------------------

GENOTYPE_PATH="qc_imputed/outputs_eur_imputed/merge_recodeA.raw"
OUTDIR="training_analysis/ML_johnsen"
DISEASE_NAME="diseaseX"
REF_HG_TARGET="hg38"
N_JOBS_PRIO="46"
GPU_INTERACTION="yes"
INTER_BATCH_SIZE=5 #(For 40, it takes 86GB. For 10, 21GB, and I only have 16GB)

FEATURES_FILE_NAME="top_ML_feature.npy"
MODEL_FILE="model_final.json"
BEST_TUNED_PARAMS_FILE="model_bestparams_build.json"

FASTA_FILE="${REPO_DIR}/references/reference_p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
VEP_DIR="${REPO_DIR}/external_tools/.vep/"

TOP_NB_SNPS=30
NTOP_INTER1=10
NTOP_INTER2=20

CONSENSUS_REPEATS="auto"
TEST_SIZE=0.2
RD_STATE_TEST=42          # kept FIXED across all repeats -- same held-out test set every time
DEVICE="cuda:0"


usage() { 
cat <<EOF
Usage: run_training.sh [options]
 
Input / output:
  --genotype_path PATH         Genotype matrix from plink --recode A (.raw)
                          (default: $GENOTYPE_PATH)
  --outdir DIR            Output folder (default: ./$OUTDIR)
  --disease_name          Disease name, e.g. "CAD" for coronary artery disease (default: $DISEASE_NAME )
  --ref_hg_target BUILD   Genome build of the data, hg19 or hg38, for VEP (default: $REF_HG_TARGET)

  --n_jobs_prio           Number of threads to compute prioritization on CPU (default: $N_JOBS_PRIO)
  --gpu_interaction bool  Compute SHAP interaction with GPU (faster but can lead to memory crash), default: $GPU_INTERACTION
  --inter_batch_size N    Compute SHAP interaction per batch (GPU). The number of batch does not influence the output values.
                          SHAP interaction is computed per samples. Default: $INTER_BATCH_SIZE

Prioritization:
  --top_nb_snps N         Number of top SNPs based on the SHAP values consensus between the repeat models (default: $TOP_NB_SNPS)
  --ntop_inter1 N         Number of top interacting SNPs based on SHAP interaction values, interaction consensus between the repeat models (default: $NTOP_INTER1)
  --ntop_inter2 N         Number of top SNPs that are interacting to the top_nb_snps (default: $NTOP_INTER2)

Repeats used in the consensus:
  --consensus_repeats R   auto   : all finished repeats, from repeat_1 (default)
                          3      : repeat_1 to repeat_3
                          4-10   : repeat_4 to repeat_10
                          Results go to <outdir>/consensus_<repeats>/, e.g. ML_johnsen/consensus_1-3/ 
Model:
  --device DEV            cuda:0 (GPU) or cpu (default: $DEVICE)
  --test_size F           Fraction of samples in the test set (default: $TEST_SIZE). /!/ Need to be same value as used during training
  --rd_state_test F       Assigned seed for partitioning the test set (default: $RD_STATE_TEST)  /!/ Need to be same value as used during training
EOF
}


while [[ $# -gt 0 ]]; do
  case "$1" in
    --genotype_path)        GENOTYPE_PATH="${2:-}";  shift 2 ;;
    --outdir)               OUTDIR="${2:-}";         shift 2 ;;
    --n_jobs_prio)          N_JOBS_PRIO="${2:-}";         shift 2 ;;
    --gpu_interaction)      GPU_INTERACTION="${2:-}";         shift 2 ;;
    --inter_batch_size)     INTER_BATCH_SIZE="${2:-}";         shift 2 ;;
    --disease_name)         DISEASE_NAME="${2:-}";   shift 2 ;;
    --ref_hg_target)        REF_HG_TARGET="${2:-}";  shift 2 ;;
    --consensus_repeats)    CONSENSUS_REPEATS="${2:-}";     shift 2 ;;
    --top_nb_snps)       TOP_NB_SNPS="${2:-}";     shift 2 ;;
    --ntop_inter1)       NTOP_INTER1="${2:-}";     shift 2 ;;
    --ntop_inter2)       NTOP_INTER2="${2:-}";     shift 2 ;;
    --device)            DEVICE="${2:-}";         shift 2 ;;
    --test_size)         TEST_SIZE="${2:-}";      shift 2 ;;
    --rd_state_test)     RD_STATE_TEST="${2:-}";  shift 2 ;;
    -h|--help)           usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
  esac
done

# --------------------------------------------------------
# --------------------------------------------------------
# Checks
# --------------------------------------------------------
[[ -f "$PYTHON_SCRIPT" ]] || { echo "ERROR: $PYTHON_SCRIPT not found" >&2; exit 1; }
[[ -f "$GENOTYPE_PATH" ]] || { echo "ERROR: genotype matrix not found: $GENOTYPE_PATH (use --genotype)" >&2; exit 1; }
[[ -d "$OUTDIR" ]]        || { echo "ERROR: $OUTDIR not found: run execute_johnsen.sh first, or give --outdir" >&2; exit 1; }
 
# Reference genome for VEP, depending on the build
if [[ "$REF_HG_TARGET" == "hg19" ]]; then
    FASTA_FILE="${REPO_DIR}/references/reference_hg19/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
    ASSEMBLY="GRCh37"
elif [[ "$REF_HG_TARGET" == "hg38" ]]; then
    FASTA_FILE="${REPO_DIR}/references/reference_p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    ASSEMBLY="GRCh38"
else
    echo "ERROR: --ref_hg_target must be hg19 or hg38 (got: '$REF_HG_TARGET')" >&2; exit 1
fi
 
# The paths go into the config file: make them absolute
GENOTYPE_PATH="$(readlink -f "$GENOTYPE_PATH")"
OUTDIR="$(cd "$OUTDIR" && pwd)"

# --------------------------------------------------------
# Repeats used in the consensus -> list REPEAT_IDS (e.g. 1 2 3)
# --------------------------------------------------------
REPEAT_IDS=()
if [[ "$CONSENSUS_REPEATS" == "auto" ]]; then
    # all finished repeats from repeat_1, without gaps
    r=1
    while [[ -f "${OUTDIR}/repeat_${r}/${MODEL_FILE}" ]]; do
        REPEAT_IDS+=( "$r" )
        r=$(( r + 1 ))
    done
    (( ${#REPEAT_IDS[@]} > 0 )) || { echo "ERROR: no finished repeat found in $OUTDIR (expected repeat_1/${MODEL_FILE})" >&2; exit 1; }
    label="${#REPEAT_IDS[@]}"                   # auto -> number of repeats, e.g. consensus_10
 
    # warn if finished repeats exist after a gap: they are left out
    for d in "$OUTDIR"/repeat_*; do
        num="${d##*repeat_}"
        if [[ "$num" =~ ^[0-9]+$ ]] && (( num > r )) && [[ -f "$d/${MODEL_FILE}" ]]; then
            echo "WARNING: repeat_${r} is not finished, so only repeat_1 ... repeat_$(( r - 1 )) are used." >&2
            echo "WARNING: later finished repeats (e.g. repeat_${num}) are left out. Use e.g. --consensus_repeats 1,2,4 to choose." >&2
            break
        fi
    done
 
elif [[ "$CONSENSUS_REPEATS" =~ ^[1-9][0-9]*$ ]]; then
    # "3" -> 1 2 3
    REPEAT_IDS=( $(seq 1 "$CONSENSUS_REPEATS") )
    label="$CONSENSUS_REPEATS"                  # 3 -> consensus_3
 
elif [[ "$CONSENSUS_REPEATS" =~ ^([1-9][0-9]*)-([1-9][0-9]*)$ ]]; then
    # "4-10" -> 4 5 ... 10   (BASH_REMATCH holds the two numbers matched in ( ))
    first="${BASH_REMATCH[1]}"; last="${BASH_REMATCH[2]}"
    (( first <= last )) || { echo "ERROR: --consensus_repeats $CONSENSUS_REPEATS: the first number must be <= the second" >&2; exit 1; }
    REPEAT_IDS=( $(seq "$first" "$last") )
    label="${first}-${last}"                    # 4-10 -> consensus_4-10
 
else
    echo "ERROR: --consensus_repeats must be auto, a number (3) or a range (4-10) (got '$CONSENSUS_REPEATS')" >&2
    exit 1
fi
 
(( ${#REPEAT_IDS[@]} >= 2 )) || { echo "ERROR: the consensus needs at least 2 repeats (got: ${REPEAT_IDS[*]})" >&2; exit 1; }
 
# Every repeat used must have its 3 files
for r in "${REPEAT_IDS[@]}"; do
    for f in "$MODEL_FILE" "$FEATURES_FILE_NAME" "$BEST_TUNED_PARAMS_FILE"; do
        [[ -f "${OUTDIR}/repeat_${r}/${f}" ]] || { echo "ERROR: missing ${OUTDIR}/repeat_${r}/${f} (did repeat_${r} finish?)" >&2; exit 1; }
    done
done
 
# Results folder (label was set above, with the repeats)
RESULT_OUT="${OUTDIR}/consensus_${label}"
 
# JSON list for the config: 1 2 3 -> [1, 2, 3]
REPEAT_IDS_JSON="[$(IFS=,; echo "${REPEAT_IDS[*]}")]"
 
# Folder for the selected SNPs (the Python script writes into it)
mkdir -p "${RESULT_OUT}/selected_snps"
 
config_path="${RESULT_OUT}/config_meanmodels.json" 

# --------------------------------------------------------
   
cat > "$config_path" << EOF
{

  "script": "${PYTHON_SCRIPT}",
  "nb_repeats": ${#REPEAT_IDS[@]},
  "repeat_ids": ${REPEAT_IDS_JSON},
  "n_jobs_prio": ${N_JOBS_PRIO},
  "gpu_interaction": "${GPU_INTERACTION}",
  "inter_batch_size": ${INTER_BATCH_SIZE},
  "vep_dir": "${VEP_DIR}",
  "assembly": "${ASSEMBLY}",
  "consensus_repeats": "${CONSENSUS_REPEATS}",
  "genotype_path": "${GENOTYPE_PATH}",
  "param_out": "${OUTDIR}/",
  "result_out": "${RESULT_OUT}/",
  "top_nb_snps": "${TOP_NB_SNPS}",
  "ntop_inter1": "${NTOP_INTER1}",
  "ntop_inter2": "${NTOP_INTER2}",
  "disease_name": "${DISEASE_NAME}",
  "features_file": "${FEATURES_FILE_NAME}",
  "model_file": "${MODEL_FILE}",
  "best_tuned_params_file": "${BEST_TUNED_PARAMS_FILE}",
  "fasta_file": "${FASTA_FILE}",
  "test_size": ${TEST_SIZE},
  "rd_state_test": ${RD_STATE_TEST},
  "device": "${DEVICE}"
}
EOF

echo "============================================================"
echo "Running selection of impactful SNPs"
echo "Models  : repeat_$(IFS=,; echo "${REPEAT_IDS[*]}") in ${OUTDIR}"
echo "Results : ${RESULT_OUT}"
echo "Config  : ${config_path}"
echo "============================================================"
 
"$PYTHON_BIN" "$PYTHON_SCRIPT" --config "$config_path"
 
echo " ... JOB COMPLETED ... "
echo " Results in: ${RESULT_OUT}/"


