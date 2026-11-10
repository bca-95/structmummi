#!/bin/bash

# This script lauches the pre-processing pipelines
# This script takes as inputs raw datasets of cases and control datasets.
# Proceed to harmonization/standardization, stratification, imputation, QC and prepare
# the final ouptut for training.

# The preprocessing pipelines can be started from later steps if needed more manipulation, 
# such as european delimitation during stratification, MAF threshold,... without repeating the
# whole process.
# However, the previous step should have been performed at list once, to create the ouptuts files
# where the name of these files can not be manipulated.
set -euo pipefail

export TMPDIR="${TMPDIR:-/tmp}"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$TMPDIR"

# --- Locate the repo, wherever it was cloned
REPO_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)" # ../structmummi (different from the base script)
SCRIPT_DIR="${REPO_DIR}/preprocessing"

# --- Default variables
input_control=""
input_case=""
ref_hg_target=""
ref_fasta=""
convert_to_rsid="yes"
liftover="yes"

eur_marginX="0.3"
eur_marginY="0.3"

cpu_threads_phasing="30"
cpu_threads_imputation="16"
hc_threshold="0.1"

maf="0.01"
perform_imputation="yes"
cleanup="yes"
start_from=1


usage() {
cat <<EOF
Usage: pre_structmummi.sh --input_control <plink_basename> --input_case <plink_basename> --ref_hg_target <hg19|hg38> [options]
 
Required:
  --input_control           PLINK basename of the control dataset (without .bed/.bim/.fam)
  --input_case              PLINK basename of the case dataset
  --ref_hg_target           Target genome build: hg19 or hg38
 
Options:
  --perform_imputation      yes|no: run phasing + imputation (default: $perform_imputation)
  --eur_marginX             Margin around EUR samples on PC1, fraction of EUR range (default: $eur_marginX)
  --eur_marginY             Margin around EUR samples on PC2, fraction of EUR range (default: $eur_marginY)
  --cpu_threads_phasing     Threads per chromosome for Eagle (default: $cpu_threads_phasing)
  --cpu_threads_imputation  Threads per chromosome for Minimac4 (default: $cpu_threads_imputation)
  --hc_threshold            Hard-call threshold after imputation (default: $hc_threshold)
  --maf                     Minor allele frequency threshold (default: $maf)
  --cleanup                 yes|no: delete intermediate files at the end (default: $cleanup)
  --start_from              Step to start from, 1-5, to resume after a failure (default: 1)
EOF
}


# Parsing argument - 
# Bash giveds to the script a list of words, where each word is considered
# as argument --> for instance : -i (arg 1) input.txt (arg 2) -o ...
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --input_control)
      input_control="$2"
      shift 2
      ;;
    --input_case)
      input_case="$2"
      shift 2
      ;;

    --convert_to_rsid)
      convert_to_rsid="$2"
      shift 2
      ;;

    --liftover)
      liftover="$2"
      shift 2
      ;;

    --ref_hg_target)
      ref_hg_target="$2"
      shift 2
      ;;

    --eur_marginX)
      eur_marginX="$2"
      shift 2
      ;;

    --eur_marginY)
      eur_marginY="$2"
      shift 2
      ;;

    --cpu_threads_phasing)
      cpu_threads_phasing="$2"
      shift 2
      ;;

    --cpu_threads_imputation)
      cpu_threads_imputation="$2"
      shift 2
      ;;

    --hc_threshold)
      hc_threshold="$2"
      shift 2
      ;;

    --maf)
      maf="$2"
      shift 2
      ;;

    --perform_imputation)
      perform_imputation="$2"
      shift 2
      ;;

    --cleanup)
      cleanup="$2"
      shift 2
      ;;

    --start_from)
      start_from="$2"
      shift 2
      ;;

    -h|--help)                usage; exit 0 ;;

    *)
      echo "Unknown option: $1" >&2; usage >&2; exit 1 ;;
#      echo "Usage: pre_structmummi.sh --input_control plink_basename (required, eg. WTCCC_control (without the extension .bim .bed .fam))
#                                      --input_case plink_basename (required, eg. WTCCC_diabetes1 (without the extension .bim .bed .fam))
#                                      --convert_to_rsid <yes or no (optional, default: yes) #the variants name is standardized to chr:pos:ref:alt.
#                                      --liftover <yes or no (optional, default: yes) #position is corrected using reference hg and allele harmonized using reference panel"
#      exit 1
#      ;;
  esac
done

# --- Check the arguments
[[ "$ref_hg_target" == "hg19" || "$ref_hg_target" == "hg38" ]] || { echo "ERROR: --ref_hg_target must be hg19 or hg38" >&2; usage >&2; exit 1; }
[[ "$perform_imputation" == "yes" || "$perform_imputation" == "no" ]] || { echo "ERROR: --perform_imputation must be yes or no" >&2; exit 1; }
[[ "$start_from" =~ ^[1-6]$ ]] || { echo "ERROR: --start_from must be a number from 1 to 6" >&2; exit 1; }
if (( start_from == 1 )); then
    [[ -n "$input_control" && -n "$input_case" ]] || { echo "ERROR: --input_control and --input_case are required" >&2; usage >&2; exit 1; }
    for f in "${input_control}.bed" "${input_case}.bed"; do
        [[ -f "$f" ]] || { echo "ERROR: $f not found" >&2; exit 1; }
    done
fi


# --- Keep a copy of everything printed on screen in a log file
log_file="pre_structmummi.log"
exec > >(tee -a "$log_file") 2>&1

# --- Helper: run one step, print its name and duration
run_step() {
    local n="$1"; shift
    local name="$1"; shift
    if (( n < start_from )); then
        echo ">>> Step $n/6 skipped ($name)"
        return 0
    fi
    echo
    echo "=================================================================="
    echo ">>> Step $n/6: $name   [$(date '+%Y-%m-%d %H:%M:%S')]"
    echo "=================================================================="
    local t0=$SECONDS
    if ! bash "${SCRIPT_DIR}/${name}" "$@"; then
        echo
        echo "!!! Step $n ($name) FAILED. Fix the problem, then resume with:" >&2
        echo "!!!   pre_structmummi.sh --ref_hg_target $ref_hg_target --start_from $n [same options as before]" >&2
        exit 1
    fi
    echo ">>> Step $n done in $(( (SECONDS - t0) / 60 )) min"
}

# step 1: Merge the control and case datasets
# This command also harmonized the variants coordinates to 1000G
run_step 1 run_merge.sh \
    --input_control $input_control \
    --input_case    $input_case  \
    --ref_hg_target $ref_hg_target \
    --convert_to_rsid $convert_to_rsid \
    --liftover $liftover \


# step2: This script aim to extract only european samples
# using the reference panel 1000G
run_step 2 simplified_run_stratification.sh \
    --ref_hg_target $ref_hg_target \
    --eur_marginX $eur_marginX \
    --eur_marginY $eur_marginY

# step 3 and 4: This script perform imputation
# Then, a simple quality control is performed on the imputed variants
# This steps take a lot of time
if [[ "$perform_imputation" == "yes" ]]; then
    run_step 3 run_imputation.sh \
            --ref_hg_target "$ref_hg_target" \
            --cpu_threads_phasing "$cpu_threads_phasing" \
            --cpu_threads_imputation "$cpu_threads_imputation"
    
    run_step 4 run_after_imputation.sh \
            --ref_hg_target "$ref_hg_target" \
            --hc_threshold "$hc_threshold"
else
    echo "NO imputation is performed"
fi

# step 5 : Perform quality control
run_step 5 qc_filtering_snpid_imputed.sh \
    --ref_hg_target "$ref_hg_target" \
    --maf "$maf" \
    --perform_imputation "$perform_imputation"
 
# step 6: Delete intermediate files
# If imputation performed the phasing and imputation file is not deleted
# and can be reused. They take a large amount of time to generated.
run_step 6 create_output_eur_imputed.sh \
    --perform_imputation "$perform_imputation" \
    --cleanup "$cleanup"

