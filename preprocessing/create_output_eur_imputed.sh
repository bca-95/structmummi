#!/bin/bash
set -euo pipefail

export TMPDIR="${TMPDIR:-/tmp}"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$TMPDIR"

# --- Find conda automatically
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found in PATH" >&2; exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
 
# conda's activation scripts use unset variables, which breaks "set -u"
set +u
conda deactivate
conda activate gcta_env
set -u

# For the imputed data, I can only prepare the european sample
#as the imputed variants were only computed on the european
current_dir=$(pwd)
output_dir=""
perform_imputation="yes"

# - For parsing argument - 
# Bash giveds to the script a list of words, where each word is considered
# as argument --> for instance : -i (arg 1) input.txt (arg 2) -o ...
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in

    --perform_imputation)
      perform_imputation="$2"
      shift 2
      ;;

    *)
      echo "Unknown option: $1"
      echo "Usage: /pathtoscript/create_output_eur_imputed.sh --perform_imputation <yes or no, default: yes>"
      exit 1
      ;;
  esac
done



if [[ "$perform_imputation" == "yes" ]]; then
    output_dir="qc_imputed/outputs_eur_imputed"
    mkdir -p $output_dir
    cd $output_dir
elif [[ "$perform_imputation" == "no" ]]; then 
    output_dir="output_eur"
    mkdir -p $output_dir
    cd $output_dir
fi

#cp ../../outputs_european/plink_12* .
#cp ../plink_12* .
cp "${current_dir}/qc_imputed/"plink_12.* .

## Generate matrix
plink2 --bfile plink_12 --recode A --out merge_recodeA

## delete intermediate file
## /!\ DON'T DELETE PHASING DATA
## /1\ DON'T DELETE POST-IMPUTATION DATA
## They take times to generate

find "$current_dir" -maxdepth 1 -type f \
  ! -name "*.png" \
  ! -name "*.pdf" \
  ! -name "*.sh" \
  ! -name "*.tped*" \
  ! -name "*.tfam*" \
  ! -name "*.prune.in" \
  ! -name ".localhistory" \
  ! -name "plink_imp2_*" \
  ! -name "plink_imp4_*" \
   -delete

find "$current_dir/qc_imputed" -maxdepth 1 -type f \
  ! -name "*.png" \
  ! -name "*.pdf" \
  ! -name "*.sh" \
  ! -name "*.tped*" \
  ! -name "*.tfam*" \
  ! -name "*.prune.in" \
  ! -name ".localhistory" \
  ! -name "plink_imp2_*" \
  ! -name "plink_imp4_*" \
   -delete


