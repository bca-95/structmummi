#!/bin/bash
set -euo pipefail

export TMPDIR=/scratch/SCRATCH_SAS/carene/tmp

# --- Locate the repo, wherever it was cloned
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)" # ../structmummi/preprocessing
REPO_DIR="$(dirname "$SCRIPT_DIR")" # dirname remove the last part of the path # ../structmummi

# --- Find conda automatically
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found in PATH" >&2; exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"
 
# conda's activation scripts use unset variables, which breaks "set -u"
set +u
conda deactivate
conda activate liftover_env
set -u

# Use compiled bcftools, the conda version raise error for the lifted data
# hwv, the conda version does not bring any error on 1000G data
export BCFTOOLS="${REPO_DIR}/external_tools/bcftools-1.23.1/bcftools"
export BCFTOOLS_PLUGINS="${REPO_DIR}/external_tools/bcftools-1.23.1/plugins"

plink2="${REPO_DIR}/external_tools/plink/plink2"


# This script aims to liftover from hg36 to hg37
# The WTCCC2 data were aligned based on hg36
# For the purpose of population stratification and genes mapping, WTCCC2 will be liftover to hg37

# For the liftover with bcftools
# the chain file name the chromosome with chr, 
# for easier data reading, the reference should be change to chr as well
# thus the input for liftover

### --> By using liftover, only the position of the variants is modified, the alleles remains the same
# this leads to mismatch in further step like during the imputation steps
# Therefore, in the subsequent steps, we performed allele harmonization.


# Set the variable, so input can take case or control
input_=""
output_=""
ref_hg_target=""
CHAIN=""
ref_seq=""

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --input_)
      input_="$2"
      shift 2
      ;;
    --output_)
      output_="$2"
      shift 2
      ;;
    --ref_hg_target)
      ref_hg_target="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: /pathtoscript/run_liftover.sh --input_ plink_basename --output_ basename_for_output"
      exit 1
      ;;
  esac
done

if [[ "$ref_hg_target" == "hg19" ]]; then
    echo " ---- Lifting from hg18 to hg19 ----"
    CHAIN="${SCRIPT_DIR}/structmummi/external_tools/liftover/hg18ToHg19.over.chain.gz"
    ref_seq="/wibicomfs/STBS/carene/tuto_gpt/tuto2/reference_hg19/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
elif [[ "$ref_hg_target" == "hg38" ]]; then
    echo " ---- Lifting from hg18 to hg38"
    CHAIN="${SCRIPT_DIR}/structmummi/external_tools/liftover/hg18ToHg38.over.chain.gz"
    ref_seq="/wibicomfs/STBS/carene/tuto_gpt/tuto2/reference_p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
fi


# ---- First lift the position using liftOver tool
# covert .bim file to ucsc .bed file for liftover
awk '{print "chr"$1, $4-1, $4, $2}' OFS="\t" "$input_".bim > snps_hg18.bed
sort -k1,1 -k2,2n snps_hg18.bed > snps_hg18.sorted.bed

### cmd : liftOver oldFile map.chain newFile unMapped
liftOver snps_hg18.sorted.bed $CHAIN snps_${ref_hg_target}.bed unmapped.bed

# convert back the lifted.bed to a .txt file that will be used to update the plink coordinate
# allvarainst
awk 'BEGIN{OFS="\t"} {print $4, $3}' snps_${ref_hg_target}.bed > snps_${ref_hg_target}.txt

if command -v conda >/dev/null 2>&1; then
    # Initialize conda
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda deactivate
    conda activate gcta_env
else
    echo "ERROR: conda not found in PATH" >&2; exit 1
fi

awk '{print $4}' snps_${ref_hg_target}.bed > lifted.snps
plink --bfile $input_ --extract lifted.snps  --update-map snps_${ref_hg_target}.txt  --make-bed --out "${input_}_${ref_hg_target}"
plink --bfile "${input_}_${ref_hg_target}" --make-bed --out lifted_${output_}_1

echo " ------------------------------"
echo " - First step of lifting, done "
echo " - The variants coordinate lifted to $ref_hg_target "
echo " - Next step, the reference allele need to match the reference build "
echo " ------------------------------"

# --- liftOver tool (from UCSC) only lift the position and do not check the allele plausibility
# convert to vcf, then bcf
$plink2 --bfile lifted_${output_}_1 --recode vcf id-delim=':' --out lifted_${output_}_2 --threads 50

echo "------------"
$BCFTOOLS --version
echo "------------"

$BCFTOOLS view lifted_${output_}_2.vcf -Ob -o lifted_${output_}_1.bcf

# Fix reference allele / strand issues
$BCFTOOLS +fixref lifted_${output_}_1.bcf \
    -Ob \
    -o ${output_}.fixed.bcf \
    -- -f $ref_seq -m flip

# Normalize (indels + multiallelics)
$BCFTOOLS norm \
    -m -any \
    -c x -f $ref_seq \
    -Ob \
    -o ${output_}.norm.bcf \
    ${output_}.fixed.bcf

# Index
$BCFTOOLS index -f ${output_}.norm.bcf

conda deactivate
conda activate gcta_env

# Convert back to vcf then to plink format
plink2 \
  --bcf ${output_}.norm.bcf \
  --double-id \
  --fa $ref_seq \
  --make-bed \
  --out ${output_}


echo "--------------------------------"
echo "----- Your dataset has been successfully lifted ... "
echo "----- $output_ created !!!"
echo "--------------------------------"


## The liftOver tool (from UCSC) only lift based on the position and do not check the allele plausibility
## If I use the tool bcftool +liftover, it check the allele of WTCCC to match with the genome-buil
## But the reference allele of the genome-build and WTCCC does not match as WTCCC is based on SNPs array, and not reference based
## And to match the alleles of WTCCC to the genome-build, I need to align WTCCC to reference panel, BUT reference panel is only based on hg37 or 38.
## convert plink to vcf in order to run liftover with bcftools
## add chr to the chromosome ID in order to match with the chain file
#plink --bfile "${input_}" --output-chr chrM --recode vcf bgz --out "${input_}"
#
#conda deactivate
#conda activate s2_mapping
#
#bcftools index -t ${input_}.vcf.gz
#
#bcftools +liftover $input_.vcf.gz -Ou \
#-- -s $REF_ORIGINAL -f $REF -c $CHAIN --reject rejected_variants.vcf | bcftools sort -Oz -o ${input_}_hg19 -W=tbi
#
#
