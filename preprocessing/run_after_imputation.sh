#!/bin/bash
#https://cloufield.github.io/GWASTutorial/Imputation/#quality-control
# The previous script run_imputation.sh perform phasing and imputation based on the
# reference panel 1000G. This outputedi .VCF files including variants from the input
# datasets (to be imputed) and the reference panel.
# This current script is to perform after imputation. It will retrieve variants of the input datasets,
# perform simple quality control filtering on the imputed variants and convert the files to PLINK format. 

set -euo pipefail

export TMPDIR="${TMPDIR:-/tmp}"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
mkdir -p "$TMPDIR"

# --- Locate the repo, wherever it was cloned
SCRIPT_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)" # ../structmummi/preprocessing
REPO_DIR="$(dirname "$SCRIPT_DIR")" # dirname remove the last part of the path # ../structmummi

# --- Find conda automatically
if ! command -v conda >/dev/null 2>&1; then
    echo "ERROR: conda not found in PATH" >&2; exit 1
fi
source "$(conda info --base)/etc/profile.d/conda.sh"

set +u
conda deactivate
conda activate gcta_env
set -u 


# Use compiled bcftools, the conda version raise error for the lifted data
# hwv, the conda version does not bring any error on 1000G data
export BCFTOOLS="${REPO_DIR}/external_tools/bcftools-1.23.1/bcftools"
export BCFTOOLS_PLUGINS="${REPO_DIR}/external_tools/bcftools-1.23.1/plugins"

ref_fasta="${REPO_DIR}/references/reference_p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
script_path="${SCRIPT_DIR}"

# Initialize variable
hc_threshold=0.1

# - For parsing argument - 
# Bash giveds to the script a list of words, where each word is considered
# as argument --> for instance : -i (arg 1) input.txt (arg 2) -o ...
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --hc_threshold)
      hc_threshold="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: /pathtoscript/run_after_imputation.sh --hc_threshold <float (default: 0.1)>"
      exit 1
      ;;
  esac
done



## I stored all variants (input+reference panel)
## I need to only focus on input variants
## MODIFY THE SCRIPT LATER
cut -f2 plink_13.bim > original.snps

seq 1 22 | parallel -j 5 '
    $BCFTOOLS view \
        -i "ID=@original.snps" \
        -Oz -o "plink_imp4_{}.vcf.gz" "plink_imp3_{}.vcf.gz"
'

# Analyze the R2 score of the imputed variants
rm r2_maf_imputed.txt

# header
echo -e "STATUS\tR2\tMAF" > r2_maf_imputed.txt

for chr in {1..22}; do
  $BCFTOOLS query \
    -f '%INFO/TYPED\t%INFO/IMPUTED\t%INFO/R2\t%INFO/MAF\n' \
    plink_imp4_${chr}.vcf.gz \
  | awk '{
      if ($1=="1" && $2=="1") status="TYPED+IMPUTED";
      else if ($1=="1") status="TYPED";
      else if ($2=="1") status="IMPUTED";
      else status="UNKNOWN";
      print status "\t" $3 "\t" $4
    }' \
  >> r2_maf_imputed.txt
  echo "--- chr $chr done..."

done

set +u
conda deactivate
conda activate jupyter_env2
set -u 

python3 $script_path/plot_R2_imputation.py

set +u
conda deactivate
conda activate gcta_env
set -u 


# The Rsq metric efficiently discriminates well-imputed variants at a threshold of 0.7 (statement in the tutorial)

echo "--- filter out bad imputed variant, with R2 < 0.7"
for chr in {1..22}; do
    $BCFTOOLS filter -i 'INFO/R2>=0.7' "plink_imp4_${chr}.vcf.gz" -Oz -o "plink_imp5_${chr}.vcf.gz"
    $BCFTOOLS index -f "plink_imp5_${chr}.vcf.gz"
    echo "--- chr $chr done..."
done

#### convert to plink using best-guess genotypes (0/1/2) 
#### But it is suggested to use "dosage" but it is continous value
#### I tried to to dosage BUT I have more missing genotype than my original data
###
###echo "--- convert the imputed .vcf back to plink"
###for chr in {1..22}; do
###    plink2 \
###    --vcf "plink_imp5_${chr}.vcf.gz" \
###    --make-pgen \
###    --out plink_imp5_${chr}.dosage
###    echo "--- chr $chr done..."
###done
##
#### convert to plink using dosage 
#### But it is suggested to use "dosage" but it is continous value
###echo "--- convert the imputed .vcf back to plink"
###for chr in {1..22}; do
###    plink2 --vcf plink_imp5_${chr}.vcf.gz dosage=DS  --make-pgen --out plink_imp5_${chr}.dosage #--double-id --import-dosage-certainty 0.8
###done

# Use genotype probability to generate hard-calls
for chr in {1..22}; do
    plink2 --vcf plink_imp5_${chr}.vcf.gz \
    --hard-call-threshold $hc_threshold \
    --make-pgen --out plink_imp5_${chr}.dosage 
done


# Merge all the chromosome into one file
# The .fam ID has been modified from going to .vcf to plink --> rest
rm merge_imputed_chr.txt

for chr in {2..22}; do
  echo "plink_imp5_${chr}.dosage" >> merge_imputed_chr.txt
done

plink2 \
  --pfile plink_imp5_1.dosage \
  --pmerge-list merge_imputed_chr.txt \
  --make-pgen \
  --out plink_imp5_all

# convert to bed
plink2 --pfile "plink_imp5_all" --fa $ref_fasta  --make-bed --out "plink_imp6_all"

# update fam
cp plink_13.fam plink_imp6_all.fam

plink2 --bfile "plink_imp6_all" --fa $ref_fasta --make-bed --out "plink_imp7_all"

# This data contains variants located in the exons + intron
# only select the exons variants
# and perform the full quality control 
