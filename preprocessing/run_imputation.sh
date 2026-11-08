#!/bin/bash
set -euo pipefail

export TMPDIR=/scratch/SCRATCH_SAS/carene/tmp
export TEMP=/scratch/SCRATCH_SAS/carene/tmp
export TMP=/scratch/SCRATCH_SAS/carene/tmp

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
conda activate gcta_env
set -u

# Use compiled bcftools, the conda version raise error for the lifted data
# hwv, the conda version does not bring any error on 1000G data
export BCFTOOLS="${REPO_DIR}/external_tools/bcftools-1.23.1/bcftools"
export BCFTOOLS_PLUGINS="${REPO_DIR}/external_tools/bcftools-1.23.1/plugins"

export run_minimac="${REPO_DIR}/external_tools/Minimac4/build/minimac4"

run_eagle="${REPO_DIR}/external_tools/Eagle_v2.4.1/eagle"
run_shapeit="${REPO_DIR}/external_tools/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin/shapeit"

# Tutorial
#https://cloufield.github.io/GWASTutorial/Phasing/#phasing-using-eagle

input_="plink_13" # This only contain EUR sample
geneticmap=""
ref_fasta=""
vcf_ref_panel=""
# Reference panel same for hg19 and hg38
ref_panel="${REPO_DIR}/reference_1000G_hg19/integrated_call_samples_v3.20130502.ALL.panel"


ref_hg_target=""
cpu_threads_phasing=40
cpu_threads_imputation=14

# - For parsing argument - 
# Bash giveds to the script a list of words, where each word is considered
# as argument --> for instance : -i (arg 1) input.txt (arg 2) -o ...
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --ref_hg_target)

      if [[ -z "$2" ]]; then
        echo "Error: --ref_hg_target requires a value (hg19 or hg38)"
        exit 1
      fi

      if [[ "$2" != "hg19" && "$2" != "hg38" ]]; then
        echo "Error: --ref_hg_target must be 'hg19' or 'hg38'"
        exit 1
      fi

      ref_hg_target="$2"
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

    *)
      echo "Unknown option: $1"
      echo "Usage: /pathtoscript/run_imputation.sh --ref_hg_target <hg19 or hg38> "
      exit 1
      ;;
  esac
done

if [[ "$ref_hg_target" == "hg19" ]]; then
    echo " ---- Imputation using reference panel built on hg19 ----"
    ref_fasta="${REPO_DIR}/references/reference_hg19/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"
    geneticmap="${REPO_DIR}/external_tools/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz"
    vcf_ref_panel="${REPO_DIR}/references/reference_1000G_hg19/prepare_impute/"
elif [[ "$ref_hg_target" == "hg38" ]]; then
    echo " ---- Imputation using reference panl built on hg38"
    ref_fasta="${REPO_DIR}/references/reference_p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    geneticmap="${REPO_DIR}/external_tools/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz"
    vcf_ref_panel="${REPO_DIR}/references/reference_1000G_hg38/prepare_impute/"
fi

# Make sure that allele match with reference base on genome assembly GRch37/GRCh38
plink2 \
    --bfile ${input_} \
    --fa ${ref_fasta} \
    --ref-from-fa \
    --make-bed \
    --out plink_imp1

# split data by chromosome 
for chr in {1..22}; do
#    plink2 --bfile ${input_} \
    plink2 --bfile plink_imp1 \
          --fa ${ref_fasta} \
          --chr $chr \
          --make-bed \
          --out plink_imp1_${chr}
done

# convert data to vcf then to bcf (bcf format allow better time performance)
# Convert ecah chromosome to vcf
for chr in {1..22}; do
    plink2 --bfile plink_imp1_${chr} \
          --fa ${ref_fasta} \
          --recode vcf bgz id-delim=':' \
          --out plink_imp1_${chr} \
          --threads 40
done

conda deactivate
conda activate s2_mapping

# Perform the phasing using eagle
# For optimal phasing, in the reference panel only contain EUR population

# I need to normalize as well the wtccc, same as what have been done in the reference panel
for chr in {1..22}; do
  $BCFTOOLS index -f  "plink_imp1_${chr}.vcf.gz"

  $BCFTOOLS view --no-version -Ou -c 2  plink_imp1_${chr}.vcf.gz | $BCFTOOLS norm --no-version -Ou -m -any | $BCFTOOLS norm --no-version -Ob -o plink_imp1_${chr}.bcf -d none -c x -f $ref_fasta && \
  $BCFTOOLS index -f plink_imp1_${chr}.bcf   
    echo "chr $chr done ..." 
done

process_eagle(){
  chr="$1"

  echo "--- Phasing of chr$chr started ... "
  $run_eagle \
    --vcfTarget="plink_imp1_${chr}.bcf" \
    --vcfRef="${vcf_ref_panel}/eur_ref_chr${chr}.bcf" \
    --geneticMapFile=${geneticmap} \
    --outPrefix="plink_imp2_${chr}" \
    --numThreads=$cpu_threads_phasing \
    --chrom=${chr}
 
    $BCFTOOLS index -f plink_imp2_${chr}.vcf.gz    
  
    echo "chr$chr done ..." 
}

#Use shapeit2 to convert .hap and .sample to VCF
# I don't need shapeit2 as the previous output is already a vcf

# ---- Run imputation itself
process_minimac() {
    chr="$1"

    echo "--- chr $chr started..."

    $run_minimac \
        "${vcf_ref_panel}/eur_ref_chr${chr}.msav" "plink_imp2_${chr}.vcf.gz" \
        -o "plink_imp3_${chr}.vcf.gz" \
        --format GT,DS,HDS,GP,SD \
        --cpu $cpu_threads_imputation

    $BCFTOOLS index -f "plink_imp4_${chr}.vcf.gz"

    echo "--- chr $chr done..."
}

export -f process_minimac
export -f process_eagle
export run_eagle
export geneticmap
export run_minimac
export vcf_ref_panel
export BCFTOOLS

# run 3 chromosomes in parallel
parallel -j 4 process_eagle ::: {1..22}
parallel -j 4 process_minimac ::: {1..22}

# - Parameter estimation




