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
conda activate gcta_env
set -u

# In this script, were are going to update the allele in WTCCC to a reference such as 1000 genomes. 
# --> This is called allele harmonization. 
# Allele harmonization is the process of aligning allele labels, strands, and reference conventions across datasets so that the same genetic variant is represented consistently.
# The allele harmonization will allow to harmonize the allele (especially) the REF allele across different and extrenal datasets than WTCCC.
# 1. For instance, eeventhough WTCCC is aligned based on genome-build hg18 (or grch36), the reference allele does not match to the genome-build. Because WTCCC dataset originates from SNP arrays and not sequence-based. It uses assay-defined alleles, not reference genome alleles
# 2 . This could generate mismatch when wanting to use external data such as when performing IMPUTATION
# 3 . Also, if we train the data on WTCCC specific allele, it will learn dataset-specific allele. By harmonizing, we can generaliwe the testing across different datasets 

# -- /!\ --- Prior allele harmonization, wtccc position has to be lifted to hg19 in order to match with the 1000 genomes  

# Initialize variables
# The input should be the basename of plink file, eg. "data".bim
input_data=""
output_name="output_name"
g1000_processed="g1000_data_3.bim"

# - For parsing argument - 
# Bash giveds to the script a list of words, where each word is considered
# as argument --> for instance : -i (arg 1) input.txt (arg 2) -o ...
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --input_data)
      input_data="$2"
      shift 2
      ;;
    --output_name)
      output_name="$2"
      shift 2
      ;;
    *)
      echo "Unknown option: $1"
      echo "Usage: /pathtoscript/qc_filtering.sh --input_control plink_basename --input_case plink_basename --convert_to_rsid <yes or no (default: yes) --liftover <yes or no (default: yes)"
      exit 1
      ;;
  esac
done


# --- We conduct the harmonization in 3 steps : 
# 1) Make sure the reference genome is similar in the WTCCC and the 1000 Genomes Project datasets.
# 2) Resolve strand issues.
# 3) Remove the SNPs which after the previous two steps still differ between datasets.

# The following steps are maybe quite technical in terms of commands, but we just compare the two data sets and make sure they correspond.

# 1) set reference genome
# Here the reference rsid is set to chr:pos to match with WTCCC data

awk '{print$2,$5}' $g1000_processed > g1000_ref-list.txt
plink --bfile $input_data --reference-allele g1000_ref-list.txt --make-bed --out plink_data_1
#plink --bfile $input_data --reference-allele "1kG_${dataset}_SNPs.txt" --make-bed --out plink_data_1

# The case and the control have the same reference genome for all SNPs.
# This command will generate some warnings for impossible A1 allele assignment.

# 2) Resolve strand issues.
# Check for potential strand issues.
awk '{print$2,$5,$6}' $g1000_processed > plink_1000g_1_tmp
awk '{print$2,$5,$6}' plink_data_1.bim > plink_data_1_tmp
sort plink_1000g_1_tmp plink_data_1_tmp |uniq -u > all_differences.txt
# 1624 differences between the files, some of these might be due to strand issues.

## Flip SNPs for resolving strand issues.
# Print SNP-identifier and remove duplicates.
awk '{print$1}' all_differences.txt | sort -u > flip_list.txt

# Generates a file of 812 SNPs. These are the non-corresponding SNPs between the two files.
# Flip the 812 non-corresponding SNPs.
plink --bfile plink_data_1 --flip flip_list.txt --reference-allele g1000_ref-list.txt --make-bed --out plink_data_2

# Check for SNPs which are still problematic after they have been flipped.
awk '{print$2,$5,$6}' plink_data_2.bim > plink_data_2_tmp
sort plink_1000g_1_tmp plink_data_2_tmp |uniq -u  > uncorresponding_SNPs.txt
# This file demonstrates that there are 84 differences between the files.

# 3) Remove problematic SNPs from HapMap and 1000 Genomes.
awk '{print$1}' uncorresponding_SNPs.txt | sort -u > SNPs_for_exlusion.txt
# The command above generates a list of the 42 SNPs which caused the 84 differences between the HapMap and the 1000 Genomes data sets after flipping and setting of the reference genome.

# Remove the 42 problematic SNPs from both datasets.
plink --bfile plink_data_2 --exclude SNPs_for_exlusion.txt --make-bed --out $output_name



