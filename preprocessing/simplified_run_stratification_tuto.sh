#!/bin/bash

##############################################################
############### START ANALISIS ###############################
##############################################################
# Run stratification on the whole variants, not only on the exons
################################################################
export TMPDIR=/wibicomfs/STBS/carene/tmp

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

plink2="${REPO_DIR}/external_tools/plink/plink2"
script_path="${SCRIPT_DIR}"

ref_1kG=""
ref_1kG_hg19="${REPO_DIR}/references/reference_1000G_hg19/1kG_MDS5"
ref_1kG_hg38="${REPO_DIR}/references/reference_1000G_hg38/allchr"

# The panel for hg19 and hg38 is the same
panel_1000g="${REPO_DIR}/references/reference_1000G_hg19/integrated_call_samples_v3.20130502.ALL.panel"

ref_fasta=""

# Margin to set the delimitation around european samples, based on the reference panel.
# It is the fraction of eur range, with 15% wider than the panel.
# A figure is generated, have a look if it works with you
eur_margin="0.15"

# - For parsing argument - 
# Bash giveds to the script a list of words, where each word is considered
# as argument --> for instance : -i (arg 1) input.txt (arg 2) -o ...
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --ref_hg_target)
      ref_hg_target="$2"
      shift 2
      ;;

    --eur_margin)
      eur_margin="$2"
      shift 2
      ;;

    *)
      echo "Unknown option: $1"
      echo "Usage: /pathtoscript/run_stratification.sh  --ref_hg_target <hg19 or hg38 > [--eur_margin 0.1]"
      exit 1
      ;;
  esac
done

if [[ "$ref_hg_target" == "hg19" ]]; then
    echo "---- The stratification will be performed using 1000G data aligned on $ref_hg_target ---"
    input_1kG="$ref_1kG_hg19"
    ref_fasta="${REPO_DIR}/references/reference_hg19/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz"

elif [[ "$ref_hg_target" == "hg38" ]]; then
    echo "---- The stratification will be performed using 1000G data aligned on $ref_hg_target ---"
    input_1kG="$ref_1kG_hg38"
    ref_fasta="${REPO_DIR}/references/reference_p14/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

fi


# with the lifted variants, the 1000 genomes should has the same refe as my data
# no need to do all the flip steps for 1000G to be merged to my data
plink2 --bfile $input_1kG   --fa $ref_fasta --snps-only just-acgt  --rm-dup force-first --make-bed --out "g1000_strat_1"


# - Keep only shared SNPs between WTCCC and 1000G
comm -12   <(cut -f2 allvariants.bim | sort) <(cut -f2 g1000_strat_1.bim | sort)  > shared_allvar_1000strat.snps
plink2 --bfile "allvariants" --fa $ref_fasta  --extract shared_allvar_1000strat.snps --make-bed --out "allvariants_2"
plink2 --bfile "g1000_strat_1" --fa $ref_fasta  --extract shared_allvar_1000strat.snps --make-bed --out "g1000_strat_2"

# Merge wtccc with 1000G
# After merging, some variants are flipped -> realign them to match the reference build
plink --bfile "allvariants_2" --bmerge g1000_strat_2.bed g1000_strat_2.bim g1000_strat_2.fam  --allow-no-sex --make-bed --out MDS_merge2 # --reference-allele g1000_strat_ref-list.txt
plink2 --bfile "MDS_merge2" --fa $ref_fasta --ref-from-fa --make-bed --out "MDS_merge3"

#https://cloufield.github.io/GWASTutorial/05_PCA/#sample-codes
# - generate the .prune.in file on the whole variants
plink2 --bfile "MDS_merge3" --fa $ref_fasta --exclude range  "$script_path/inversion.txt" --indep-pairwise 50 5 0.2 --out indepSNP_allvar

# Remove related samples using king-cutoff
plink2  --bfile "MDS_merge3" --fa $ref_fasta --extract indepSNP_allvar.prune.in --king-cutoff 0.0884 --threads 50 --out kingcuttof

# Projection (related and unrelated samples)
plink2  --bfile "MDS_merge3" --fa $ref_fasta --keep kingcuttof.king.cutoff.in.id  --extract indepSNP_allvar.prune.in --freq counts --threads 50 --pca allele-wts 10 --out "MDS_merge3"

plink2  --bfile "MDS_merge3" --threads 50 --read-freq MDS_merge3.acount --score MDS_merge3.eigenvec.allele 2 6 header-read no-mean-imputation variance-standardize --score-col-nums 7-16 --out MDS_merge3_projected  # 7-16


# Convert population codes into superpopulation codes (i.e., AFR,AMR,ASN, and EUR).
awk '{print$1,$1,$2}' $panel_1000g > race_1kG.txt
sed 's/JPT/ASN/g' race_1kG.txt>race_1kG2.txt
sed 's/ASW/AFR/g' race_1kG2.txt>race_1kG3.txt
sed 's/CEU/EUR/g' race_1kG3.txt>race_1kG4.txt
sed 's/CHB/ASN/g' race_1kG4.txt>race_1kG5.txt
sed 's/CHD/ASN/g' race_1kG5.txt>race_1kG6.txt
sed 's/YRI/AFR/g' race_1kG6.txt>race_1kG7.txt
sed 's/LWK/AFR/g' race_1kG7.txt>race_1kG8.txt
sed 's/TSI/EUR/g' race_1kG8.txt>race_1kG9.txt
sed 's/MXL/AMR/g' race_1kG9.txt>race_1kG10.txt
sed 's/GBR/EUR/g' race_1kG10.txt>race_1kG11.txt
sed 's/FIN/EUR/g' race_1kG11.txt>race_1kG12.txt
sed 's/CHS/ASN/g' race_1kG12.txt>race_1kG13.txt
sed 's/PUR/AMR/g' race_1kG13.txt>race_1kG14.txt

# Create a racefile of your own data.
awk '{print$1,$2,"OWN"}' allvariants_2.fam>racefile_own.txt

# Concatenate racefiles.
cat race_1kG14.txt racefile_own.txt | sed -e '1i\FID IID race' > racefile.txt

set +u
conda deactivate
conda activate jupyter_env2
set -u

python3 $script_path/plot_PCA.py --margin "$eur_margin"

set +u
conda deactivate
conda activate gcta_env
set -u 

## Keep samples within european group
## Samples belonging to outliers or other ancestries are filtered out
## The cut-off are based on the reference panels (1kG) on the first two dimensional of PCA.
## User are responsible for modifying the cut-off if necessary 

# Read the EUR boundaries computed by plot_PCA.py
read -r xmin xmax ymin ymax < <(tail -n 1 EUR_boundaries.txt)
echo "---- Keeping samples inside PC1 [$xmin, $xmax] and PC2 [$ymin, $ymax]"

# Keep the samples inside the EUR box.
awk -v xmin="$xmin" -v xmax="$xmax" -v ymin="$ymin" -v ymax="$ymax" '
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "PC1_AVG") c1 = i
        if ($i == "PC2_AVG") c2 = i
    }
    if (!c1 || !c2) { print "ERROR: PC1_AVG/PC2_AVG not found in header" > "/dev/stderr"; exit 1 }
    next
}
$c1 >= xmin && $c1 <= xmax && $c2 >= ymin && $c2 <= ymax { print $1, $2 }
' MDS_merge3_projected.sscore > EUR_samples

##awk 'NR>1 && ($6 > 0  && $7 < 0.036) {print $1, $2}' MDS_merge3_projected.sscore > EUR_samples # for UC
##the previous command I used before automated delimitation

### Extract these individuals in HapMap data.
plink2 --bfile "allvariants" --fa $ref_fasta --keep EUR_samples --make-bed --out "plink_13"
## Note, since our HapMap data did include any ethnic outliers, no individuls were removed at this step. However, if our data would have included individuals outside of the thresholds we set, then these individuals would have been removed.
## MDS is applied to Identify population structure
## Create covariates based on MDS.
## Perform an MDS ONLY on HapMap data without ethnic outliers. The values of the 10 MDS dimensions are subsequently used as covariates in the association analysis in the third tutorial.
plink2 --bfile "plink_13" --fa $ref_fasta --extract indepSNP_allvar.prune.in --pca 10 --out "plink_13"
#
# Change the format of the .mds file into a plink covariate file.
awk '{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' plink_13.eigenvec > covar_pca.txt
## The values in covar_mds.txt will be used as covariates, to adjust for remaining population stratification, in the third tutorial where we will perform a genome-wide association analysis.
