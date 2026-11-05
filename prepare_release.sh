#!/bin/bash
# =====================================================================
# FOR THE MAINTAINER ONLY (you), not for users.
# Packs the external tools, the references and the conda environments
# into .tar.gz files to upload to Zenodo, and writes their checksums.
#
# Usage (from anywhere):  bash prepare_release.sh
# Output:                 structmummi/release_files/*.tar.gz  -> upload these to Zenodo
#                         structmummi/checksums.md5           -> commit this to GitHub
# =====================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
OUT="${REPO_DIR}/release_files"
mkdir -p "$OUT"
cd "$REPO_DIR"

# --- conda-pack is needed to pack the environments
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! command -v conda-pack >/dev/null 2>&1; then #with every package at exactly my version.
    echo "ERROR: conda-pack not found. Install it once with:" >&2
    echo "       conda install -n base -c conda-forge conda-pack" >&2
    exit 1
fi

echo "--- 1/4 Packing external tools"
tar -czf "$OUT/external_tools.tar.gz" external_tools/plink external_tools/bcftools-1.23.1

#LATER
echo "--- 2/4 Packing references"
tar -czf "$OUT/reference_1000G_hg19.tar.gz" references/reference_1000G_hg19
tar -czf "$OUT/reference_1000G_hg38.tar.gz" references/reference_1000G_hg38

echo "--- 3/4 Packing conda environments"
for env in gcta_env liftover_env; do
    conda-pack -n "$env" -o "$OUT/${env}.tar.gz" --ignore-missing-files --force
    # Also keep a readable list of the exact versions in the repo
    conda env export -n "$env" --no-builds > "${REPO_DIR}/${env}.yml"
done

echo "--- 4/4 Writing checksums"
(cd "$OUT" && md5sum *.tar.gz) > "${REPO_DIR}/checksums.md5"

echo
echo "Done. Files to upload to Zenodo:"
ls -lh "$OUT"
echo
echo "Then commit checksums.md5 and the .yml files to GitHub."
