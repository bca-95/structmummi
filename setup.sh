#!/bin/bash
# =====================================================================
# FOR USERS: downloads tools, packages, environments, and references 
# Usage:  bash setup.sh
# Needs:  Linux x86_64, conda, curl, ~<XX> GB of free disk space
# =====================================================================
set -euo pipefail

REPO_DIR="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"

# --- Where the files are stored (Zenodo record)
ZENODO_RECORD="10.5072/zenodo.633483"      
BASE_URL="${BASE_URL:-https://zenodo.org/records/${ZENODO_RECORD}/files}"

FILES="external_tools.tar.gz  gcta_env.tar.gz liftover_env.tar.gz reference_1000G_hg19.tar.gz reference_1000G_hg38.tar.gz"
DL_DIR="${REPO_DIR}/downloads"

# --- Checks
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found in PATH" >&2; exit 1; }
command -v curl  >/dev/null 2>&1 || { echo "ERROR: curl not found" >&2; exit 1; }

# Conda environments are installed where 'conda activate <name>' finds them.
# If you can't write there (shared conda on a cluster), run:
#   CONDA_ENVS_DIR=~/.conda/envs bash setup.sh
ENVS_DIR="${CONDA_ENVS_DIR:-$(conda info --base)/envs}"
mkdir -p "$ENVS_DIR"
[[ -w "$ENVS_DIR" ]] || { echo "ERROR: cannot write to $ENVS_DIR. Run: CONDA_ENVS_DIR=~/.conda/envs bash setup.sh" >&2; exit 1; }

# --- 1. Download (resumes if interrupted)
mkdir -p "$DL_DIR"
cd "$DL_DIR"
for f in $FILES; do
    echo "--- Downloading $f"
    curl -L -C - --fail -o "$f" "${BASE_URL}/${f}?download=1"
done

# --- 2. Check the files are complete and not corrupted
echo "--- Checking files"
md5sum -c "${REPO_DIR}/checksums.md5"

# --- 3. Unpack tools and references into the repo
echo "--- Installing external tools and references"
for f in external_tools.tar.gz reference_1000G_hg19.tar.gz reference_1000G_hg38.tar.gz; do
    tar -xzf "$f" -C "$REPO_DIR"
done

# --- 4. Unpack the conda environments
for env in gcta_env liftover_env; do
    target="${ENVS_DIR}/${env}"
    if [[ -d "$target" ]]; then
        echo "--- Environment $env already exists in $target, skipped"
        continue
    fi
    echo "--- Installing conda environment $env"
    mkdir -p "$target"
    tar -xzf "${env}.tar.gz" -C "$target"
    "${target}/bin/conda-unpack"     # fixes the paths inside the environment
done

# --- 5. Clean up
cd "$REPO_DIR"
rm -rf "$DL_DIR"

echo
echo "Setup complete. Test with:"
echo "  conda activate gcta_env && plink --version"
echo "  ${REPO_DIR}/external_tools/plink/plink2 --version"

