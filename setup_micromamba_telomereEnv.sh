#!/usr/bin/env bash
set -euo pipefail

module purge
module load Micromamba/2.0.2-0

ENV_NAME="telomereEnv"

# Remove existing environment if it exists
if micromamba env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
    echo "Removing existing environment '${ENV_NAME}'..."
    micromamba env remove -n "${ENV_NAME}" -y
fi

# Clean package caches to avoid using old cached packages
echo "Cleaning micromamba/conda caches..."
micromamba clean --all -y

# Create fresh environment
# Versions are pinned (major.minor, or exact where behavior-sensitive) so re-running this
# script later doesn't silently pull newer package versions and change pipeline behavior.
# Re-verify/update these deliberately rather than relaxing them.
echo "Creating new environment '${ENV_NAME}'..."
micromamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
    python=3.11 \
    r-base=4.4 \
    r-optparse=1.7.5 \
    r-data.table=1.15 \
    r-stringr=1.5 \
    bioconductor-genomicalignments=1.40 \
    pysam=0.22 \
    matplotlib=3.9 \
    pandas=2.2 \
    snakemake=8.25 \
    samtools=1.20

echo "Installing TelomereHunter2..."
# Not version-pinned: exact PyPI release numbers drift too fast to hardcode reliably here.
# The lock file written below records whatever version actually got installed.
micromamba run -n "${ENV_NAME}" pip install telomerehunter2

echo "Environment '${ENV_NAME}' is ready."
echo "Recording exact resolved versions to telomereEnv.lock.yaml for reproducibility..."
micromamba env export -n "${ENV_NAME}" > "$(dirname "${BASH_SOURCE[0]}")/telomereEnv.lock.yaml"
echo "Run Python with:"
echo "  micromamba run -n ${ENV_NAME} python your_script.py"
echo "Run R with:"
echo "  micromamba run -n ${ENV_NAME} Rscript your_script.R"