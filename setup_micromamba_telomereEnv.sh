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
echo "Creating new environment '${ENV_NAME}'..."
micromamba create -y -n "${ENV_NAME}" -c conda-forge -c bioconda \
    python=3.11 \
    r-base=4.4 \
    r-optparse \
    r-data.table \
    r-stringr \
    bioconductor-genomicalignments \
    bioconductor-bsgenome.hsapiens.ucsc.hg19 \
    pysam \
    matplotlib \
    pandas \
    snakemake \
    samtools

echo "Environment '${ENV_NAME}' is ready."
echo "Run Python with:"
echo "  micromamba run -n ${ENV_NAME} python your_script.py"
echo "Run R with:"
echo "  micromamba run -n ${ENV_NAME} Rscript your_script.R"