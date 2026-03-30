#!/usr/bin/env bash
# Stage 1: Raw QC on all FASTQ files

set -euo pipefail

INDIR="data/raw"
OUTDIR="results/qc/raw_fastqc"
mkdir -p "$OUTDIR"

echo "Running FastQC on raw reads..."

fastqc "$INDIR"/*.fastq.gz \
    --outdir "$OUTDIR" \
    --threads 4

echo "Running MultiQC..."
multiqc "$OUTDIR" \
    --outdir results/qc \
    --filename multiqc_raw \
    --force

echo "Raw QC complete"
