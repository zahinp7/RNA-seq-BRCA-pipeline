#!/usr/bin/env bash
# Stage 0: Download raw FASTQ files from SRA
# Dataset:  human breast tumor vs normalGSE183947 
# 3 tumor + 3 normal, smallest files, different donors

set -euo pipefail

OUTDIR="data/raw"
mkdir -p "$OUTDIR" logs

TUMOR=(SRR15852416 SRR15852407 SRR15852399)
NORMAL=(SRR15852426 SRR15852427 SRR15852428)
ALL=("${TUMOR[@]}" "${NORMAL[@]}")

echo "Download started" | tee logs/download.log

for SRR in "${ALL[@]}"; do
    echo "Downloading $SRR..." | tee -a logs/download.log
    prefetch "$SRR" --output-directory "$OUTDIR"
    fastq-dump --outdir "$OUTDIR" \
               --gzip \
               --skip-technical \
               --readids \
               --read-filter pass \
               --dumpbase \
               --split-3 \
               --clip \
               "$OUTDIR/$SRR/$SRR.sra"
    rm -rf "$OUTDIR/$SRR"
    echo "$SRR done" | tee -a logs/download.log
done

echo "All downloads complete" | tee -a logs/download.log
