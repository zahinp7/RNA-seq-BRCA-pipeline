#!/usr/bin/env bash
# Stage 3b: Build STAR genome index

set -euo pipefail

REFDIR="reference"
INDEXDIR="reference/star_index"
mkdir -p "$INDEXDIR"

echo "Building STAR index..."
STAR \
    --runMode genomeGenerate \
    --genomeDir "$INDEXDIR" \
    --genomeFastaFiles "$REFDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa" \
    --sjdbGTFfile "$REFDIR/Homo_sapiens.GRCh38.109.gtf" \
    --sjdbOverhang 149 \
    --runThreadN 4

echo "STAR index complete"
