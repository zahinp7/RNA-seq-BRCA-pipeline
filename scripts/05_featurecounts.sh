#!/usr/bin/env bash
# Stage 5: Count reads per gene using featureCounts

set -euo pipefail

ALIGNDIR="data/aligned"
COUNTDIR="data/counts"
GTF="reference/Homo_sapiens.GRCh38.109.gtf"
mkdir -p "$COUNTDIR"

BAM_FILES=$(ls "$ALIGNDIR"/*_Aligned.sortedByCoord.out.bam | tr '\n' ' ')

echo "Running featureCounts..."
featureCounts \
    -T 16 \
    -p \
    --countReadPairs \
    -s 2 \
    -a "$GTF" \
    -o "$COUNTDIR/counts.txt" \
    $BAM_FILES

echo "featureCounts complete"