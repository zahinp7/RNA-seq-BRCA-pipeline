#!/usr/bin/env bash
# Stage 5: Count reads per gene using featureCounts
# Run locally: bash scripts/05_featurecounts.sh
# Run on HPC (Bridges-2): submit as SLURM job with --cpus-per-task=8 --mem-per-cpu=2000M

set -euo pipefail

ALIGNDIR="data/aligned"
COUNTDIR="data/counts"
GTF="reference/Homo_sapiens.GRCh38.109.gtf"
THREADS=${1:-4}
mkdir -p "$COUNTDIR"

BAM_FILES=$(ls "$ALIGNDIR"/*_Aligned.sortedByCoord.out.bam | tr '\n' ' ')

echo "Running featureCounts..."
featureCounts \
    -T "$THREADS" \
    -p \
    --countReadPairs \
    -s 2 \
    -a "$GTF" \
    -o "$COUNTDIR/counts.txt" \
    $BAM_FILES

echo "featureCounts complete"