#!/usr/bin/env bash
# Stage 4: Align trimmed reads to GRCh38 using STAR
# Run locally: bash scripts/04_align.sh
# Run on HPC (Bridges-2): submit as SLURM job with --cpus-per-task=16 --mem-per-cpu=2000M

set -euo pipefail

TRIMDIR="data/trimmed"
ALIGNDIR="data/aligned"
INDEXDIR="reference/star_index"
THREADS=${1:-4}  # default 4 threads locally, override with: bash 04_align.sh 16
mkdir -p "$ALIGNDIR"

SAMPLES=(
    SRR15852399
    SRR15852407
    SRR15852416
    SRR15852426
    SRR15852432
    SRR15852438
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Aligning $SAMPLE..."
    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$INDEXDIR" \
        --readFilesIn "$TRIMDIR/${SAMPLE}_R1_trimmed.fastq.gz" \
                      "$TRIMDIR/${SAMPLE}_R2_trimmed.fastq.gz" \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS NM MD \
        --outFileNamePrefix "$ALIGNDIR/${SAMPLE}_" \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000
    echo "$SAMPLE done"
done

echo "Alignment complete"
