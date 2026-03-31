#!/usr/bin/env bash
# Stage 4: Align trimmed reads to GRCh38 using STAR

set -euo pipefail

TRIMDIR="data/trimmed"
ALIGNDIR="data/aligned"
INDEXDIR="reference/star_index"
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
        --runThreadN 16 \
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
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000
    echo "$SAMPLE aligned"
done

echo "Alignment complete"