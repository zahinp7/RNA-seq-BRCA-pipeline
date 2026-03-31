#!/usr/bin/env bash
# Stage 2: Adapter trimming and quality filtering with fastp

set -euo pipefail

INDIR="data/raw"
OUTDIR="data/trimmed"
QCDIR="results/qc/fastp"
mkdir -p "$OUTDIR" "$QCDIR"

SAMPLES=(
    SRR15852399
    SRR15852407
    SRR15852416
    SRR15852426
    SRR15852432
    SRR15852438
)

for SAMPLE in "${SAMPLES[@]}"; do
    echo "Trimming $SAMPLE..."
    fastp \
        --in1  "$INDIR/${SAMPLE}_pass_1.fastq.gz" \
        --in2  "$INDIR/${SAMPLE}_pass_2.fastq.gz" \
        --out1 "$OUTDIR/${SAMPLE}_R1_trimmed.fastq.gz" \
        --out2 "$OUTDIR/${SAMPLE}_R2_trimmed.fastq.gz" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 36 \
        --thread 4 \
        --json "$QCDIR/${SAMPLE}_fastp.json" \
        --html "$QCDIR/${SAMPLE}_fastp.html"
    echo "$SAMPLE done"
done

echo "Running MultiQC on fastp reports..."
multiqc "$QCDIR" \
    --outdir results/qc \
    --filename multiqc_trimmed \
    --force

echo "Trimming complete"
