#!/usr/bin/env bash
# Stage 3a: Download human reference genome and GTF annotation

set -euo pipefail

REFDIR="reference"
mkdir -p "$REFDIR"

echo "Downloading GRCh38 genome (primary assembly)..."
wget -P "$REFDIR" https://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

echo "Downloading GTF annotation..."
wget -P "$REFDIR" https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz

echo "Decompressing..."
gunzip "$REFDIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gunzip "$REFDIR/Homo_sapiens.GRCh38.109.gtf.gz"

echo "Reference download complete"
