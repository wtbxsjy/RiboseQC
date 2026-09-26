#!/usr/bin/env bash
# Download the nf-core/riboseq test dataset (human chr20) used for benchmarking.
# Files come from the nf-core/test-datasets "modules" branch.
set -euo pipefail
DEST=${1:-"$(dirname "$0")/data"}
BASE=https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/riboseq_expression
mkdir -p "$DEST"
cd "$DEST"
for f in Homo_sapiens.GRCh38.111_chr20.gtf Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz \
         aligned_reads/SRX11780887_chr20.bam aligned_reads/SRX11780887_chr20.bam.bai \
         aligned_reads/SRX11780888_chr20.bam aligned_reads/SRX11780888_chr20.bam.bai; do
    [ -s "$(basename "$f")" ] || curl -sSfL -o "$(basename "$f")" "$BASE/$f"
done
[ -s Homo_sapiens.GRCh38.dna.chromosome.20.fa ] || gunzip -k Homo_sapiens.GRCh38.dna.chromosome.20.fa.gz
ls -la
