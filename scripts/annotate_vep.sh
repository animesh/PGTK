#!/bin/bash
# Annotate merged SNV VCFs with VEP using local GTF + genome
# PASS variants only; extract gene consequences for downstream comparison

set -euo pipefail

INDIR="/mnt/z/Download/TK/results/mutect2_merged"
OUTDIR="/mnt/z/Download/TK/results/vep"
THREADS=4

mkdir -p "$OUTDIR"

annotate_one() {
    local vcf="$1"
    local group
    group=$(basename "$vcf" .merged.snv.vcf.gz)
    local out="$OUTDIR/${group}.vep.tsv"

    [[ -s "$out" ]] && { echo "[SKIP] $group already done"; return 0; }

    echo "[$(date '+%H:%M:%S')] Annotating $group ..."

    bcftools view -f PASS "$vcf" | \
    docker run --rm -i \
        -v /mnt/z/Download/TK:/data \
        ensemblorg/ensembl-vep:latest \
        vep \
            --input_file STDIN \
            --output_file STDOUT \
            --format vcf \
            --tab \
            --no_stats \
            --gtf /data/Homo_sapiens.GRCh38.110.gtf.gz \
            --fasta /data/genome.fa \
            --fields "Gene,SYMBOL,Consequence,IMPACT,Feature,HGVSc,HGVSp" \
            --hgvs \
            --numbers \
            --canonical \
            --symbol \
            --no_intergenic \
            --buffer_size 5000 \
            --fork "$THREADS" \
        2>/dev/null \
        > "$out"

    echo "[$(date '+%H:%M:%S')] Done $group"
}

export -f annotate_one
export THREADS

ls "$INDIR"/*.merged.snv.vcf.gz | parallel -j 3 annotate_one {}

echo "All annotation done."
