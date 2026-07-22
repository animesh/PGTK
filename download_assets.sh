#!/bin/bash

set -euo pipefail

WORKDIR="${1:-$PWD}"
CONTAINER_DIR="$WORKDIR/singularity_cache"
REFERENCE_DIR="$WORKDIR/reference_downloads"
DOWNLOAD_TMP_DIR="$WORKDIR/download_tmp"

export TMPDIR="$DOWNLOAD_TMP_DIR"

#export http_proxy= "http://proxy.saga:3128/"
#export https_proxy= "http://proxy.saga:3128/"
#export HTTP_PROXY="$http_proxy"
#export HTTPS_PROXY="$https_proxy"
#export no_proxy="localhost,127.0.0.1"
#export NO_PROXY="$no_proxy"

mkdir -p "$REFERENCE_DIR" "$DOWNLOAD_TMP_DIR"

command -v curl >/dev/null
command -v gzip >/dev/null
command -v tar >/dev/null
command -v singularity >/dev/null

resume_download() {
    local url="$1"
    local output_path="$2"
    local description="$3"
    local attempt=0
    local max_attempts=500
    local current_size=0

    while (( attempt < max_attempts )); do
        attempt=$((attempt + 1))
        current_size=$(stat -c '%s' "$output_path" 2>/dev/null || echo 0)

        echo
        echo "$description"
        echo "Attempt ${attempt}/${max_attempts}; existing bytes: ${current_size}"

        if curl \
            --fail \
            --location \
            --continue-at - \
            --connect-timeout 60 \
            --speed-time 120 \
            --speed-limit 1024 \
            --output "$output_path" \
            "$url"; then
            echo "Completed: $output_path"
            return 0
        fi

        current_size=$(stat -c '%s' "$output_path" 2>/dev/null || echo 0)
        echo "Transfer interrupted; retained bytes: ${current_size}"
        sleep 20
    done

    echo "ERROR: failed to complete $url after ${max_attempts} resumable attempts" >&2
    return 1
}

validate_gzip() {
    local file="$1"
    echo "Validating gzip archive: $file"
    gzip -t "$file"
}

validate_tar() {
    local file="$1"
    echo "Validating tar.gz archive: $file"
    tar -tzf "$file" >/dev/null
}

validate_vep_structure() {
    local archive="$1"
    local listing="$DOWNLOAD_TMP_DIR/vep_contents.txt"

    echo "Validating VEP archive structure"
    tar -tzf "$archive" > "$listing"

    if ! grep -m 1 '^homo_sapiens/111_GRCh38/' "$listing" >/dev/null; then
        echo "ERROR: expected homo_sapiens/111_GRCh38 directory is absent from $archive" >&2
        return 1
    fi
}

validate_arriba_structure() {
    local archive="$1"
    local listing="$DOWNLOAD_TMP_DIR/arriba_contents.txt"

    echo "Validating Arriba archive structure"
    tar -tzf "$archive" > "$listing"

    grep -m 1 'blacklist_hg38_GRCh38.*\.tsv\.gz$' "$listing" >/dev/null || {
        echo "ERROR: Arriba blacklist is absent from $archive" >&2
        return 1
    }

    grep -m 1 'known_fusions_hg38_GRCh38.*\.tsv\.gz$' "$listing" >/dev/null || {
        echo "ERROR: Arriba known-fusions file is absent from $archive" >&2
        return 1
    }

    grep -m 1 'protein_domains_hg38_GRCh38.*\.gff3$' "$listing" >/dev/null || {
        echo "ERROR: Arriba protein-domain file is absent from $archive" >&2
        return 1
    }
}

resume_download \
    'https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
    'Downloading GRCh38 primary assembly'
validate_gzip "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

resume_download \
    'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz' \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.111.gtf.gz" \
    'Downloading Ensembl release 111 GTF'
validate_gzip "$REFERENCE_DIR/Homo_sapiens.GRCh38.111.gtf.gz"

resume_download \
    'https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz' \
    "$REFERENCE_DIR/homo_sapiens_vep_111_GRCh38.tar.gz" \
    'Downloading Ensembl VEP 111 GRCh38 cache'
validate_tar "$REFERENCE_DIR/homo_sapiens_vep_111_GRCh38.tar.gz"
validate_vep_structure "$REFERENCE_DIR/homo_sapiens_vep_111_GRCh38.tar.gz"

resume_download \
    'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29' \
    "$REFERENCE_DIR/human_reviewed_isoforms.fasta.gz" \
    'Downloading reviewed UniProt human proteome with isoforms'
validate_gzip "$REFERENCE_DIR/human_reviewed_isoforms.fasta.gz"

resume_download \
    'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' \
    "$REFERENCE_DIR/arriba_v2.4.0.tar.gz" \
    'Downloading Arriba 2.4.0 resources'
validate_tar "$REFERENCE_DIR/arriba_v2.4.0.tar.gz"
validate_arriba_structure "$REFERENCE_DIR/arriba_v2.4.0.tar.gz"

required_images=(
    'quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img'
    'quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img'
    'quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img'
    'quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img'
    'quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img'
    'quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img'
    'quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img'
    'quay.io-biocontainers-bcftools-1.21--h8b25389_0.img'
)

echo "Validating Singularity images"
for image_name in "${required_images[@]}"; do
    image_path="$CONTAINER_DIR/$image_name"

    if [[ ! -s "$image_path" ]]; then
        echo "ERROR: missing container image: $image_path" >&2
        exit 1
    fi

    singularity inspect "$image_path" >/dev/null
    echo "Validated container: $image_name"
done

echo
echo "All reference archives and container images are complete and valid."
echo "References: $REFERENCE_DIR"
echo "Containers: $CONTAINER_DIR"
