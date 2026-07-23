#!/bin/bash
set -euo pipefail

WORKDIR="${1:-$(pwd -P)}"
CONTAINER_DIR="$WORKDIR/singularity_cache"
REFERENCE_DIR="$WORKDIR/reference_downloads"
TMP_DIR="$WORKDIR/download_tmp"

#export http_proxy="http://proxy.saga:3128/"
#export https_proxy="http://proxy.saga:3128/"
#export HTTP_PROXY="$http_proxy"
#export HTTPS_PROXY="$https_proxy"

export SINGULARITY_CACHEDIR="$CONTAINER_DIR/oci_cache"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
export SINGULARITY_TMPDIR="$WORKDIR/singularity_tmp"
export APPTAINER_TMPDIR="$SINGULARITY_TMPDIR"
export APPTAINER_QUIET=true
export SINGULARITY_QUIET=true
export TMPDIR="$TMP_DIR"

mkdir -p \
    "$CONTAINER_DIR" \
    "$SINGULARITY_CACHEDIR" \
    "$SINGULARITY_TMPDIR" \
    "$REFERENCE_DIR" \
    "$TMP_DIR"

command -v singularity >/dev/null
command -v curl >/dev/null
command -v gzip >/dev/null
command -v tar >/dev/null
command -v python >/dev/null

pull_image() {
    local name="$1"
    local uri="$2"
    local output="$CONTAINER_DIR/$name"
    local temporary="$output.tmp"
    local log_file="$TMP_DIR/${name}.pull.log"

    if [[ -s "$output" ]] && singularity inspect "$output" >/dev/null 2>&1; then
        echo "Using validated image: $output"
        return 0
    fi

    rm -f "$output" "$temporary"

    for attempt in 1 2 3 4 5; do
        echo "Pull attempt ${attempt}/5: $uri"
        rm -f "$temporary"
        : > "$log_file"

        if singularity --quiet pull --force "$temporary" "$uri" >"$log_file" 2>&1 && \
           singularity inspect "$temporary" >/dev/null 2>&1; then
            mv "$temporary" "$output"
            rm -f "$log_file"
            echo "Validated image: $output"
            return 0
        fi

        echo "Pull attempt ${attempt}/5 failed; preserving OCI cache for retry" >&2
        tail -n 40 "$log_file" >&2 || true
        rm -f "$temporary"
        sleep $((attempt * 30))
    done

    echo "ERROR: failed to pull $uri" >&2
    echo "Last pull log: $log_file" >&2
    return 1
}

resolve_quay_tag() {
    local repository="$1"
    local version_prefix="$2"
    local fallback_tag="${3:-}"
    local response_file="$TMP_DIR/${repository}-${version_prefix}.quay-tags.json"
    local resolved_tag=""

    for attempt in 1 2 3 4 5; do
        echo "Resolve Quay tag attempt ${attempt}/5: ${repository} ${version_prefix}" >&2
        rm -f "$response_file"

        if curl --fail --silent --show-error --location \
            --retry 5 \
            --retry-all-errors \
            --retry-delay 10 \
            --connect-timeout 60 \
            --speed-time 120 \
            --speed-limit 1024 \
            --output "$response_file" \
            "https://quay.io/api/v1/repository/biocontainers/${repository}/tag/?onlyActiveTags=true&limit=100"; then

            resolved_tag=$(
                python - "$version_prefix" "$response_file" <<'PYTAG'
import json
import sys

prefix = sys.argv[1] + "--"
with open(sys.argv[2], encoding="utf-8") as handle:
    payload = json.load(handle)

tags = [
    item["name"]
    for item in payload.get("tags", [])
    if item.get("name", "").startswith(prefix)
]

if not tags:
    raise SystemExit(1)

print(sorted(tags)[-1])
PYTAG
            ) || resolved_tag=""

            if [[ -n "$resolved_tag" ]]; then
                rm -f "$response_file"
                printf '%s\n' "$resolved_tag"
                return 0
            fi
        fi

        sleep $((attempt * 20))
    done

    rm -f "$response_file"

    if [[ -n "$fallback_tag" ]]; then
        echo "WARNING: Quay API unavailable; using pinned fallback tag: $fallback_tag" >&2
        printf '%s\n' "$fallback_tag"
        return 0
    fi

    echo "ERROR: failed to resolve a Quay tag for ${repository} ${version_prefix}" >&2
    return 1
}

pull_biocontainer_stable() {
    local stable_name="$1"
    local repository="$2"
    local version_prefix="$3"
    local fallback_tag="$4"
    local output="$CONTAINER_DIR/$stable_name"
    local resolved_tag

    if [[ -s "$output" ]] && singularity inspect "$output" >/dev/null 2>&1; then
        echo "Using validated image: $output"
        return 0
    fi

    resolved_tag=$(resolve_quay_tag "$repository" "$version_prefix" "$fallback_tag")
    pull_image "$stable_name" "docker://quay.io/biocontainers/${repository}:${resolved_tag}"
}

validate_gzip() {
    gzip -t "$1" >/dev/null 2>&1
}

validate_tar_gzip() {
    tar -tzf "$1" >/dev/null 2>&1
}

download_asset() {
    local url="$1"
    local output="$2"
    local label="$3"
    local validator="$4"
    local resumable="$5"
    local temporary="${output}.download"
    local status=0

    if [[ -s "$output" ]] && "$validator" "$output"; then
        echo "Using validated reference: $output"
        return 0
    fi

    if [[ -e "$output" ]]; then
        echo "Existing $label is incomplete or invalid"
    fi

    for attempt in $(seq 1 500); do
        if [[ "$resumable" == true ]]; then
            echo "$label, attempt ${attempt}/500, existing bytes: $(stat -c %s "$output" 2>/dev/null || echo 0)"
            status=0

            curl --fail --location --continue-at - \
                --retry 5 \
                --retry-all-errors \
                --retry-delay 10 \
                --connect-timeout 60 \
                --speed-time 120 \
                --speed-limit 1024 \
                --output "$output" \
                "$url" || status=$?

            if (( status == 0 )) && "$validator" "$output"; then
                echo "Downloaded and validated: $output"
                return 0
            fi

            if (( status == 33 )); then
                echo "$label server does not support byte ranges; restarting without resume"
                rm -f "$output"
                resumable=false
            fi
        else
            echo "$label, fresh attempt ${attempt}/500"
            rm -f "$temporary"
            status=0

            curl --fail --location \
                --retry 5 \
                --retry-all-errors \
                --retry-delay 10 \
                --connect-timeout 60 \
                --speed-time 120 \
                --speed-limit 1024 \
                --output "$temporary" \
                "$url" || status=$?

            if (( status == 0 )) && "$validator" "$temporary"; then
                mv "$temporary" "$output"
                echo "Downloaded and validated: $output"
                return 0
            fi

            rm -f "$temporary"
        fi

        sleep 20
    done

    echo "ERROR: failed to download and validate $url" >&2
    return 1
}

pull_image \
    'quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img' \
    'docker://quay.io/biocontainers/sra-tools:3.2.1--h4304569_0'

pull_image \
    'quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img' \
    'docker://quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'

pull_image \
    'quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img' \
    'docker://quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'

pull_image \
    'quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img' \
    'docker://quay.io/biocontainers/star:2.7.11b--h43eeafb_1'

pull_image \
    'quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img' \
    'docker://quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

pull_image \
    'quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img' \
    'docker://quay.io/biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0'

pull_image \
    'quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img' \
    'docker://quay.io/biocontainers/pypgatk:0.0.24--pyhdfd78af_0'

pull_image \
    'quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img' \
    'docker://quay.io/biocontainers/arriba:2.4.0--h0033a41_2'

pull_image \
    'quay.io-biocontainers-bcftools-1.21--h8b25389_0.img' \
    'docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0'

pull_image \
    'quay.io-biocontainers-samtools-1.21--h96c455f_1.img' \
    'docker://quay.io/biocontainers/samtools:1.21--h96c455f_1'

pull_image \
    'quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img' \
    'docker://quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'

pull_biocontainer_stable \
    'stringtie-3.0.3.img' \
    'stringtie' \
    '3.0.3' \
    '3.0.3--h29c0135_0'

pull_biocontainer_stable \
    'transdecoder-6.0.0.img' \
    'transdecoder' \
    '6.0.0' \
    '6.0.0--pl5321hdfd78af_0'

pull_biocontainer_stable \
    'gffcompare-0.12.10.img' \
    'gffcompare' \
    '0.12.10' \
    '0.12.10--h9948957_0'

pull_image \
    'pvactools-7.1.1.img' \
    'docker://griffithlab/pvactools:7.1.1'

download_asset \
    'https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz' \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" \
    'GRCh38 genome' \
    validate_gzip \
    true

download_asset \
    'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz' \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.111.gtf.gz" \
    'Ensembl GTF' \
    validate_gzip \
    true

download_asset \
    'https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz' \
    "$REFERENCE_DIR/Homo_sapiens.GRCh38.cdna.all.fa.gz" \
    'Ensembl cDNA transcripts' \
    validate_gzip \
    true

download_asset \
    'https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz' \
    "$REFERENCE_DIR/homo_sapiens_vep_111_GRCh38.tar.gz" \
    'VEP cache' \
    validate_tar_gzip \
    true

# UniProt streaming responses do not support HTTP range requests.
download_asset \
    'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29' \
    "$REFERENCE_DIR/human_reviewed_isoforms.fasta.gz" \
    'UniProt proteome' \
    validate_gzip \
    false

download_asset \
    'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz' \
    "$REFERENCE_DIR/arriba_v2.4.0.tar.gz" \
    'Arriba resources' \
    validate_tar_gzip \
    true

echo "All 15 containers and all reference archives are valid."
