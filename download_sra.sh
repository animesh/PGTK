#!/bin/bash

set -euo pipefail

WORKDIR="${1:-$(pwd -P)}"
SAMPLESHEET="$WORKDIR/samples.csv"
SRA_DIR="$WORKDIR/sra_cache"
SRA_IMAGE="$WORKDIR/singularity_cache/quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img"
TMP_BASE="$WORKDIR/download_tmp"

#export http_proxy= #"http://proxy.saga:3128/"
#export https_proxy= #"http://proxy.saga:3128/"
#export HTTP_PROXY="$http_proxy"
#export HTTPS_PROXY="$https_proxy"
#export no_proxy="localhost,127.0.0.1"
#export NO_PROXY="$no_proxy"
export SINGULARITY_TMPDIR="$WORKDIR/singularity_tmp"
export APPTAINER_TMPDIR="$SINGULARITY_TMPDIR"
export TMPDIR="$TMP_BASE"

mkdir -p "$SRA_DIR" "$SINGULARITY_TMPDIR" "$TMP_BASE"

test -s "$SAMPLESHEET" || {
    echo "ERROR: samplesheet not found: $SAMPLESHEET" >&2
    exit 1
}

test -s "$SRA_IMAGE" || {
    echo "ERROR: SRA Tools image not found: $SRA_IMAGE" >&2
    exit 1
}

command -v singularity >/dev/null
singularity inspect "$SRA_IMAGE" >/dev/null

mapfile -t SRRS < <(
    awk -F, 'NR > 1 {gsub(/\r/, "", $2); if ($2 != "") print $2}' "$SAMPLESHEET" |
    sort -u
)

if (( ${#SRRS[@]} == 0 )); then
    echo "ERROR: no SRR accessions found in column 2 of $SAMPLESHEET" >&2
    exit 1
fi

for srr in "${SRRS[@]}"; do
    final_dir="$SRA_DIR/$srr"
    final_file="$final_dir/$srr.sra"
    attempt_dir="$SRA_DIR/.${srr}.download"

    if [[ -s "$final_file" ]]; then
        echo "Using existing SRA: $final_file"
        continue
    fi

    echo "Prefetching $srr"
    mkdir -p "$final_dir"

    attempt=0
    while [[ ! -s "$final_file" ]]; do
        attempt=$((attempt + 1))
        if (( attempt > 20 )); then
            echo "ERROR: failed to prefetch $srr after 20 attempts" >&2
            exit 1
        fi

        rm -rf "$attempt_dir"
        mkdir -p "$attempt_dir"

        if singularity exec \
            --bind "$WORKDIR:$WORKDIR" \
            --pwd "$WORKDIR" \
            --env "http_proxy=$http_proxy" \
            --env "https_proxy=$https_proxy" \
            --env "HTTP_PROXY=$HTTP_PROXY" \
            --env "HTTPS_PROXY=$HTTPS_PROXY" \
            "$SRA_IMAGE" \
            prefetch \
                --max-size u \
                --output-directory "$attempt_dir" \
                "$srr"; then

            downloaded_file=$(find "$attempt_dir" -type f -name "$srr.sra" -print -quit)

            if [[ -n "$downloaded_file" && -s "$downloaded_file" ]]; then
                mv "$downloaded_file" "$final_file"
                rm -rf "$attempt_dir"
                break
            fi
        fi

        echo "Retrying $srr after failed or incomplete download, attempt $attempt/20"
        sleep 30
    done

    if ! singularity exec \
        --bind "$WORKDIR:$WORKDIR" \
        --pwd "$WORKDIR" \
        "$SRA_IMAGE" \
        vdb-validate "$final_file"; then
        echo "ERROR: vdb-validate failed for $final_file" >&2
        rm -f "$final_file"
        exit 1
    fi

    echo "Downloaded and validated: $final_file"
done

echo
echo "All SRA archives are available under: $SRA_DIR"
find "$SRA_DIR" -type f -name '*.sra' -printf '%p %s bytes\n' | sort
