#!/usr/bin/env bash
set -euo pipefail

# download_fastq_from_project.sh
# Find SRR run accessions for a project (NCBI/SRA) and download FASTQ.gz files via ENA
# Usage: download_fastq_from_project.sh -p PRJNA1176350 -o data/fastq

print_usage() {
    cat <<USAGE
Usage: $0 [options]

Options:
  -p PROJECT     Project accession (e.g. PRJNA1176350) or a comma-separated list of SRR accessions
  -f FILE        File containing one SRR per line (alternative to -p)
  -o OUTDIR      Output directory (default: data/fastq)
  -n             Dry run: only print discovered fastq URLs
  -t THREADS     Number of parallel downloads for aria2c (default: 4)
  -h             Show this help

Examples:
  $0 -p PRJNA1176350 -o data/fastq
  $0 -f srrs.txt -o data/fastq
  $0 -p SRR31089070,SRR31089071 -n
USAGE
}

OUTDIR="data/fastq"
PROJECT=""
SRR_FILE=""
DRY_RUN=0
THREADS=4

while getopts ":p:f:o:t:nh" opt; do
  case ${opt} in
    p ) PROJECT="$OPTARG" ;; 
    f ) SRR_FILE="$OPTARG" ;; 
    o ) OUTDIR="$OPTARG" ;; 
    n ) DRY_RUN=1 ;; 
    t ) THREADS="$OPTARG" ;; 
    h ) print_usage; exit 0 ;; 
    \? ) echo "Invalid option: -$OPTARG" >&2; print_usage; exit 1 ;;
  esac
done

if [[ -z "$PROJECT" && -z "$SRR_FILE" ]]; then
  echo "Error: either -p or -f must be provided" >&2
  print_usage
  exit 1
fi

mkdir -p "$OUTDIR"

tmp_runinfo="$(mktemp)"
tmp_srrs="$(mktemp)"
tmp_urls="$(mktemp)"
trap 'rm -f "$tmp_runinfo" "$tmp_srrs" "$tmp_urls"' EXIT

# If SRR file provided, read those
if [[ -n "$SRR_FILE" ]]; then
  awk 'NF' "$SRR_FILE" | sed 's/\r$//' > "$tmp_srrs"
fi

# If project provided and looks like a project accession (PRJ), fetch run accessions from NCBI SRA
if [[ -n "$PROJECT" ]]; then
  # If PROJECT looks like a comma-separated list of SRRs, split and write
  if [[ "$PROJECT" =~ ^SRR ]]; then
    IFS=',' read -r -a arr <<< "$PROJECT"
    for v in "${arr[@]}"; do echo "$v" >> "$tmp_srrs"; done
  else
    echo "Querying NCBI for project $PROJECT ..."
    # Use esearch+efetch to get runinfo CSV and extract Run column
    curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=${PROJECT}&usehistory=y" \
      | grep -oP '(?<=<QueryKey>)\d+|(?<=<WebEnv>)[^<]+' \
      | { read qk; read we; \
          curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=$qk&WebEnv=$we&rettype=runinfo&retmode=text" -o "$tmp_runinfo"; }
    if [[ ! -s "$tmp_runinfo" ]]; then
      echo "Failed to fetch runinfo from NCBI for project $PROJECT" >&2
      exit 1
    fi
    # extract Run column (first column)
    tail -n +2 "$tmp_runinfo" | cut -d',' -f1 | sed 's/\r$//' | awk 'NF' >> "$tmp_srrs"
  fi
fi

SRR_COUNT=$(wc -l < "$tmp_srrs" | tr -d ' ')
if [[ "$SRR_COUNT" -eq 0 ]]; then
  echo "No SRR accessions found" >&2
  exit 1
fi

echo "Found $SRR_COUNT SRR accessions. Querying ENA for fastq URLs..."

while read -r srr; do
  srr_trim=$(echo "$srr" | tr -d '\r' | sed 's/^\s*//;s/\s*$//')
  [[ -z "$srr_trim" ]] && continue
  # filereport returns header + rows; get fastq_ftp field
  out=$(curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${srr_trim}&result=read_run&fields=fastq_ftp") || true
  # skip header, extract field(s)
  echo "$out" | tail -n +2 | awk -F'\t' '{print $NF}' | tr ';' '\n' | tr ',' '\n' | sed 's/^[[:space:]]*//; s/[[:space:]]*$//' | grep -v '^$' >> "$tmp_urls" || true
done < "$tmp_srrs"

# Normalize URLs to https and clean
awk '!/^$/{u=$0; sub(/^ftp:\/\//, "https://", u); if(u !~ /^https?:\/\//) u="https://" u; print u}' "$tmp_urls" \
  | sed 's/\r//g' | sort -u > "$OUTDIR/ena_fastq_urls.txt"

URL_COUNT=$(wc -l < "$OUTDIR/ena_fastq_urls.txt" | tr -d ' ')
if [[ "$URL_COUNT" -eq 0 ]]; then
  echo "No fastq URLs found via ENA for the given SRRs" >&2
  exit 1
fi

echo "Discovered $URL_COUNT FASTQ URLs. URLs saved to $OUTDIR/ena_fastq_urls.txt"

if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "Dry run requested — printing URLs:" 
  sed -n '1,200p' "$OUTDIR/ena_fastq_urls.txt"
  exit 0
fi

echo "Starting downloads into $OUTDIR"

# Use aria2c if available (parallel), else fallback to wget
if command -v aria2c > /dev/null 2>&1; then
  echo "aria2c detected — using parallel downloads (threads=$THREADS)"
  aria2c -x "$THREADS" -s "$THREADS" -j "$THREADS" -i "$OUTDIR/ena_fastq_urls.txt" -d "$OUTDIR"
else
  echo "aria2c not found — using wget (resumable)"
  wget -c -P "$OUTDIR" -i "$OUTDIR/ena_fastq_urls.txt"
fi

echo "Download finished. Listing files in $OUTDIR:"
ls -lh "$OUTDIR" | sed -n '1,200p'

exit 0
