#!/bin/bash
#SBATCH --account=nn9036k
#SBATCH --job-name=resultsTKvep
#SBATCH --time=168:00:00
#SBATCH --partition=normal
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --mail-user=animesh.sharma@ntnu.no
#SBATCH --mail-type=ALL
#SBATCH --output=resultsTKvep-%j.log

set -euo pipefail

PROJECT_ID="nn9036k"
NEXTFLOW="/cluster/home/ash022/scripts/nextflow"
WORKDIR="${SLURM_SUBMIT_DIR:-$PWD}"
RUN_ID="$SLURM_JOB_ID"
CONTAINER_DIR="$WORKDIR/singularity_cache"
REFERENCE_DIR="$WORKDIR/reference_downloads"
SRA_DIR="$WORKDIR/sra_cache"

cd "$WORKDIR"

#export http_proxy= "http://proxy.saga:3128/"
#export https_proxy= "http://proxy.saga:3128/"
#export HTTP_PROXY="$http_proxy"
#export HTTPS_PROXY="$https_proxy"
#export no_proxy="localhost,127.0.0.1"
#export NO_PROXY="$no_proxy"

export NXF_HOME="$WORKDIR/.nextflow_home"
export NXF_SINGULARITY_CACHEDIR="$CONTAINER_DIR"
export NXF_WORK="$WORKDIR/work"
export NXF_OPTS="-Xms1g -Xmx6g"
export SINGULARITY_CACHEDIR="$CONTAINER_DIR/oci_cache"
export APPTAINER_CACHEDIR="$SINGULARITY_CACHEDIR"
export SINGULARITY_TMPDIR="$WORKDIR/singularity_tmp"
export APPTAINER_TMPDIR="$SINGULARITY_TMPDIR"
export TMPDIR="$WORKDIR/tmp"

mkdir -p "$NXF_HOME" "$SINGULARITY_TMPDIR" "$TMPDIR" "$NXF_WORK" results

module purge
module load Java/21

command -v singularity
singularity --version
java --version
"$NEXTFLOW" -version

test -s main.nf
test -s samples.csv

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
for image_name in "${required_images[@]}"; do
  image_path="$CONTAINER_DIR/$image_name"
  test -s "$image_path" || { echo "ERROR: missing $image_path" >&2; exit 1; }
  singularity inspect "$image_path" >/dev/null
  echo "Validated container: $image_name"
done

GENOME_ARCHIVE="$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
GTF_ARCHIVE="$REFERENCE_DIR/Homo_sapiens.GRCh38.111.gtf.gz"
VEP_ARCHIVE="$REFERENCE_DIR/homo_sapiens_vep_111_GRCh38.tar.gz"
PROTEOME_ARCHIVE="$REFERENCE_DIR/human_reviewed_isoforms.fasta.gz"
ARRIBA_ARCHIVE="$REFERENCE_DIR/arriba_v2.4.0.tar.gz"
for archive in "$GENOME_ARCHIVE" "$GTF_ARCHIVE" "$VEP_ARCHIVE" "$PROTEOME_ARCHIVE" "$ARRIBA_ARCHIVE"; do
  test -s "$archive" || { echo "ERROR: missing $archive" >&2; exit 1; }
done
echo "Confirmed all local reference archives are present"

mapfile -t REQUIRED_SRRS < <(awk -F, 'NR>1 {gsub(/\r/,"",$2); if ($2!="") print $2}' samples.csv | sort -u)
for srr in "${REQUIRED_SRRS[@]}"; do
  sra_file="$SRA_DIR/$srr/$srr.sra"
  test -s "$sra_file" || {
    echo "ERROR: missing SRA archive: $sra_file" >&2
    echo "Run ./download_sra.sh on the login node first." >&2
    exit 1
  }
done
echo "Validated all local SRA archives"

"$NEXTFLOW" run main.nf \
  -ansi-log false \
  -with-singularity \
  -process.executor slurm \
  "-process.clusterOptions=--account=$PROJECT_ID" \
  -work-dir "$NXF_WORK" \
  --sra_dir "$SRA_DIR" \
  --genome_url "file://$GENOME_ARCHIVE" \
  --gtf_url "file://$GTF_ARCHIVE" \
  --vep_cache_url "file://$VEP_ARCHIVE" \
  --proteome_url "file://$PROTEOME_ARCHIVE" \
  --arriba_url "file://$ARRIBA_ARCHIVE" \
  -with-report "results/pipeline_report-${RUN_ID}.html" \
  -with-timeline "results/pipeline_timeline-${RUN_ID}.html" \
  -with-trace "results/pipeline_trace-${RUN_ID}.tsv" \
  -resume
