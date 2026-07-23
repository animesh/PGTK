#!/bin/bash
set -uo pipefail

WORKDIR="${1:-$(pwd -P)}"
CONTAINER_DIR="$WORKDIR/singularity_cache"
STAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="$WORKDIR/cli_probe_$STAMP"
REPORT="$WORKDIR/cli_probe_${STAMP}.txt"
mkdir -p "$OUTDIR"
exec > >(tee "$REPORT") 2>&1

run_probe() {
    local label="$1"
    local image="$2"
    shift 2
    local outfile="$OUTDIR/${label}.txt"
    printf '\n===== %s =====\n' "$label"
    printf 'IMAGE: %s\n' "$image"
    printf 'COMMAND:'
    printf ' %q' "$@"
    printf '\n'
    set +e
    singularity exec --bind "$WORKDIR:$WORKDIR" --pwd "$WORKDIR" "$image" "$@" >"$outfile" 2>&1
    local status=$?
    set -e
    printf 'EXIT: %s\n' "$status"
    sed -n '1,240p' "$outfile"
}

SRA="$CONTAINER_DIR/quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img"
STAR="$CONTAINER_DIR/quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img"
VEP="$CONTAINER_DIR/quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img"
PYPGATK="$CONTAINER_DIR/quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img"
TRANSDECODER="$CONTAINER_DIR/transdecoder-6.0.0.img"
PVAC="$CONTAINER_DIR/pvactools-7.1.1.img"

printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Workdir: %s\n' "$WORKDIR"
printf 'Raw outputs: %s\n' "$OUTDIR"
printf 'Combined report: %s\n' "$REPORT"

run_probe fasterq_help "$SRA" fasterq-dump --help
run_probe star_help "$STAR" STAR --help
run_probe star_parameters "$STAR" bash -lc "STAR --help 2>&1 | grep -E 'twopassMode|chimOutType|outSAMtype|sjdbGTFfile|sjdbOverhang'"

run_probe vep_help "$VEP" vep --help
run_probe vep_source_options "$VEP" bash -lc "grep -Rho -- '--dir_cache\|--fasta\|--canonical\|--protein\|--hgvs\|--fork' /usr/local/bin/vep /opt/conda/share/ensembl-vep-* 2>/dev/null | sort -u"

run_probe pypgatk_top_help "$PYPGATK" pypgatk --help
run_probe pypgatk_command_help_long "$PYPGATK" pypgatk vcf-to-proteindb --help
run_probe pypgatk_command_help_short "$PYPGATK" pypgatk vcf-to-proteindb -h
run_probe pypgatk_cli_locations "$PYPGATK" bash -lc "command -v pypgatk; command -v pypgatk_cli; command -v pypgatk_cli.py; compgen -c | grep -E '^pypgatk' | sort -u"
run_probe pypgatk_click_source "$PYPGATK" bash -lc "python - <<'PY'
import inspect
import pypgatk.commands.vcf_to_proteindb as m
print(inspect.getsource(m))
PY"

run_probe transdecoder_locations "$TRANSDECODER" bash -lc "command -v TransDecoder; command -v TransDecoder.LongOrfs; command -v TransDecoder.Predict; command -v gtf_genome_to_cdna_fasta.pl"
run_probe transdecoder_wrapper_help "$TRANSDECODER" TransDecoder --help
run_probe transdecoder_longorfs_noargs "$TRANSDECODER" TransDecoder.LongOrfs
run_probe transdecoder_predict_noargs "$TRANSDECODER" TransDecoder.Predict
run_probe transdecoder_gtf_noargs "$TRANSDECODER" gtf_genome_to_cdna_fasta.pl

run_probe pvacfuse_help "$PVAC" pvacfuse generate_protein_fasta --help

printf '\n===== KEY OPTION SUMMARY =====\n'
printf '\nSTAR:\n'
grep -Eh 'twopassMode|chimOutType|outSAMtype|sjdbGTFfile|sjdbOverhang' "$OUTDIR"/star_*.txt | sort -u || true
printf '\nVEP:\n'
grep -Eh -- '--dir_cache|--fasta|--canonical|--protein|--hgvs|--fork' "$OUTDIR"/vep_*.txt | sort -u || true
printf '\npypgatk:\n'
grep -Eh -- '--vcf|--input[_-]vcf|--input[_-]fasta|--protein[_-]db[_-]fasta|--gene[_-]annotations[_-]gtf|--output[_-]proteindb|--annotation[_-]field[_-]name|--af[_-]field|consequence' "$OUTDIR"/pypgatk_*.txt | sort -u || true
printf '\nTransDecoder:\n'
grep -Eh -- '(^|[[:space:]])-t|--transcripts|--output_dir|(^|[[:space:]])-O|Usage|usage' "$OUTDIR"/transdecoder_*.txt | sort -u || true
printf '\npVACfuse:\n'
grep -Eh -- '--input-type|Arriba|AGFusion|downstream-sequence-length' "$OUTDIR"/pvacfuse_help.txt | sort -u || true

printf '\nREPORT=%s\n' "$REPORT"
printf 'RAW_OUTPUT_DIR=%s\n' "$OUTDIR"
