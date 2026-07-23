#!/bin/bash
set -uo pipefail

WORKDIR="${1:-$(pwd -P)}"
CONTAINER_DIR="$WORKDIR/singularity_cache"
STAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="$WORKDIR/cli_probe_$STAMP"
REPORT="$WORKDIR/cli_probe_${STAMP}.txt"

mkdir -p "$OUTDIR"
exec > >(tee "$REPORT") 2>&1

probe() {
    local label="$1"
    local image="$2"
    shift 2
    local output="$OUTDIR/${label}.txt"

    printf '\n===== %s =====\n' "$label"
    printf 'IMAGE: %s\n' "$image"
    printf 'COMMAND:'
    printf ' %q' "$@"
    printf '\n'

    if [[ ! -s "$image" ]]; then
        printf 'EXIT: 127\nERROR: missing image: %s\n' "$image"
        return 0
    fi

    set +e
    singularity exec \
        --bind "$WORKDIR:$WORKDIR" \
        --pwd "$WORKDIR" \
        "$image" "$@" >"$output" 2>&1
    local status=$?
    set -e

    printf 'EXIT: %s\n' "$status"
    sed -n '1,260p' "$output"
}

SRA="$CONTAINER_DIR/quay.io-biocontainers-sra-tools-3.2.1--h4304569_0.img"
TRIM="$CONTAINER_DIR/quay.io-biocontainers-trim-galore-0.6.10--hdfd78af_0.img"
FASTQC="$CONTAINER_DIR/quay.io-biocontainers-fastqc-0.12.1--hdfd78af_0.img"
STAR="$CONTAINER_DIR/quay.io-biocontainers-star-2.7.11b--h43eeafb_1.img"
SAMTOOLS="$CONTAINER_DIR/quay.io-biocontainers-samtools-1.21--h96c455f_1.img"
GATK="$CONTAINER_DIR/quay.io-biocontainers-gatk4-4.6.1.0--py310hdfd78af_0.img"
BCFTOOLS="$CONTAINER_DIR/quay.io-biocontainers-bcftools-1.21--h8b25389_0.img"
VEP="$CONTAINER_DIR/quay.io-biocontainers-ensembl-vep-111.0--pl5321h2a3209d_0.img"
PYPGATK="$CONTAINER_DIR/quay.io-biocontainers-pypgatk-0.0.24--pyhdfd78af_0.img"
ARRIBA="$CONTAINER_DIR/quay.io-biocontainers-arriba-2.4.0--h0033a41_2.img"
MULTIQC="$CONTAINER_DIR/quay.io-biocontainers-multiqc-1.35--pyhdfd78af_1.img"
STRINGTIE="$CONTAINER_DIR/stringtie-3.0.3.img"
TRANSDECODER="$CONTAINER_DIR/transdecoder-6.0.0.img"
PVAC="$CONTAINER_DIR/pvactools-7.1.1.img"
GFFCOMPARE="$CONTAINER_DIR/gffcompare-0.12.10.img"

printf 'Started: %s\n' "$(date --iso-8601=seconds)"
printf 'Workdir: %s\n' "$WORKDIR"
printf 'Raw outputs: %s\n' "$OUTDIR"
printf 'Combined report: %s\n' "$REPORT"

probe fasterq_help "$SRA" fasterq-dump --help
probe trim_galore_help "$TRIM" trim_galore --help
probe fastqc_help "$FASTQC" fastqc --help
probe star_help "$STAR" STAR --help
probe samtools_help "$SAMTOOLS" samtools --help
probe samtools_dict_help "$SAMTOOLS" samtools dict --help
probe samtools_faidx_help "$SAMTOOLS" samtools faidx --help
probe gatk_tools "$GATK" gatk --list
probe bcftools_help "$BCFTOOLS" bcftools --help
probe vep_help "$VEP" vep --help
probe pypgatk_help "$PYPGATK" pypgatk vcf-to-proteindb --help
probe pypgatk_source "$PYPGATK" bash --noprofile --norc -c "python - <<'PY'
import inspect
import pypgatk.commands.vcf_to_proteindb as module
print(inspect.getsource(module))
PY"
probe arriba_help "$ARRIBA" arriba -h
probe stringtie_help "$STRINGTIE" stringtie --help
probe transdecoder_files "$TRANSDECODER" find /usr/local/opt/transdecoder -maxdepth 3 -type f
probe transdecoder_longorfs "$TRANSDECODER" /usr/local/opt/transdecoder/util/TransDecoder.LongOrfs
probe transdecoder_predict "$TRANSDECODER" /usr/local/opt/transdecoder/util/TransDecoder.Predict
probe transdecoder_gtf "$TRANSDECODER" /usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl
probe pvacfuse_help "$PVAC" pvacfuse generate_protein_fasta --help
probe gffcompare_help "$GFFCOMPARE" gffcompare --help
probe gffcompare_version "$GFFCOMPARE" gffcompare --version
probe multiqc_help "$MULTIQC" multiqc --help

printf '\n===== KEY SUMMARY =====\n'
printf '\npypgatk:\n'
grep -Eh -- '--input_fasta|--vcf|--gene_annotations_gtf|--output_proteindb|--annotation_field_name|--af_field|--include_consequences' "$OUTDIR"/pypgatk_*.txt | sort -u || true
printf '\nTransDecoder:\n'
grep -Eh -- '-t <string>|-m <int>|--output_dir|gtf_genome_to_cdna_fasta' "$OUTDIR"/transdecoder_*.txt | sort -u || true
printf '\npVACfuse:\n'
grep -Eh -- '--input-type|Arriba|AGFusion|downstream-sequence-length' "$OUTDIR"/pvacfuse_help.txt | sort -u || true
printf '\ngffcompare:\n'
grep -Eh -- 'gffcompare| -r | -o |version' "$OUTDIR"/gffcompare_*.txt | sort -u || true
printf '\nREPORT=%s\nRAW_OUTPUT_DIR=%s\n' "$REPORT" "$OUTDIR"
