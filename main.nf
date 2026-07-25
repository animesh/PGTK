#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplesheet = "${projectDir}/samples.csv"
params.outdir = "${projectDir}/results"
params.read_length = 150
params.skip_trimming = false
params.sra_dir = "${projectDir}/sra_cache"
params.fusion_flank_aa = 50
params.splice_min_coverage = 2.5
params.splice_min_junction_reads = 3
params.splice_min_isoform_fraction = 0.05
params.splice_min_protein_aa = 60
params.splice_class_codes = 'j,u'

process DOWNLOAD_REFERENCES {
    tag 'GRCh38_Ensembl111'
    cpus 8; memory '32 GB'; time '12h'; disk '100 GB'; queue 'normal'
    output:
    path 'refs/genome.fa', emit: genome
    path 'refs/genes.gtf', emit: gtf
    path 'refs/cdna.fa', emit: cdna
    path 'refs/human_reviewed_isoforms.fasta', emit: proteome
    path 'refs/vep_cache', emit: vep_cache
    path 'refs/arriba_blacklist.tsv.gz', emit: blacklist
    path 'refs/arriba_known_fusions.tsv.gz', emit: known
    path 'refs/arriba_protein_domains.gff3', emit: domains
    script:
    """
    set -euo pipefail

    REFERENCE_DOWNLOADS='${projectDir}/reference_downloads'

    test -s "\${REFERENCE_DOWNLOADS}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    test -s "\${REFERENCE_DOWNLOADS}/Homo_sapiens.GRCh38.111.gtf.gz"
    test -s "\${REFERENCE_DOWNLOADS}/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    test -s "\${REFERENCE_DOWNLOADS}/human_reviewed_isoforms.fasta.gz"
    test -s "\${REFERENCE_DOWNLOADS}/homo_sapiens_vep_111_GRCh38.tar.gz"
    test -s "\${REFERENCE_DOWNLOADS}/arriba_v2.4.0.tar.gz"

    mkdir -p refs/vep_cache refs/arriba_unpack

    gzip -dc "\${REFERENCE_DOWNLOADS}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" > refs/genome.fa
    gzip -dc "\${REFERENCE_DOWNLOADS}/Homo_sapiens.GRCh38.111.gtf.gz" > refs/genes.gtf
    gzip -dc "\${REFERENCE_DOWNLOADS}/Homo_sapiens.GRCh38.cdna.all.fa.gz" > refs/cdna.fa
    gzip -dc "\${REFERENCE_DOWNLOADS}/human_reviewed_isoforms.fasta.gz" > refs/human_reviewed_isoforms.fasta

    tar -xzf "\${REFERENCE_DOWNLOADS}/homo_sapiens_vep_111_GRCh38.tar.gz" -C refs/vep_cache
    tar -xzf "\${REFERENCE_DOWNLOADS}/arriba_v2.4.0.tar.gz" -C refs/arriba_unpack

    cp \$(find refs/arriba_unpack -type f -name 'blacklist_hg38_GRCh38*.tsv.gz' -print -quit) refs/arriba_blacklist.tsv.gz
    cp \$(find refs/arriba_unpack -type f -name 'known_fusions_hg38_GRCh38*.tsv.gz' -print -quit) refs/arriba_known_fusions.tsv.gz
    cp \$(find refs/arriba_unpack -type f -name 'protein_domains_hg38_GRCh38*.gff3' -print -quit) refs/arriba_protein_domains.gff3

    test -s refs/genome.fa
    test -s refs/genes.gtf
    test -s refs/cdna.fa
    test -s refs/human_reviewed_isoforms.fasta
    test -d refs/vep_cache/homo_sapiens/111_GRCh38
    test -s refs/arriba_blacklist.tsv.gz
    test -s refs/arriba_known_fusions.tsv.gz
    test -s refs/arriba_protein_domains.gff3
    """
}

process SRA_TO_FASTQ {
    tag "${meta.sample}:${srr}"
    cpus 16; memory '32 GB'; time '24h'; disk '150 GB'; queue 'normal'
    container 'quay.io/biocontainers/sra-tools:3.2.1--h4304569_0'
    input: tuple val(meta), val(srr), path(sra_file)
    output: tuple val(meta), path("${srr}_1.fastq.gz"), path("${srr}_2.fastq.gz")
    script:
    """
    set -euo pipefail

    mkdir -p fasterq_tmp

    fasterq-dump \
        --verbose \
        --details \
        --split-files \
        --threads ${task.cpus} \
        --temp fasterq_tmp \
        --outdir . \
        ${sra_file}

    test -s ${srr}_1.fastq
    test -s ${srr}_2.fastq

    gzip -1 ${srr}_1.fastq ${srr}_2.fastq

    test -s ${srr}_1.fastq.gz
    test -s ${srr}_2.fastq.gz

    rm -rf fasterq_tmp
    """
}

process CAT_FASTQ {
    tag "${meta.sample}"
    cpus 2; memory '8 GB'; time '4h'; disk '200 GB'; queue 'normal'
    input: tuple val(meta), path(r1s), path(r2s)
    output: tuple val(meta), path("${meta.sample}_R1.fastq.gz"), path("${meta.sample}_R2.fastq.gz")
    script:
    def r1 = r1s instanceof List ? r1s.sort().join(' ') : r1s
    def r2 = r2s instanceof List ? r2s.sort().join(' ') : r2s
    """
    cat ${r1} > ${meta.sample}_R1.fastq.gz
    cat ${r2} > ${meta.sample}_R2.fastq.gz
    """
}

process FASTQC_RAW {
    tag "${meta.sample}:raw"
    cpus 8; memory '16 GB'; time '8h'; disk '100 GB'; queue 'normal'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc/fastqc_raw", mode:'copy'
    input: tuple val(meta), path(r1), path(r2)
    output:
    path "${meta.sample}.raw_fastqc", emit: qc
    script:
    """
    mkdir ${meta.sample}.raw_fastqc
    fastqc --threads ${task.cpus} --outdir ${meta.sample}.raw_fastqc ${r1} ${r2}
    """
}

process TRIM_GALORE {
    tag "${meta.sample}"
    cpus 8; memory '16 GB'; time '12h'; disk '150 GB'; queue 'normal'
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
    publishDir "${params.outdir}/qc/trim_galore", mode:'copy', pattern:'*_trimming_report.txt'
    input: tuple val(meta), path(r1), path(r2)
    output:
    tuple val(meta), path("${meta.sample}_R1.trimmed.fastq.gz"), path("${meta.sample}_R2.trimmed.fastq.gz"), emit: reads
    path '*_trimming_report.txt', emit: reports
    script:
    """
    trim_galore --paired --quality 20 --length 36 --cores ${task.cpus} --gzip --basename ${meta.sample} ${r1} ${r2}
    mv ${meta.sample}_val_1.fq.gz ${meta.sample}_R1.trimmed.fastq.gz
    mv ${meta.sample}_val_2.fq.gz ${meta.sample}_R2.trimmed.fastq.gz
    """
}

process FASTQC_TRIMMED {
    tag "${meta.sample}:trimmed"
    cpus 8; memory '16 GB'; time '8h'; disk '100 GB'; queue 'normal'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc/fastqc_trimmed", mode:'copy'
    input: tuple val(meta), path(r1), path(r2)
    output:
    path "${meta.sample}.trimmed_fastqc", emit: qc
    script:
    """
    mkdir ${meta.sample}.trimmed_fastqc
    fastqc --threads ${task.cpus} --outdir ${meta.sample}.trimmed_fastqc ${r1} ${r2}
    """
}

process STAR_INDEX {
    tag 'GRCh38_Ensembl111'
    cpus 20; memory '64 GB'; time '12h'; disk '320 GB'; queue 'normal'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'
    input: path genome; path gtf
    output: path 'star_index'
    script:
    """
    mkdir star_index
    STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${genome} --sjdbGTFfile ${gtf} --sjdbOverhang ${params.read_length-1} --runThreadN ${task.cpus}
    """
}

process STAR_ALIGN {
    tag "${meta.sample}"
    cpus 32; memory '256 GB'; time '24h'; disk '320 GB'; queue 'bigmem'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'
    input: tuple val(meta), path(r1), path(r2); path index
    output:
    tuple val(meta), path("${meta.sample}.Aligned.out.bam"), emit: bam
    path "${meta.sample}.Log.final.out", emit: logs
    script:
    """
    STAR \
        --genomeDir ${index} \
        --readFilesIn ${r1} ${r2} \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --twopassMode Basic \
        --outFileNamePrefix ${meta.sample}. \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within \
        --outBAMcompression 0 \
        --outSAMattributes NH HI AS NM MD RG \
        --outSAMattrRGline ID:${meta.sample} SM:${meta.sample} PL:ILLUMINA LB:${meta.sample} \
        --outSAMmapqUnique 60 \
        --outFilterMultimapNmax 50 \
        --outFilterMismatchNmax 6 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --peOverlapNbasesMin 10 \
        --alignSplicedMateMapLminOverLmate 0.5 \
        --alignSJstitchMismatchNmax 5 -1 5 5 \
        --chimSegmentMin 10 \
        --chimOutType WithinBAM HardClip \
        --chimJunctionOverhangMin 10 \
        --chimScoreDropMax 30 \
        --chimScoreJunctionNonGTAG 0 \
        --chimScoreSeparation 1 \
        --chimSegmentReadGapMax 3 \
        --chimMultimapNmax 50
    """
}

process SORT_INDEX_BAM {
    tag "${meta.sample}"
    cpus 20; memory '64 GB'; time '24h'; disk '200 GB'; queue 'normal'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/bam/star", mode:'copy', pattern:'*.Aligned.sortedByCoord.out.bam*'
    input: tuple val(meta), path(bam)
    output: tuple val(meta), path("${meta.sample}.Aligned.sortedByCoord.out.bam"), path("${meta.sample}.Aligned.sortedByCoord.out.bam.bai")
    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.sample}.Aligned.sortedByCoord.out.bam ${bam}
    samtools index -@ ${task.cpus} ${meta.sample}.Aligned.sortedByCoord.out.bam
    """
}

process SAMTOOLS_FLAGSTAT {
    tag "${meta.sample}"
    cpus 4; memory '16 GB'; time '4h'; disk '20 GB'; queue 'normal'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    publishDir "${params.outdir}/qc/flagstat", mode:'copy'
    input: tuple val(meta), path(bam), path(bai)
    output: path "${meta.sample}.flagstat.txt"
    script:
    """
    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.sample}.flagstat.txt
    """
}

process REF_INDEX {
    tag 'GRCh38'; cpus 4; memory '16 GB'; time '4h'; disk '30 GB'; queue 'normal'
    container 'quay.io/biocontainers/samtools:1.21--h96c455f_1'
    input: path genome
    output: tuple path('genome.fa'), path('genome.fa.fai'), path('genome.dict')
    script:
    """
    set -euo pipefail
    samtools faidx genome.fa
    samtools dict -o genome.dict genome.fa
    test -s genome.fa
    test -s genome.fa.fai
    test -s genome.dict
    """
}
process MARK_DUPLICATES {
    tag "${meta.sample}"
    cpus 20; memory '64 GB'; time '24h'; disk '200 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path("${meta.sample}.markdup.bam"), path("${meta.sample}.markdup.bam.bai"), emit: bam
    path "${meta.sample}.metrics.txt", emit: metrics
    script:
    """
    set -euo pipefail
    mkdir -p gatk_tmp
    trap 'rm -rf gatk_tmp' EXIT
    gatk --java-options "-Xms4g -Xmx56g -XX:ParallelGCThreads=20 -Djava.io.tmpdir=\${PWD}/gatk_tmp" \
        MarkDuplicates \
        -I ${bam} \
        -O ${meta.sample}.markdup.bam \
        -M ${meta.sample}.metrics.txt \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY LENIENT \
        --MAX_RECORDS_IN_RAM 1000000
    mv ${meta.sample}.markdup.bai ${meta.sample}.markdup.bam.bai
    test -s ${meta.sample}.markdup.bam
    test -s ${meta.sample}.markdup.bam.bai
    test -s ${meta.sample}.metrics.txt
    """
}

process SPLIT_N_CIGAR {
    tag "${meta.sample}"; cpus 20; memory '64 GB'; time '24h'; disk '200 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: tuple val(meta), path(bam), path(bai); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.split.bam"), path("${meta.sample}.split.bam.bai")
    script:
    """
    set -euo pipefail
    mkdir -p gatk_tmp
    trap 'rm -rf gatk_tmp' EXIT
    df -h .
    gatk --java-options "-Xms4g -Xmx56g -Djava.io.tmpdir=\${PWD}/gatk_tmp" SplitNCigarReads -R ${genome} -I ${bam} -O ${meta.sample}.split.bam --create-output-bam-index true
    if [[ -s ${meta.sample}.split.bai && ! -e ${meta.sample}.split.bam.bai ]]; then
        mv ${meta.sample}.split.bai ${meta.sample}.split.bam.bai
    fi
    test -s ${meta.sample}.split.bam
    test -s ${meta.sample}.split.bam.bai
    """
}

process HAPLOTYPE_CALLER {
    tag "${meta.sample}"; cpus 20; memory '64 GB'; time '48h'; disk '120 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/gvcf", mode:'copy'
    input: tuple val(meta), path(bam), path(bai); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.g.vcf.gz"), path("${meta.sample}.g.vcf.gz.tbi")
    script:
    """
    set -euo pipefail
    mkdir -p gatk_tmp
    trap 'rm -rf gatk_tmp' EXIT
    gatk --java-options "-Xms4g -Xmx56g -Djava.io.tmpdir=\${PWD}/gatk_tmp" HaplotypeCaller -R ${genome} -I ${bam} -O ${meta.sample}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20 --native-pair-hmm-threads ${task.cpus}
    """
}

process GENOTYPE_FILTER {
    tag "${meta.sample}"; cpus 8; memory '32 GB'; time '12h'; disk '40 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_pass", mode:'copy', pattern:'*.pass.vcf.gz*'
    input: tuple val(meta), path(gvcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.pass.vcf.gz"), path("${meta.sample}.pass.vcf.gz.tbi")
    script:
    """
    set -euo pipefail
    mkdir -p gatk_tmp
    trap 'rm -rf gatk_tmp' EXIT
    gatk --java-options "-Xms4g -Xmx28g -Djava.io.tmpdir=\${PWD}/gatk_tmp" GenotypeGVCFs -R ${genome} -V ${gvcf} -O ${meta.sample}.raw.vcf.gz
    gatk --java-options "-Xms4g -Xmx28g -Djava.io.tmpdir=\${PWD}/gatk_tmp" VariantFiltration -R ${genome} -V ${meta.sample}.raw.vcf.gz --window 35 --cluster 3 --filter-expression 'QD < 2.0' --filter-name QD2 --filter-expression 'FS > 30.0' --filter-name FS30 --filter-expression 'MQ < 40.0' --filter-name MQ40 --filter-expression 'MQRankSum < -12.5' --filter-name MQRankSum-12.5 --filter-expression 'ReadPosRankSum < -8.0' --filter-name ReadPos-8 -O ${meta.sample}.filtered.vcf.gz
    gatk --java-options "-Xms4g -Xmx28g -Djava.io.tmpdir=\${PWD}/gatk_tmp" SelectVariants -R ${genome} -V ${meta.sample}.filtered.vcf.gz --exclude-filtered -O ${meta.sample}.pass.vcf.gz
    """
}

process BCFTOOLS_STATS {
    tag "${meta.sample}"
    cpus 2; memory '8 GB'; time '4h'; disk '20 GB'; queue 'normal'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    publishDir "${params.outdir}/qc/bcftools", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi)
    output: path "${meta.sample}.bcftools.stats.txt"
    script:
    """
    bcftools stats ${vcf} > ${meta.sample}.bcftools.stats.txt
    """
}

process VEP_ANNOTATE {
    tag "${meta.sample}"; cpus 20; memory '64 GB'; time '24h'; disk '60 GB'; queue 'normal'
    container 'quay.io/biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0'
    publishDir "${params.outdir}/vep", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); tuple path(genome), path(fai), path(dict); path cache
    output: tuple val(meta), path("${meta.sample}.vep.vcf.gz"), path("${meta.sample}.vep.vcf.gz.tbi")
    script:
    """
    vep --input_file ${vcf} --output_file ${meta.sample}.vep.vcf --format vcf --vcf --cache --offline --cache_version 111 --dir_cache ${cache} --species homo_sapiens --assembly GRCh38 --fasta ${genome} --canonical --protein --symbol --numbers --biotype --total_length --hgvs --fork ${task.cpus} --force_overwrite
    bgzip ${meta.sample}.vep.vcf && tabix -p vcf ${meta.sample}.vep.vcf.gz
    """
}

process PYPGATK_FASTA {
    tag "${meta.sample}"; cpus 8; memory '32 GB'; time '12h'; disk '80 GB'; queue 'normal'
    container 'quay.io/biocontainers/pypgatk:0.0.24--pyhdfd78af_0'
    publishDir "${params.outdir}/variant_fasta", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); path gtf; path cdna
    output: tuple val(meta), path("${meta.sample}.variant_proteins.fasta")
    script:
    """
    set -euo pipefail

    gzip -t ${vcf}
    gzip -dc ${vcf} > ${meta.sample}.pypgatk.vcf
    test -s ${meta.sample}.pypgatk.vcf
    grep -q '^##fileformat=VCF' ${meta.sample}.pypgatk.vcf

    if ! grep -qv '^#' ${meta.sample}.pypgatk.vcf; then
        echo "Warning: no VCF records for ${meta.sample}" >&2
        : > ${meta.sample}.variant_proteins.fasta
        exit 0
    fi

    pypgatk vcf-to-proteindb \
        --vcf ${meta.sample}.pypgatk.vcf \
        --input_fasta ${cdna} \
        --gene_annotations_gtf ${gtf} \
        --annotation_field_name CSQ \
        --af_field AF \
        --include_consequences missense_variant,frameshift_variant,stop_gained,stop_lost,start_lost,splice_donor_variant,splice_acceptor_variant,inframe_insertion,inframe_deletion \
        --output_proteindb ${meta.sample}.variant_proteins.fasta

    if [[ ! -s ${meta.sample}.variant_proteins.fasta ]]; then
        echo "Warning: no variant proteins generated for ${meta.sample}" >&2
        : > ${meta.sample}.variant_proteins.fasta
    fi
    sed -i 's/^>/>${meta.sample}|/' ${meta.sample}.variant_proteins.fasta
    """
}

process ARRIBA {
    tag "${meta.sample}"; cpus 8; memory '32 GB'; time '12h'; disk '60 GB'; queue 'normal'
    container 'quay.io/biocontainers/arriba:2.4.0--h0033a41_2'
    publishDir "${params.outdir}/fusions", mode:'copy'
    input: tuple val(meta), path(bam); path genome; path gtf; path blacklist; path known; path domains
    output:
    tuple val(meta), path("${meta.sample}.fusions.tsv"), emit: accepted
    path "${meta.sample}.fusions.discarded.tsv", emit: discarded
    script:
    """
    arriba -x ${bam} -a ${genome} -g ${gtf} -b ${blacklist} -k ${known} -p ${domains} -o ${meta.sample}.fusions.tsv -O ${meta.sample}.fusions.discarded.tsv
    """
}

process FUSION_FASTA {
    tag "${meta.sample}"
    cpus 4; memory '16 GB'; time '8h'; disk '20 GB'; queue 'normal'
    container "${projectDir}/singularity_cache/pvactools-7.1.1.img"
    publishDir "${params.outdir}/fusion_fasta", mode:'copy'
    input: tuple val(meta), path(fusions)
    output: tuple val(meta), path("${meta.sample}.fusion_proteins.fasta")
    script:
    """
    set -euo pipefail
    if awk 'END { exit !(NR > 1) }' ${fusions}; then
        pvacfuse generate_protein_fasta \\
            --downstream-sequence-length full \\
            ${fusions} \\
            ${params.fusion_flank_aa} \\
            ${meta.sample}.fusion_proteins.fasta
        sed -i 's/^>/>${meta.sample}|FUSION|/' ${meta.sample}.fusion_proteins.fasta
    else
        : > ${meta.sample}.fusion_proteins.fasta
    fi
    """
}

process STRINGTIE_ASSEMBLY {
    tag "${meta.sample}"
    cpus 20; memory '64 GB'; time '24h'; disk '80 GB'; queue 'normal'
    container "${projectDir}/singularity_cache/stringtie-3.0.3.img"
    publishDir "${params.outdir}/splicing/stringtie", mode:'copy'
    input: tuple val(meta), path(bam), path(bai); path gtf
    output: tuple val(meta), path("${meta.sample}.assembled.gtf")
    script:
    """
    set -euo pipefail
    stringtie ${bam} \\
        -G ${gtf} \\
        -p ${task.cpus} \\
        -c ${params.splice_min_coverage} \\
        -f ${params.splice_min_isoform_fraction} \\
        -j ${params.splice_min_junction_reads} \\
        -o ${meta.sample}.assembled.gtf
    test -s ${meta.sample}.assembled.gtf
    """
}

process GFFCOMPARE_NOVEL {
    tag "${meta.sample}"
    cpus 4; memory '16 GB'; time '8h'; disk '30 GB'; queue 'normal'
    container "${projectDir}/singularity_cache/gffcompare-0.12.10.img"
    publishDir "${params.outdir}/splicing/gffcompare", mode:'copy'
    input:
    tuple val(meta), path(assembled_gtf)
    path reference_gtf
    output:
    tuple val(meta), path("${meta.sample}.novel.gtf"), emit: novel
    tuple val(meta), path("${meta.sample}.gffcompare.annotated.gtf"), emit: annotated
    path "${meta.sample}.gffcompare.stats", emit: stats
    script:
    """
    set -euo pipefail

    prefix=${meta.sample}.gffcompare

    gffcompare \
        -r ${reference_gtf} \
        -o \${prefix} \
        ${assembled_gtf}

    if [[ ! -s \${prefix}.annotated.gtf ]]; then
        echo "ERROR: gffcompare did not create a non-empty annotated GTF" >&2
        find . -maxdepth 1 -type f -printf '%f %s bytes\n' | sort >&2
        exit 1
    fi

    if [[ -s \${prefix}.stats ]]; then
        cp \${prefix}.stats ${meta.sample}.gffcompare.stats
    else
        printf 'No non-empty gffcompare statistics report was produced for %s\n' '${meta.sample}' \
            > ${meta.sample}.gffcompare.stats
    fi

    awk -v allowed='${params.splice_class_codes}' '
        BEGIN {
            count = split(allowed, values, ",")
            for (i = 1; i <= count; i++) keep_code[values[i]] = 1
        }
        \$3 == "transcript" {
            transcript = ""
            code = ""
            if (match(\$0, /transcript_id "[^"]+"/))
                transcript = substr(\$0, RSTART + 15, RLENGTH - 16)
            if (match(\$0, /class_code "[^"]+"/))
                code = substr(\$0, RSTART + 12, RLENGTH - 13)
            selected[transcript] = keep_code[code]
        }
        {
            transcript = ""
            if (match(\$0, /transcript_id "[^"]+"/))
                transcript = substr(\$0, RSTART + 15, RLENGTH - 16)
            if (selected[transcript]) print
        }
    ' ${meta.sample}.gffcompare.annotated.gtf > ${meta.sample}.novel.gtf

    test -e ${meta.sample}.novel.gtf
    test -s ${meta.sample}.gffcompare.annotated.gtf
    test -s ${meta.sample}.gffcompare.stats
    """
}

process SPLICE_PROTEIN_FASTA {
    tag "${meta.sample}"
    cpus 20; memory '64 GB'; time '24h'; disk '60 GB'; queue 'normal'
    container "${projectDir}/singularity_cache/transdecoder-6.0.0.img"
    publishDir "${params.outdir}/splice_fasta", mode:'copy'
    input: tuple val(meta), path(novel_gtf); path genome
    output: tuple val(meta), path("${meta.sample}.splice_proteins.fasta")
    script:
    """
    set -euo pipefail

    if [[ ! -s ${novel_gtf} ]]; then
        : > ${meta.sample}.splice_proteins.fasta
        exit 0
    fi

    /usr/local/opt/transdecoder/util/gtf_genome_to_cdna_fasta.pl \\
        ${novel_gtf} \\
        ${genome} \\
        > ${meta.sample}.transcripts.fasta

    if [[ ! -s ${meta.sample}.transcripts.fasta ]]; then
        : > ${meta.sample}.splice_proteins.fasta
        exit 0
    fi

    /usr/local/opt/transdecoder/util/TransDecoder.LongOrfs \\
        -t ${meta.sample}.transcripts.fasta \\
        -m ${params.splice_min_protein_aa} \\
        --output_dir transdecoder

    /usr/local/opt/transdecoder/util/TransDecoder.Predict \\
        -t ${meta.sample}.transcripts.fasta \\
        --output_dir transdecoder

    peptide_file=\$(find transdecoder -type f -name '*.transdecoder.pep' -print -quit)

    if [[ -z "\${peptide_file}" || ! -s "\${peptide_file}" ]]; then
        : > ${meta.sample}.splice_proteins.fasta
        exit 0
    fi

    awk -v sample='${meta.sample}' -v min_aa='${params.splice_min_protein_aa}' '
        /^>/ {
            if (sequence != "" && length(sequence) >= min_aa && !seen[sequence]++) {
                print ">" sample "|SPLICE|" header
                print sequence
            }
            header = substr(\$0, 2)
            sequence = ""
            next
        }
        { sequence = sequence \$0 }
        END {
            if (sequence != "" && length(sequence) >= min_aa && !seen[sequence]++) {
                print ">" sample "|SPLICE|" header
                print sequence
            }
        }
    ' "\${peptide_file}" > ${meta.sample}.splice_proteins.fasta
    """
}

process COMBINE_PROTEIN_FASTA {
    tag "${meta.sample}"
    cpus 2; memory '8 GB'; time '4h'; disk '30 GB'; queue 'normal'
    publishDir "${params.outdir}/combined_fasta", mode:'copy'
    input:
    tuple val(meta), path(variant_fasta), path(fusion_fasta), path(splice_fasta)
    path proteome
    output: tuple val(meta), path("${meta.sample}.exploratory_proteogenomics.fasta")
    script:
    """
    set -euo pipefail
    cat ${proteome} ${variant_fasta} ${fusion_fasta} ${splice_fasta} \\
        > ${meta.sample}.combined.raw.fasta

    awk '
        /^>/ {
            if (sequence != "" && !seen[sequence]++) {
                print header
                print sequence
            }
            header = \$0
            sequence = ""
            next
        }
        { sequence = sequence \$0 }
        END {
            if (sequence != "" && !seen[sequence]++) {
                print header
                print sequence
            }
        }
    ' ${meta.sample}.combined.raw.fasta \\
        > ${meta.sample}.exploratory_proteogenomics.fasta

    test -s ${meta.sample}.exploratory_proteogenomics.fasta
    """
}

process PROGRESSION_SUBTRACT {
    tag "${meta.sample}"; cpus 2; memory '8 GB'; time '4h'; disk '20 GB'; queue 'normal'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    publishDir "${params.outdir}/progression_vcf", mode:'copy'
    input: tuple val(meta), path(pv), path(pt), path(bv), path(bt)
    output: tuple val(meta), path("${meta.sample}.progression.vep.vcf.gz"), path("${meta.sample}.progression.vep.vcf.gz.tbi")
    script:
    """
    bcftools isec -C -w 1 -O z -o ${meta.sample}.progression.vep.vcf.gz ${pv} ${bv}
    bcftools index --tbi ${meta.sample}.progression.vep.vcf.gz
    """
}

process PROGRESSION_FASTA {
    tag "${meta.sample}"; cpus 8; memory '32 GB'; time '12h'; disk '80 GB'; queue 'normal'
    container 'quay.io/biocontainers/pypgatk:0.0.24--pyhdfd78af_0'
    publishDir "${params.outdir}/progression_fasta", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); path gtf; path cdna
    output: path "${meta.sample}.progression_proteins.fasta"
    script:
    """
    set -euo pipefail

    gzip -t ${vcf}
    gzip -dc ${vcf} > ${meta.sample}.progression.pypgatk.vcf
    test -s ${meta.sample}.progression.pypgatk.vcf
    grep -q '^##fileformat=VCF' ${meta.sample}.progression.pypgatk.vcf

    if ! grep -qv '^#' ${meta.sample}.progression.pypgatk.vcf; then
        echo "Warning: no progression VCF records for ${meta.sample}" >&2
        : > ${meta.sample}.progression_proteins.fasta
        exit 0
    fi

    pypgatk vcf-to-proteindb \
        --vcf ${meta.sample}.progression.pypgatk.vcf \
        --input_fasta ${cdna} \
        --gene_annotations_gtf ${gtf} \
        --annotation_field_name CSQ \
        --af_field AF \
        --include_consequences missense_variant,frameshift_variant,stop_gained,stop_lost,start_lost,splice_donor_variant,splice_acceptor_variant,inframe_insertion,inframe_deletion \
        --output_proteindb ${meta.sample}.progression_proteins.fasta

    if [[ ! -s ${meta.sample}.progression_proteins.fasta ]]; then
        echo "Warning: no progression proteins generated for ${meta.sample}" >&2
        : > ${meta.sample}.progression_proteins.fasta
    fi
    sed -i 's/^>/>${meta.sample}|PROGRESSION|/' ${meta.sample}.progression_proteins.fasta
    """
}

process MULTIQC {
    tag 'final_report'
    cpus 8; memory '32 GB'; time '4h'; disk '40 GB'; queue 'normal'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/multiqc", mode:'copy'
    input: path qc_files
    output:
    path 'multiqc_report.html'
    path 'multiqc_report_data'
    script:
    """
    multiqc . \
        --force \
        --title 'PGTK RNA Proteogenomics QC' \
        --filename multiqc_report.html \
        --data-dir \
        --data-format tsv
    """
}

workflow {
    samples = channel.fromPath(params.samplesheet, checkIfExists:true).splitCsv(header:true).map { row ->
        if (!row.sample || !row.srr || !row.TK || !row.Group) error 'samples.csv requires sample,srr,TK,Group,baseline'
        def meta=[sample:row.sample.trim(),tk:row.TK.trim(),group:row.Group.trim(),baseline:(row.baseline?:'').trim().toLowerCase()]
        def srr = row.srr.trim()
        def sraFile = file("${params.sra_dir}/${srr}/${srr}.sra", checkIfExists: true)
        tuple(meta, srr, sraFile)
    }
    refs=DOWNLOAD_REFERENCES()
    ref=REF_INDEX(refs.genome)
    staridx=STAR_INDEX(refs.genome,refs.gtf)
    downloaded=SRA_TO_FASTQ(samples)
    reads=downloaded.map { m,r1,r2 -> tuple(m.sample,m,r1,r2) }.groupTuple(by:0).map { id,ms,r1s,r2s -> tuple(ms[0],r1s,r2s) } | CAT_FASTQ
    raw_qc=FASTQC_RAW(reads)
    trim_result=TRIM_GALORE(reads)
    trimmed=params.skip_trimming ? reads : trim_result.reads
    trimmed_qc=FASTQC_TRIMMED(trimmed)
    star_result=STAR_ALIGN(trimmed,staridx)
    arriba_result=ARRIBA(star_result.bam,refs.genome,refs.gtf,refs.blacklist,refs.known,refs.domains)
    fusion_fasta=FUSION_FASTA(arriba_result.accepted)
    sortedbam=SORT_INDEX_BAM(star_result.bam)
    assembled=STRINGTIE_ASSEMBLY(sortedbam,refs.gtf)
    novel_result=GFFCOMPARE_NOVEL(assembled,refs.gtf)
    splice_fasta=SPLICE_PROTEIN_FASTA(novel_result.novel,refs.genome)
    flagstat=SAMTOOLS_FLAGSTAT(sortedbam)
    md_result=MARK_DUPLICATES(sortedbam)
    split=SPLIT_N_CIGAR(md_result.bam,ref)
    gvcf=HAPLOTYPE_CALLER(split,ref)
    pass=GENOTYPE_FILTER(gvcf,ref)
    variant_stats=BCFTOOLS_STATS(pass)
    annotated=VEP_ANNOTATE(pass,ref,refs.vep_cache)
    variant_fasta=PYPGATK_FASTA(annotated,refs.gtf,refs.cdna)

    variant_keyed=variant_fasta.map { m,f -> tuple(m.sample,m,f) }
    fusion_keyed=fusion_fasta.map { m,f -> tuple(m.sample,m,f) }
    splice_keyed=splice_fasta.map { m,f -> tuple(m.sample,m,f) }
    combined_inputs=variant_keyed
        .join(fusion_keyed)
        .join(splice_keyed)
        .map { sample,m1,vf,m2,ff,m3,sf -> tuple(m1,vf,ff,sf) }
    COMBINE_PROTEIN_FASTA(combined_inputs,refs.proteome)
    groups=annotated.branch { m,v,t -> baseline:m.baseline=='true'; progression:m.baseline=='false'; other:true }
    bases=groups.baseline.map { m,v,t -> tuple(m.tk,v,t) }
    pairs=groups.progression.map { m,v,t -> tuple(m.tk,m,v,t) }.combine(bases,by:0).map { k,m,pv,pt,bv,bt -> tuple(m,pv,pt,bv,bt) }
    prog=PROGRESSION_SUBTRACT(pairs)
    PROGRESSION_FASTA(prog,refs.gtf,refs.cdna)

    qc_files = raw_qc.qc
        .mix(trimmed_qc.qc, trim_result.reports, star_result.logs,
             flagstat, md_result.metrics, variant_stats)
        .collect()
    MULTIQC(qc_files)
}
