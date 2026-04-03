#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ─────────────────────────────────────────────────────────────────────────────
// RNA-seq Somatic Variant Calling Pipeline  (FASTQ → SNV/SV → Group comparison)
//
// Input  : samples.csv  (sample, fastq_1, fastq_2, strandedness, TK, Group)
// Groups : resistant  vs  sensitive
//
// Steps:
//   1. STAR_GENOMEGENERATE  — build STAR index from genome.fa + GTF (once)
//   2. TRIM_GALORE          — adapter trimming (auto-detect) + quality/length filter
//   3. STAR_ALIGN           — align trimmed paired-end reads, 2-pass mode
//   4. SAMTOOLS_SORT_INDEX  — coordinate-sort + index the STAR BAM
//   5. MARKDUPLICATES       — mark optical/PCR duplicates (Picard/GATK)
//   6. REMAP_MAPQ           — remap STAR's MAPQ 255 → 60 for GATK compatibility
//   7. BQSR_RECAL           — build per-sample base quality recalibration table (dbSNP)
//   8. BQSR_APPLY           — apply recalibrated quality scores to BAM
//   9. MUTECT2              — tumor-only SNV/indel calling
//  10. LEARN_ORIENTATION    — F1R2 read orientation model
//  11. FILTER_MUTECT2       — apply Mutect2 filters with orientation priors
//  12. MERGE_SNV_GROUP      — merge per-sample VCFs into per-group VCF
//  13. DELLY_CALL           — per-sample structural variant calling
//  14. MERGE_SV_SITES       — merge SV site lists per group
//  15. DELLY_GENOTYPE       — re-genotype each sample against group sites
//  16. MERGE_GENO_GROUP     — merge genotyped samples into final per-group SV VCF
// ─────────────────────────────────────────────────────────────────────────────


// ── 1. Build STAR genome index (runs once) ────────────────────────────────────

process STAR_GENOMEGENERATE {
    tag "star_index"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'

    input:
    path fasta
    path gtf

    output:
    path "star_index"

    script:
    """
    mkdir -p star_index
    STAR \\
        --runMode genomeGenerate \\
        --genomeDir star_index \\
        --genomeFastaFiles ${fasta} \\
        --sjdbGTFfile ${gtf} \\
        --sjdbOverhang ${params.read_length - 1} \\
        --runThreadN ${task.cpus} \\
        --genomeSAindexNbases 14
    """
}


// ── 2. Adapter trimming (Trim Galore, paired-end, 150 bp) ─────────────────────

process TRIM_GALORE {
    tag "${meta.sample}"
    label 'process_medium'
    container 'wave.seqera.io/wt/a56a73714737/wave/build:trim-galore--f48a39e5cc6e0d0d'

    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*_trimming_report.txt'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("${meta.sample}_R1_trimmed.fq.gz"), path("${meta.sample}_R2_trimmed.fq.gz"), emit: reads
    path "*_trimming_report.txt",                                                                       emit: reports

    script:
    """
    # --detect_adapter_for_pe: Trim Galore scans the first 1M read-pairs to
    # identify the adapter automatically (Illumina universal, Nextera, etc.)
    # rather than us hard-coding a single adapter sequence across all samples.
    trim_galore \\
        --paired \\
        --quality 20 \\
        --length 36 \\
        --cores ${task.cpus} \\
        --gzip \\
        --basename ${meta.sample} \\
        ${fastq_1} ${fastq_2}

    # Trim Galore paired output names: <basename>_R1_val_1.fq.gz / <basename>_R2_val_2.fq.gz
    mv ${meta.sample}_R1_val_1.fq.gz ${meta.sample}_R1_trimmed.fq.gz
    mv ${meta.sample}_R2_val_2.fq.gz ${meta.sample}_R2_trimmed.fq.gz
    """
}


// ── 3. STAR 2-pass alignment ──────────────────────────────────────────────────

process STAR_ALIGN {
    tag "${meta.sample}"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path star_index

    output:
    tuple val(meta), path("${meta.sample}.Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR \\
        --runMode alignReads \\
        --genomeDir ${star_index} \\
        --readFilesIn ${fastq_1} ${fastq_2} \\
        --readFilesCommand zcat \\
        --outSAMtype BAM SortedByCoordinate \\
        --outSAMattributes NH HI AS NM MD RG \\
        --outSAMattrRGline ID:${meta.sample} SM:${meta.sample} PL:ILLUMINA LB:${meta.sample} \\
        --twopassMode Basic \\
        --runThreadN ${task.cpus} \\
        --outFileNamePrefix ${meta.sample}. \\
        --outSAMmapqUnique 60 \\
        --outFilterMultimapNmax 1 \\
        --outFilterMismatchNmax 6 \\
        --outFilterMismatchNoverReadLmax 0.04
    """
}


// ── 3. Sort + index BAM ───────────────────────────────────────────────────────

process SAMTOOLS_SORT_INDEX {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.sample}.sorted.bam"), path("${meta.sample}.sorted.bam.bai")

    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.sample}.sorted.bam ${bam}
    samtools index ${meta.sample}.sorted.bam
    """
}


// ── 4. Mark duplicates ────────────────────────────────────────────────────────

process MARKDUPLICATES {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.markdup.bam*'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.markdup.bam"), path("${meta.sample}.markdup.bam.bai"), emit: bam
    path "${meta.sample}.markdup.metrics.txt",                                                    emit: metrics

    script:
    """
    gatk MarkDuplicates \\
        -I ${bam} \\
        -O ${meta.sample}.markdup.bam \\
        -M ${meta.sample}.markdup.metrics.txt \\
        --CREATE_INDEX true \\
        --VALIDATION_STRINGENCY LENIENT
    mv ${meta.sample}.markdup.bai ${meta.sample}.markdup.bam.bai
    """
}


// ── 5. Remap MAPQ 255 → 60 ────────────────────────────────────────────────────
//    STAR's --outSAMmapqUnique 60 above already handles this for new alignments,
//    but we keep this step as a safety net to catch any residual 255s.

process REMAP_MAPQ {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.remapq.bam"), path("${meta.sample}.remapq.bam.bai")

    script:
    """
    samtools view -h ${bam} \\
        | awk 'BEGIN{OFS="\\t"} /^@/{print; next} \$5==255{\$5=60} {print}' \\
        | samtools sort -@ ${task.cpus} -o ${meta.sample}.remapq.bam
    samtools index ${meta.sample}.remapq.bam
    """
}


// ── 6. Base Quality Score Recalibration (BQSR) ───────────────────────────────
//    BaseRecalibrator builds a per-sample recalibration table by comparing
//    observed base qualities at known-variant sites (dbSNP) against expected
//    qualities.  ApplyBQSR rewrites the BAM with corrected quality scores.
//    This corrects systematic per-cycle and per-context Q-score errors that
//    would otherwise cause Mutect2 to mis-estimate variant confidence.
//
//    Requires: params.dbsnp (path to VCF.gz; .tbi index must sit alongside it)
//    Skip:     params.skip_bqsr = true (e.g. when no dbSNP is available)

process BQSR_RECAL {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)
    path dbsnp
    path dbsnp_tbi

    output:
    tuple val(meta), path(bam), path(bai), path("${meta.sample}.recal.table")

    script:
    """
    gatk BaseRecalibrator \\
        -R ${fasta} \\
        -I ${bam} \\
        --known-sites ${dbsnp} \\
        -O ${meta.sample}.recal.table
    """
}

process BQSR_APPLY {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai), path(recal_table)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.bqsr.bam"), path("${meta.sample}.bqsr.bam.bai")

    script:
    """
    gatk ApplyBQSR \\
        -R ${fasta} \\
        -I ${bam} \\
        --bqsr-recal-file ${recal_table} \\
        -O ${meta.sample}.bqsr.bam
    samtools index ${meta.sample}.bqsr.bam
    """
}

// ── Reference indexing (runs once) ───────────────────────────────────────────

process SAMTOOLS_FAIDX {
    tag "faidx"
    label 'process_low'
    container 'quay.io/biocontainers/samtools:1.21--h50ea8bc_0'

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta}.fai")

    script:
    """
    samtools faidx ${fasta}
    """
}

process GATK_DICT {
    tag "dict"
    label 'process_low'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    input:
    tuple path(fasta), path(fai)

    output:
    tuple path(fasta), path(fai), path("${fasta.baseName}.dict")

    script:
    """
    gatk CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict
    """
}


// ── 6. Mutect2 SNV/indel calling ─────────────────────────────────────────────

process MUTECT2 {
    tag "${meta.sample}"
    label 'process_high'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    publishDir "${params.outdir}/mutect2/${meta.group}", mode: 'copy', pattern: '*.{vcf.gz,vcf.gz.tbi,stats}'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.unfiltered.vcf.gz"), path("${meta.sample}.unfiltered.vcf.gz.tbi"), path("${meta.sample}.f1r2.tar.gz"), path("${meta.sample}.unfiltered.vcf.gz.stats")

    script:
    """
    gatk Mutect2 \\
        -R ${fasta} \\
        -I ${bam} \\
        -tumor ${meta.sample} \\
        --f1r2-tar-gz ${meta.sample}.f1r2.tar.gz \\
        --native-pair-hmm-threads ${task.cpus} \\
        --dont-use-soft-clipped-bases \\
        -O ${meta.sample}.unfiltered.vcf.gz
    """
}


// ── 7. Learn read orientation model ──────────────────────────────────────────

process LEARN_ORIENTATION {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    input:
    tuple val(meta), path(vcf), path(tbi), path(f1r2), path(stats)

    output:
    tuple val(meta), path(vcf), path(tbi), path(stats), path("${meta.sample}.artifact-priors.tar.gz")

    script:
    """
    gatk LearnReadOrientationModel \\
        -I ${f1r2} \\
        -O ${meta.sample}.artifact-priors.tar.gz
    """
}


// ── 8. Filter Mutect2 calls ───────────────────────────────────────────────────

process FILTER_MUTECT2 {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    publishDir "${params.outdir}/mutect2/${meta.group}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(stats), path(priors)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.filtered.vcf.gz"), path("${meta.sample}.filtered.vcf.gz.tbi")

    script:
    """
    gatk FilterMutectCalls \\
        -R ${fasta} \\
        -V ${vcf} \\
        --stats ${stats} \\
        --ob-priors ${priors} \\
        -O ${meta.sample}.filtered.vcf.gz
    """
}


// ── 9. Merge per-sample SNV VCFs into per-group VCF ──────────────────────────

process MERGE_SNV_GROUP {
    tag "${group}"
    label 'process_medium'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    publishDir "${params.outdir}/mutect2_merged", mode: 'copy'

    input:
    tuple val(group), path(vcfs), path(tbis)

    output:
    tuple val(group), path("${group}.merged.snv.vcf.gz"), path("${group}.merged.snv.vcf.gz.tbi")

    script:
    def vcf_list = (vcfs instanceof List ? vcfs : [vcfs])
    def vcf_args = vcf_list.collect { v -> v.toString() }.join(' ')
    """
    # If >1 sample: merge into a multi-sample VCF; if only 1: just copy + index
    n_vcfs=\$(echo "${vcf_args}" | wc -w)
    if [ "\${n_vcfs}" -eq 1 ]; then
        cp ${vcf_args} ${group}.merged.snv.vcf.gz
        bcftools index --tbi ${group}.merged.snv.vcf.gz
    else
        bcftools merge \\
            --force-samples \\
            --merge none \\
            -O z \\
            -o ${group}.merged.snv.vcf.gz \\
            ${vcf_args}
        bcftools index --tbi ${group}.merged.snv.vcf.gz
    fi
    """
}


// ── 10. DELLY SV calling ──────────────────────────────────────────────────────

process DELLY_CALL {
    tag "${meta.sample}"
    label 'process_high'
    container 'quay.io/biocontainers/delly:1.2.9--hd63ebec_1'

    publishDir "${params.outdir}/delly/${meta.group}", mode: 'copy', pattern: '*.bcf*'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.sv.bcf"), path("${meta.sample}.sv.bcf.csi")

    script:
    """
    delly call \\
        -g ${fasta} \\
        -o ${meta.sample}.sv.bcf \\
        ${bam}
    """
}


// ── 11. Merge SV site lists per group ─────────────────────────────────────────

process MERGE_SV_SITES {
    tag "${group}"
    label 'process_medium'
    container 'quay.io/biocontainers/delly:1.2.9--hd63ebec_1'

    input:
    tuple val(group), path(bcfs), path(csis)

    output:
    tuple val(group), path("${group}.sites.bcf"), path("${group}.sites.bcf.csi")

    script:
    def bcf_args = (bcfs instanceof List ? bcfs : [bcfs]).collect { b -> b.toString() }.join(' ')
    """
    n_bcfs=\$(echo "${bcf_args}" | wc -w)
    if [ "\${n_bcfs}" -eq 1 ]; then
        cp ${bcf_args} ${group}.sites.bcf
        bcftools index --csi ${group}.sites.bcf
    else
        delly merge \\
            -o ${group}.sites.bcf \\
            ${bcf_args}
    fi
    """
}


// ── 12. Re-genotype each sample against group-merged SV sites ─────────────────

process DELLY_GENOTYPE {
    tag "${meta.sample}"
    label 'process_high'
    container 'quay.io/biocontainers/delly:1.2.9--hd63ebec_1'

    publishDir "${params.outdir}/delly/${meta.group}", mode: 'copy', pattern: '*.geno.bcf*'

    input:
    tuple val(meta), path(bam), path(bai), val(group), path(sites_bcf), path(sites_csi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.geno.bcf"), path("${meta.sample}.geno.bcf.csi")

    script:
    """
    delly call \\
        -g ${fasta} \\
        -v ${sites_bcf} \\
        -o ${meta.sample}.geno.bcf \\
        ${bam}
    """
}


// ── 13. Merge genotyped samples into final per-group SV VCF ──────────────────

process MERGE_GENO_GROUP {
    tag "${group}"
    label 'process_medium'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    publishDir "${params.outdir}/delly_final", mode: 'copy'

    input:
    tuple val(group), path(bcfs), path(csis)

    output:
    tuple val(group), path("${group}.final.sv.vcf.gz"), path("${group}.final.sv.vcf.gz.tbi")

    script:
    def bcf_args = (bcfs instanceof List ? bcfs : [bcfs]).collect { b -> b.toString() }.join(' ')
    """
    n_bcfs=\$(echo "${bcf_args}" | wc -w)
    if [ "\${n_bcfs}" -eq 1 ]; then
        bcftools view -O z -o ${group}.final.sv.vcf.gz ${bcf_args}
        bcftools index --tbi ${group}.final.sv.vcf.gz
    else
        bcftools merge \\
            --force-samples \\
            --merge all \\
            -O z \\
            -o ${group}.final.sv.vcf.gz \\
            ${bcf_args}
        bcftools index --tbi ${group}.final.sv.vcf.gz
    fi
    """
}


// ── WORKFLOW ──────────────────────────────────────────────────────────────────

workflow {

    // ── Parse samples.csv ─────────────────────────────────────────────────────
    // Columns: sample, fastq_1, fastq_2, strandedness, TK, Group
    // meta.group = phenotype group (resistant / sensitive)
    // meta.tk    = TK patient ID
    def ch_reads = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            def meta = [
                sample : row.sample,
                group  : row.Group,
                tk     : row.TK
            ]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }

    // ── Reference channels (singletons) ──────────────────────────────────────
    def ch_fasta  = channel.value(file(params.genome))
    def ch_gtf    = channel.value(file(params.gtf))

    // ── Build STAR index (once) ───────────────────────────────────────────────
    // STAR_GENOMEGENERATE outputs a single directory — wrap as value channel
    // so it broadcasts to all STAR_ALIGN tasks without being consumed.
    def ch_star_index = params.star_index
        ? channel.value(file(params.star_index))
        : STAR_GENOMEGENERATE(ch_fasta, ch_gtf).collect().map { dirs -> dirs[0] }

    // ── Reference genome indexing (once) ─────────────────────────────────────
    def ch_fai = SAMTOOLS_FAIDX(ch_fasta)
    def ch_ref  = GATK_DICT(ch_fai)

    // ── Adapter trimming (Trim Galore, paired-end) ───────────────────────────
    def ch_trimmed = params.skip_trimming
        ? ch_reads
        : TRIM_GALORE(ch_reads).reads

    // ── Alignment → markdup → MAPQ remap ─────────────────────────────────────
    def ch_aligned  = STAR_ALIGN(ch_trimmed, ch_star_index)
    def ch_sorted   = SAMTOOLS_SORT_INDEX(ch_aligned)
    def ch_markdup  = MARKDUPLICATES(ch_sorted).bam
    def ch_remapped = REMAP_MAPQ(ch_markdup)

    // ── BQSR: correct systematic base quality score errors ───────────────────
    // When params.skip_bqsr = true (or params.dbsnp is not set), the remapped
    // BAMs are passed directly to callers without recalibration.
    def ch_bqsr_input = ch_remapped
    def ch_for_calling
    if (!params.skip_bqsr && params.dbsnp) {
        def ch_dbsnp     = channel.value(file(params.dbsnp))
        def ch_dbsnp_tbi = channel.value(file("${params.dbsnp}.tbi"))
        def ch_recal     = BQSR_RECAL(ch_bqsr_input, ch_ref, ch_dbsnp, ch_dbsnp_tbi)
        ch_for_calling   = BQSR_APPLY(ch_recal, ch_ref)
    } else {
        ch_for_calling = ch_remapped
    }

    // ── SNV/Indel calling ─────────────────────────────────────────────────────
    def ch_raw       = MUTECT2(ch_for_calling, ch_ref)
    def ch_oriented  = LEARN_ORIENTATION(ch_raw)
    def ch_filtered  = FILTER_MUTECT2(ch_oriented, ch_ref)

    // Merge per-sample filtered VCFs → per-group VCF
    ch_filtered
        | map { meta, vcf, tbi -> tuple(meta.group, vcf, tbi) }
        | groupTuple(by: 0)
        | MERGE_SNV_GROUP

    // ── SV calling ────────────────────────────────────────────────────────────
    def ch_sv_raw = DELLY_CALL(ch_for_calling, ch_ref)

    // Step 1: merge per-sample SV site lists per group
    def ch_sv_sites = ch_sv_raw
        | map { meta, bcf, csi -> tuple(meta.group, bcf, csi) }
        | groupTuple(by: 0)
        | MERGE_SV_SITES

    // Step 2: re-genotype each sample against its group's merged sites
    def ch_geno_input = ch_for_calling
        | map { meta, bam, bai -> tuple(meta.group, meta, bam, bai) }
        | combine(ch_sv_sites, by: 0)
        | map { group, meta, bam, bai, sites_bcf, sites_csi ->
            tuple(meta, bam, bai, group, sites_bcf, sites_csi)
        }

    def ch_geno = DELLY_GENOTYPE(ch_geno_input, ch_ref)

    // Step 3: merge genotyped samples into final per-group SV VCF
    ch_geno
        | map { meta, bcf, csi -> tuple(meta.group, bcf, csi) }
        | groupTuple(by: 0)
        | MERGE_GENO_GROUP
}
