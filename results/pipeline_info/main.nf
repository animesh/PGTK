#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ─────────────────────────────────────────────
// Somatic Variant Calling Pipeline
// SNVs/Indels : GATK Mutect2 (tumor-only, RNA-seq)
// SVs          : DELLY2 (tumor-only)
// Per-replicate calling → per-group VCF merge
//
// NOTE: BAMs are STAR-aligned RNA-seq.
//   STAR sets MAPQ=255 for uniquely-mapped reads;
//   GATK treats 255 as "unavailable" and drops them.
//   REMAP_MAPQ reassigns 255→60 before Mutect2.
// ─────────────────────────────────────────────


// ── Index the reference (once) ───────────────

process SAMTOOLS_FAIDX {
    tag "faidx"
    label 'process_low'
    container 'wave.seqera.io/wt/e0cf8d584fad/wave/build:gatk4-4.6.1.0_samtools-1.21--529d27cd48e49b3c'

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
    container 'wave.seqera.io/wt/e0cf8d584fad/wave/build:gatk4-4.6.1.0_samtools-1.21--529d27cd48e49b3c'

    input:
    tuple path(fasta), path(fai)

    output:
    tuple path(fasta), path(fai), path("${fasta.baseName}.dict")

    script:
    """
    gatk CreateSequenceDictionary -R ${fasta} -O ${fasta.baseName}.dict
    """
}

// ── Remap MAPQ=255 → 60 for STAR-aligned BAMs ─
//    STAR uses 255 to mean "uniquely mapped"; GATK
//    treats 255 as "unavailable" and drops all such
//    reads. Remapping to 60 (standard unique-mapping
//    MAPQ) makes them visible to Mutect2.

process REMAP_MAPQ {
    tag "${meta.sample}"
    label 'process_medium'
    container 'wave.seqera.io/wt/e0cf8d584fad/wave/build:gatk4-4.6.1.0_samtools-1.21--529d27cd48e49b3c'

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

// ── SNV/Indel calling with Mutect2 ───────────

process MUTECT2 {
    tag "${meta.sample}"
    label 'process_high'
    container 'wave.seqera.io/wt/e0cf8d584fad/wave/build:gatk4-4.6.1.0_samtools-1.21--529d27cd48e49b3c'

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

process LEARN_ORIENTATION {
    tag "${meta.sample}"
    label 'process_medium'
    container 'wave.seqera.io/wt/e0cf8d584fad/wave/build:gatk4-4.6.1.0_samtools-1.21--529d27cd48e49b3c'

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

process FILTER_MUTECT2 {
    tag "${meta.sample}"
    label 'process_medium'
    container 'wave.seqera.io/wt/e0cf8d584fad/wave/build:gatk4-4.6.1.0_samtools-1.21--529d27cd48e49b3c'

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

// ── Merge per-replicate SNV VCFs per group ────

process MERGE_SNV_GROUP {
    tag "${group}"
    label 'process_medium'
    container 'wave.seqera.io/wt/6a1618608b77/wave/build:bcftools-1.21--9d258438ff2655b9'

    publishDir "${params.outdir}/mutect2_merged", mode: 'copy'

    input:
    tuple val(group), path(vcfs), path(tbis)

    output:
    tuple val(group), path("${group}.merged.snv.vcf.gz"), path("${group}.merged.snv.vcf.gz.tbi")

    script:
    def vcf_list = (vcfs instanceof List ? vcfs : [vcfs])
    def vcf_args = vcf_list.collect { v -> v.toString() }.join(' ')
    """
    # Rename each replicate's sample name to the group name,
    # then concat into a single-sample union VCF.
    mkdir -p reheadered
    for vcf in ${vcf_args}; do
        echo "${group}" > sample_name.txt
        bcftools reheader -s sample_name.txt -o reheadered/\${vcf} \${vcf}
        bcftools index --tbi reheadered/\${vcf}
    done

    bcftools concat \\
        --allow-overlaps \\
        -O z \\
        -o ${group}.merged.snv.vcf.gz \\
        \$(ls reheadered/*.vcf.gz)
    bcftools index --tbi ${group}.merged.snv.vcf.gz
    """
}

// ── SV calling with DELLY ─────────────────────

process DELLY_CALL {
    tag "${meta.sample}"
    label 'process_high'
    container 'wave.seqera.io/wt/629bd2b58e86/wave/build:delly-1.2.9--15db17ca668b2e79'

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

// ── Merge per-replicate SV site lists per group ─

process MERGE_SV_GROUP {
    tag "${group}"
    label 'process_medium'
    container 'wave.seqera.io/wt/629bd2b58e86/wave/build:delly-1.2.9--15db17ca668b2e79'

    input:
    tuple val(group), path(bcfs), path(csis)

    output:
    tuple val(group), path("${group}.sites.bcf"), path("${group}.sites.bcf.csi")

    script:
    def bcf_args = (bcfs instanceof List ? bcfs : [bcfs]).collect { b -> b.toString() }.join(' ')
    """
    delly merge \\
        -o ${group}.sites.bcf \\
        ${bcf_args}
    """
}

// ── Re-genotype each replicate against merged SV sites ──

process DELLY_GENOTYPE {
    tag "${meta.sample}"
    label 'process_high'
    container 'wave.seqera.io/wt/629bd2b58e86/wave/build:delly-1.2.9--15db17ca668b2e79'

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

// ── Final per-group SV VCF (merge genotyped replicates) ─

process MERGE_GENO_GROUP {
    tag "${group}"
    label 'process_medium'
    container 'wave.seqera.io/wt/6a1618608b77/wave/build:bcftools-1.21--9d258438ff2655b9'

    publishDir "${params.outdir}/delly_final", mode: 'copy'

    input:
    tuple val(group), path(bcfs), path(csis)

    output:
    tuple val(group), path("${group}.final.sv.vcf.gz"), path("${group}.final.sv.vcf.gz.tbi")

    script:
    def bcf_args = (bcfs instanceof List ? bcfs : [bcfs]).collect { b -> b.toString() }.join(' ')
    """
    bcftools merge \\
        --force-samples \\
        --merge all \\
        -O z \\
        -o ${group}.final.sv.vcf.gz \\
        ${bcf_args}
    bcftools index --tbi ${group}.final.sv.vcf.gz
    """
}

// ── WORKFLOW ──────────────────────────────────

workflow {

    // Parse samplesheet → channel of [meta, bam, bai]
    def ch_samples = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            def meta = [sample: row.sample, group: row.group]
            tuple(meta, file(row.bam), file(row.bai))
        }

    // Reference genome as a singleton value channel
    def ch_fasta = channel.value(file(params.genome))

    // ── Index reference (runs once) ───────────
    def ch_fai = SAMTOOLS_FAIDX(ch_fasta)
    def ch_ref = GATK_DICT(ch_fai)

    // ── SNV/Indel pipeline ────────────────────
    // Remap MAPQ=255→60 (STAR assigns 255 to uniquely-mapped
    // reads, but GATK drops them as "MAPQ unavailable")
    def ch_remapped    = REMAP_MAPQ(ch_samples)
    def ch_mutect2_raw = MUTECT2(ch_remapped, ch_ref)
    def ch_oriented    = LEARN_ORIENTATION(ch_mutect2_raw)
    def ch_snv         = FILTER_MUTECT2(ch_oriented, ch_ref)

    // Group filtered VCFs by sample group, merge replicates
    ch_snv
        | map { meta, vcf, tbi -> tuple(meta.group, vcf, tbi) }
        | groupTuple(by: 0)
        | MERGE_SNV_GROUP

    // ── SV pipeline ───────────────────────────
    def ch_sv_raw = DELLY_CALL(ch_samples, ch_ref)

    // Step 1: merge per-replicate SV site lists per group
    def ch_sv_sites = ch_sv_raw
        | map { meta, bcf, csi -> tuple(meta.group, bcf, csi) }
        | groupTuple(by: 0)
        | MERGE_SV_GROUP

    // Step 2: re-genotype each replicate against its group's merged sites
    def ch_geno_input = ch_samples
        | map { meta, bam, bai -> tuple(meta.group, meta, bam, bai) }
        | combine(ch_sv_sites, by: 0)
        | map { group, meta, bam, bai, sites_bcf, sites_csi ->
            tuple(meta, bam, bai, group, sites_bcf, sites_csi)
        }

    def ch_geno = DELLY_GENOTYPE(ch_geno_input, ch_ref)

    // Step 3: merge genotyped replicates into final per-group SV VCF
    ch_geno
        | map { meta, bcf, csi -> tuple(meta.group, bcf, csi) }
        | groupTuple(by: 0)
        | MERGE_GENO_GROUP
}
