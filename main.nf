#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ---------------------------------------------------------------------------
// RNA-seq Somatic Variant Calling + Variant Proteomics Pipeline
//
// Input  : samples.csv  (sample, srr, TK, Group, baseline)
//          One row per SRR accession. Multiple rows with the same sample name
//          are automatically merged after download (handles triplicates).
//
//          sample   : cell line name (TK12, TK13, TK14, ...)
//          srr      : SRA run accession (SRR31089074, ...)
//          TK       : patient ID linking longitudinal samples (e.g. patient1)
//          Group    : sensitive | resistant
//          baseline : true  = earliest timepoint per patient (TK12)
//                     false = later timepoints (TK13, TK14)
//                     ""    = single-timepoint samples (TK9, TK10, TK16, TK18)
//
// Data source: BioProject PRJNA1176350
//   TK12 SRR31089074, TK13 SRR31089073, TK14 SRR31089072 (+ additional
//   triplicate accessions -- check the BioProject page and add rows to CSV)
//
// Steps:
//   0.  SRA_DOWNLOAD         -- prefetch + fasterq-dump per SRR accession
//   0b. CAT_FASTQ            -- concatenate per-sample SRR FASTQs (triplicates)
//   1.  STAR_GENOMEGENERATE  -- build STAR index (once)
//   2.  TRIM_GALORE          -- adapter trim + quality/length filter
//   3.  STAR_ALIGN           -- 2-pass alignment + chimeric read output for Arriba
//   4.  SAMTOOLS_SORT_INDEX  -- coordinate-sort + index
//   5.  MARKDUPLICATES       -- mark PCR/optical duplicates
//   6.  REMAP_MAPQ           -- remap STAR MAPQ 255 -> 60 (safety net)
//   7.  SPLIT_N_CIGAR        -- split reads at N-CIGAR splice junctions (REQUIRED for GATK RNA-seq)
//   8.  BQSR_RECAL           -- base quality score recalibration table
//   9.  BQSR_APPLY           -- apply recalibrated qualities
//  10.  HAPLOTYPE_CALLER     -- per-sample GVCF (replaces Mutect2: catches germline + somatic)
//  11.  GENOTYPE_GVCFS       -- GVCF -> VCF per sample
//  12.  VARIANT_FILTRATION   -- RNA-seq hard filters (QD, FS, MQ, ReadPosRankSum)
//  13.  SELECT_PASS          -- retain only PASS variants
//  14.  VEP_ANNOTATE         -- Ensembl VEP: consequence + HGVS + protein change
//  15.  PYPGATK_FASTA        -- variant VCF -> per-sample variant protein FASTA
//  16.  ARRIBA               -- gene fusion calling from chimeric reads (replaces DELLY)
//  17.  TK_PROGRESSION_SUB   -- pseudo-paired subtraction: remove TK12 germline/early
//                               variants from TK13/TK14 to isolate progression-specific calls
//  18.  PROGRESSION_FASTA    -- variant FASTA for progression-only variants
//  19.  MERGE_PER_GROUP      -- merge per-sample VCFs into per-group VCF
//
// Key design decisions vs previous version:
//   - HaplotypeCaller (germline + somatic diploid) replaces Mutect2 tumor-only.
//     For variant proteomics we want ALL protein-coding variants per sample,
//     not just somatic. Mutect2 without a PoN also hallucinate somatic calls
//     at germline sites and suppress heterozygous germline SNPs as contamination.
//   - SplitNCigarReads is the critical missing GATK RNA-seq step: without it
//     GATK generates false positives at every splice junction.
//   - Arriba replaces DELLY. DELLY was designed for WGS/WES; it interprets
//     normal splicing as deletions on RNA-seq BAMs. Arriba uses STAR chimeric
//     output specifically designed for RNA fusion detection.
//   - The TK12 pseudo-paired progression analysis exploits the unique
//     longitudinal design: TK12 (earliest timepoint) acts as per-patient
//     germline/early-clonal reference. Variants appearing in TK13/TK14 but
//     absent in TK12 are progression-specific and most biologically interesting
//     in the context of ATX-101 sensitivity.
// ---------------------------------------------------------------------------


// ---------------------------------------------------------------------------
// 1. Build STAR genome index (runs once)
//    Added: --sjdbOverhang for PE150 data
//    Chimeric params live in STAR_ALIGN, not here.
// ---------------------------------------------------------------------------

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


// ---------------------------------------------------------------------------
// 2. Adapter trimming
// ---------------------------------------------------------------------------

process TRIM_GALORE {
    tag "${meta.sample}"
    label 'process_medium'
    container 'wave.seqera.io/wt/a56a73714737/wave/build:trim-galore--f48a39e5cc6e0d0d'

    publishDir "${params.outdir}/trimmed", mode: 'copy', pattern: '*_trimming_report.txt'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)

    output:
    tuple val(meta), path("${meta.sample}_R1_trimmed.fq.gz"), path("${meta.sample}_R2_trimmed.fq.gz"), emit: reads
    path "*_trimming_report.txt", emit: reports

    script:
    """
    trim_galore \\
        --paired \\
        --quality 20 \\
        --length 36 \\
        --cores ${task.cpus} \\
        --gzip \\
        --basename ${meta.sample} \\
        ${fastq_1} ${fastq_2}

    mv ${meta.sample}_R1_val_1.fq.gz ${meta.sample}_R1_trimmed.fq.gz
    mv ${meta.sample}_R2_val_2.fq.gz ${meta.sample}_R2_trimmed.fq.gz
    """
}


// ---------------------------------------------------------------------------
// 3. STAR 2-pass alignment
//    CHANGED vs previous:
//    - outFilterMultimapNmax 1 -> 50: Arriba needs multimapper chimeric reads.
//      GATK ignores multimappers via MAPQ filter (unique reads = MAPQ 60).
//    - Added chimeric output params required by Arriba:
//        --chimOutType WithinBAM SeparateSAMold  (Chimeric.out.sam for Arriba)
//        --chimSegmentMin 10
//        --chimJunctionOverhangMin 10
//        --chimScoreMin 1, --chimScoreDropMax 20
//        --chimMultimapNmax 50
//        --alignSJstitchMismatchNmax 5 -1 5 5
//        --chimSegmentReadGapMax 3
//        --peOverlapNbasesMin 12
//    - outSAMmapqUnique 60 retained: unique mappers still get MAPQ 60 for GATK.
// ---------------------------------------------------------------------------

process STAR_ALIGN {
    tag "${meta.sample}"
    label 'process_high'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'

    input:
    tuple val(meta), path(fastq_1), path(fastq_2)
    path star_index

    output:
    tuple val(meta), path("${meta.sample}.Aligned.sortedByCoord.out.bam"), emit: bam
    tuple val(meta), path("${meta.sample}.Chimeric.out.sam"),              emit: chimeric

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
        --outFilterMultimapNmax 50 \\
        --outFilterMismatchNmax 6 \\
        --outFilterMismatchNoverReadLmax 0.04 \\
        --chimOutType WithinBAM SeparateSAMold \\
        --chimSegmentMin 10 \\
        --chimJunctionOverhangMin 10 \\
        --chimScoreMin 1 \\
        --chimScoreDropMax 20 \\
        --chimNonchimScoreDropMin 10 \\
        --chimMultimapNmax 50 \\
        --chimMultimapScoreRange 3 \\
        --alignSJstitchMismatchNmax 5 -1 5 5 \\
        --chimSegmentReadGapMax 3 \\
        --peOverlapNbasesMin 12 \\
        --peOverlapMMp 0.1
    """
}


// ---------------------------------------------------------------------------
// 4. Sort + index BAM
// ---------------------------------------------------------------------------

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


// ---------------------------------------------------------------------------
// 5. Mark duplicates
// ---------------------------------------------------------------------------

process MARKDUPLICATES {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    publishDir "${params.outdir}/bam", mode: 'copy', pattern: '*.markdup.bam*'

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("${meta.sample}.markdup.bam"), path("${meta.sample}.markdup.bam.bai"), emit: bam
    path "${meta.sample}.markdup.metrics.txt", emit: metrics

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


// ---------------------------------------------------------------------------
// 6. Remap MAPQ 255 -> 60 (safety net)
// ---------------------------------------------------------------------------

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


// ---------------------------------------------------------------------------
// 7. SplitNCigarReads  [NEW - was missing, critical for GATK RNA-seq]
//
//    RNA-seq alignments contain N-CIGAR operations at splice junctions.
//    GATK HaplotypeCaller/Mutect2 was not designed to handle these and will
//    generate false variant calls at every splice site without this step.
//    SplitNCigarReads splits each read spanning a junction into separate
//    segments, discarding the intronic overhang.
//
//    This must run BEFORE BQSR and variant calling.
//    Reference: GATK best practices for RNA-seq variant discovery.
// ---------------------------------------------------------------------------

process SPLIT_N_CIGAR {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.split.bam"), path("${meta.sample}.split.bam.bai")

    script:
    """
    gatk SplitNCigarReads \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.sample}.split.bam \\
        --create-output-bam-index true
    """
}


// ---------------------------------------------------------------------------
// Reference indexing (runs once)
// ---------------------------------------------------------------------------

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


// ---------------------------------------------------------------------------
// 8-9. BQSR (unchanged logic; now runs on split-N BAMs)
// ---------------------------------------------------------------------------

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


// ---------------------------------------------------------------------------
// 10. HaplotypeCaller  [REPLACES Mutect2]
//
//     Why HaplotypeCaller instead of Mutect2 for variant proteomics:
//     - Goal: a complete per-sample protein-coding variant catalogue
//       (germline SNPs + somatic mutations) to build a variant FASTA.
//       Mutect2 suppresses heterozygous germline variants at AF~0.5 as
//       apparent "normal contamination" and without a PoN calls many
//       germline sites as somatic. Both behaviours lose real protein variants.
//     - HaplotypeCaller calls diploid genotypes faithfully at all covered sites.
//     - GVCF mode (-ERC GVCF) enables GenotypeGVCFs for future joint
//       genotyping across all TK samples if needed.
//     - --dont-use-soft-clipped-bases: required for RNA-seq; soft-clipped
//       bases at junctions are alignment artefacts, not real variants.
//     - -stand-call-conf 20: relaxed from default 30 to improve sensitivity
//       at lower-covered coding regions in 20M read RNA-seq data.
// ---------------------------------------------------------------------------

process HAPLOTYPE_CALLER {
    tag "${meta.sample}"
    label 'process_high'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    publishDir "${params.outdir}/gvcf", mode: 'copy', pattern: '*.g.vcf.gz*'

    input:
    tuple val(meta), path(bam), path(bai)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.g.vcf.gz"), path("${meta.sample}.g.vcf.gz.tbi")

    script:
    """
    gatk HaplotypeCaller \\
        -R ${fasta} \\
        -I ${bam} \\
        -O ${meta.sample}.g.vcf.gz \\
        -ERC GVCF \\
        --dont-use-soft-clipped-bases \\
        -stand-call-conf 20 \\
        --native-pair-hmm-threads ${task.cpus}
    """
}


// ---------------------------------------------------------------------------
// 11. GenotypeGVCFs -- convert per-sample GVCF to called VCF
// ---------------------------------------------------------------------------

process GENOTYPE_GVCFS {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    input:
    tuple val(meta), path(gvcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.raw.vcf.gz"), path("${meta.sample}.raw.vcf.gz.tbi")

    script:
    """
    gatk GenotypeGVCFs \\
        -R ${fasta} \\
        -V ${gvcf} \\
        -O ${meta.sample}.raw.vcf.gz
    """
}


// ---------------------------------------------------------------------------
// 12. VariantFiltration -- RNA-seq hard filters
//
//     GATK recommends hard filters instead of VQSR for RNA-seq because
//     the annotation distributions differ from WGS/WES and VQSR training
//     sets are not appropriate.
//     Cluster/window filter removes tightly clustered SNPs (common artefact
//     near splice sites even after SplitNCigarReads).
//     These are GATK-recommended thresholds for RNA-seq:
//       QD < 2.0        : low variant quality relative to depth (noisy site)
//       FS > 30.0       : high strand bias (SNP threshold; WGS uses 60)
//       MQ < 40.0       : low mapping quality of supporting reads
//       MQRankSum < -12.5 : alt reads map worse than ref reads
//       ReadPosRankSum < -8.0 : alt allele near read ends (soft-clip artefact)
// ---------------------------------------------------------------------------

process VARIANT_FILTRATION {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.filtered.vcf.gz"), path("${meta.sample}.filtered.vcf.gz.tbi")

    script:
    """
    gatk VariantFiltration \\
        -R ${fasta} \\
        -V ${vcf} \\
        --window 35 \\
        --cluster 3 \\
        --filter-expression "QD < 2.0"          --filter-name "QD2" \\
        --filter-expression "FS > 30.0"          --filter-name "FS30" \\
        --filter-expression "MQ < 40.0"          --filter-name "MQ40" \\
        --filter-expression "MQRankSum < -12.5"  --filter-name "MQRankSum-12.5" \\
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPos-8" \\
        -O ${meta.sample}.filtered.vcf.gz
    """
}


// ---------------------------------------------------------------------------
// 13. SelectVariants -- retain only PASS records for downstream steps
// ---------------------------------------------------------------------------

process SELECT_PASS {
    tag "${meta.sample}"
    label 'process_low'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'

    publishDir "${params.outdir}/vcf_pass/${meta.group}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)

    output:
    tuple val(meta), path("${meta.sample}.pass.vcf.gz"), path("${meta.sample}.pass.vcf.gz.tbi")

    script:
    """
    gatk SelectVariants \\
        -R ${fasta} \\
        -V ${vcf} \\
        --exclude-filtered \\
        -O ${meta.sample}.pass.vcf.gz
    """
}


// ---------------------------------------------------------------------------
// 14. VEP annotation  [NEW]
//
//     Annotates each variant with:
//       - Predicted consequence (missense, frameshift, stop-gained, etc.)
//       - HGVS protein notation (p.Arg123His)
//       - Gene symbol, Ensembl transcript, canonical transcript flag
//       - Biotype filter: protein_coding only (--biotype)
//     The CSQ field produced by VEP is what pypgatk reads in the next step.
//
//     Requirements:
//       params.vep_cache     -- path to pre-downloaded VEP cache directory
//                              (homo_sapiens GRCh38; ~15 GB; download once with:
//                               vep_install -a cf -s homo_sapiens -y GRCh38
//                                           --CACHEDIR /path/to/vep_cache)
//       params.vep_assembly  -- "GRCh38" (default)
// ---------------------------------------------------------------------------

process VEP_ANNOTATE {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0'

    publishDir "${params.outdir}/vep/${meta.group}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    tuple path(fasta), path(fai), path(dict)
    path vep_cache

    output:
    tuple val(meta), path("${meta.sample}.vep.vcf.gz"), path("${meta.sample}.vep.vcf.gz.tbi")

    script:
    def assembly = params.vep_assembly ?: "GRCh38"
    """
    vep \\
        --input_file ${vcf} \\
        --output_file ${meta.sample}.vep.vcf \\
        --format vcf \\
        --vcf \\
        --cache \\
        --offline \\
        --dir_cache ${vep_cache} \\
        --species homo_sapiens \\
        --assembly ${assembly} \\
        --fasta ${fasta} \\
        --pick \\
        --canonical \\
        --protein \\
        --symbol \\
        --numbers \\
        --biotype \\
        --total_length \\
        --hgvs \\
        --sift b \\
        --polyphen b \\
        --fork ${task.cpus}

    bgzip ${meta.sample}.vep.vcf
    tabix -p vcf ${meta.sample}.vep.vcf.gz
    """
}


// ---------------------------------------------------------------------------
// 15. pypgatk: variant VCF -> variant protein FASTA  [NEW]
//
//     pypgatk reads the VEP-annotated VCF (CSQ field) and generates a FASTA
//     of predicted mutant protein sequences for each non-synonymous variant.
//     These are appended to the reference Uniprot/Swissprot FASTA for
//     database search in MaxQuant/MSFragger.
//
//     Consequence filter retains variants that change protein sequence:
//       missense_variant, frameshift_variant, stop_gained, stop_lost,
//       start_lost, splice_donor_variant, splice_acceptor_variant,
//       inframe_insertion, inframe_deletion
//
//     Output: ${sample}.variant_proteins.fasta
//       To be concatenated with canonical proteome FASTA for MS/MS search.
//
//     Requirements:
//       params.reference_proteome -- Uniprot/Swissprot human FASTA
//                                    (used by pypgatk as reference translation base)
// ---------------------------------------------------------------------------

process PYPGATK_FASTA {
    tag "${meta.sample}"
    label 'process_medium'
    container 'quay.io/biocontainers/pypgatk:0.0.25--pyhdfd78af_0'

    publishDir "${params.outdir}/variant_fasta/${meta.group}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path gtf
    path reference_proteome

    output:
    tuple val(meta), path("${meta.sample}.variant_proteins.fasta")

    script:
    """
    pypgatk vcf-to-proteindb \\
        --input-vcf ${vcf} \\
        --gene-annotations-gtf ${gtf} \\
        --protein-db-fasta ${reference_proteome} \\
        --af-field AF \\
        --annotation-field-name CSQ \\
        --consequence-filter missense_variant,frameshift_variant,stop_gained,stop_lost,start_lost,splice_donor_variant,splice_acceptor_variant,inframe_insertion,inframe_deletion \\
        --output-proteindb ${meta.sample}.variant_proteins.fasta

    # Append sample tag to FASTA headers to track variant origin in MS results
    sed -i "s/^>/>${meta.sample}|/" ${meta.sample}.variant_proteins.fasta
    """
}


// ---------------------------------------------------------------------------
// 16. Arriba gene fusion calling  [REPLACES DELLY]
//
//     Why Arriba instead of DELLY for RNA-seq:
//     - DELLY is a WGS/WES SV caller. On RNA-seq BAMs it interprets normal
//       constitutive splicing as deletions and generates thousands of false SVs.
//     - Arriba is purpose-built for RNA-seq fusion detection. It uses the
//       chimeric read output from STAR (Chimeric.out.sam) alongside the main
//       BAM to identify genuine inter-gene fusions with low false-positive rate.
//     - In MM, gene fusions involving IGH, FGFR3, CCND1, MAF are clinically
//       relevant and directly testable in the MS data.
//
//     Requirements:
//       params.arriba_blacklist      -- ${arriba_dir}/database/blacklist_hg38_GRCh38_v2.4.0.tsv.gz
//       params.arriba_known_fusions  -- ${arriba_dir}/database/known_fusions_hg38_GRCh38_v2.4.0.tsv.gz
//       params.arriba_protein_domains -- ${arriba_dir}/database/protein_domains_hg38_GRCh38_v2.4.0.gff3
//       (all bundled with Arriba release; download from github.com/suhrig/arriba/releases)
// ---------------------------------------------------------------------------

process ARRIBA {
    tag "${meta.sample}"
    label 'process_high'
    container 'quay.io/biocontainers/arriba:2.4.0--h0033a41_2'

    publishDir "${params.outdir}/fusions/${meta.group}", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai), path(chimeric_sam)
    tuple path(fasta), path(fai), path(dict)
    path gtf
    path blacklist
    path known_fusions
    path protein_domains

    output:
    tuple val(meta), path("${meta.sample}.fusions.tsv"),         emit: fusions
    tuple val(meta), path("${meta.sample}.fusions.discarded.tsv"), emit: discarded

    script:
    """
    arriba \\
        -x ${bam} \\
        -c ${chimeric_sam} \\
        -a ${fasta} \\
        -g ${gtf} \\
        -b ${blacklist} \\
        -k ${known_fusions} \\
        -p ${protein_domains} \\
        -o ${meta.sample}.fusions.tsv \\
        -O ${meta.sample}.fusions.discarded.tsv \\
        -T \\
        -P
    """
}


// ---------------------------------------------------------------------------
// 17. TK progression subtraction  [NEW]
//
//     Design rationale:
//     TK12, TK13, TK14 are from the same patient at increasing timepoints.
//     TK12 acts as pseudo-baseline: variants present in TK12 are either
//     germline or early clonal events present before the period of interest.
//     bcftools isec -C returns records private to file A (not in file B):
//       private to progression sample = absent in TK12 = acquired during
//       disease progression between TK12 and TK13/TK14.
//
//     These progression-specific variants are the most biologically relevant
//     for correlating with ATX-101 hypersensitivity observed in TK13/TK14.
//
//     Input:  progression sample PASS VCF + its TK-matched baseline PASS VCF
//     Output: VCF containing only variants absent from the baseline
// ---------------------------------------------------------------------------

process TK_PROGRESSION_SUB {
    tag "${meta.sample}_vs_${meta.tk}_baseline"
    label 'process_low'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    publishDir "${params.outdir}/progression/${meta.group}", mode: 'copy'

    input:
    tuple val(meta), path(prog_vcf), path(prog_tbi), path(base_vcf), path(base_tbi)

    output:
    tuple val(meta), path("${meta.sample}.progression_only.vcf.gz"), path("${meta.sample}.progression_only.vcf.gz.tbi")

    script:
    """
    # -C: output only records private to the first file (not present in baseline)
    # -w1: use the first file's coordinates (progression sample)
    # Intersection is by chromosome + position + ref + alt (exact match)
    bcftools isec \\
        -C \\
        -O z \\
        -w 1 \\
        -o ${meta.sample}.progression_only.vcf.gz \\
        ${prog_vcf} ${base_vcf}

    bcftools index --tbi ${meta.sample}.progression_only.vcf.gz

    # Report how many progression-specific variants were identified
    echo "Progression-specific variants in ${meta.sample} vs TK baseline:"
    bcftools stats ${meta.sample}.progression_only.vcf.gz | grep "^SN"
    """
}


// ---------------------------------------------------------------------------
// 18. FASTA for progression-only variants  [NEW]
//     Same pypgatk logic as step 15, applied to progression-specific VCF.
//     These FASTAs are the most informative for the ATX-101 proteomics search
//     because they contain only variants acquired during progression to the
//     ATX-101-sensitive state.
// ---------------------------------------------------------------------------

process PROGRESSION_FASTA {
    tag "${meta.sample}_progression"
    label 'process_medium'
    container 'quay.io/biocontainers/pypgatk:0.0.25--pyhdfd78af_0'

    publishDir "${params.outdir}/progression_fasta/${meta.group}", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path gtf
    path reference_proteome

    output:
    tuple val(meta), path("${meta.sample}.progression_proteins.fasta")

    script:
    """
    # VEP annotation on progression-specific VCF first (inline, lightweight)
    pypgatk vcf-to-proteindb \\
        --input-vcf ${vcf} \\
        --gene-annotations-gtf ${gtf} \\
        --protein-db-fasta ${reference_proteome} \\
        --af-field AF \\
        --annotation-field-name CSQ \\
        --consequence-filter missense_variant,frameshift_variant,stop_gained,stop_lost,start_lost,splice_donor_variant,splice_acceptor_variant,inframe_insertion,inframe_deletion \\
        --output-proteindb ${meta.sample}.progression_proteins.fasta

    sed -i "s/^>/>${meta.sample}|PROGRESSION|/" ${meta.sample}.progression_proteins.fasta
    """
}


// ---------------------------------------------------------------------------
// 19. Merge per-sample PASS VCFs into per-group VCF
//     Useful for group-level annotation and downstream cohort analysis.
// ---------------------------------------------------------------------------

process MERGE_PER_GROUP {
    tag "${group}"
    label 'process_medium'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'

    publishDir "${params.outdir}/vcf_merged", mode: 'copy'

    input:
    tuple val(group), path(vcfs), path(tbis)

    output:
    tuple val(group), path("${group}.merged.vcf.gz"), path("${group}.merged.vcf.gz.tbi")

    script:
    def vcf_list = (vcfs instanceof List ? vcfs : [vcfs]).collect { v -> v.toString() }.join(' ')
    """
    n_vcfs=\$(echo "${vcf_list}" | wc -w)
    if [ "\${n_vcfs}" -eq 1 ]; then
        cp ${vcf_list} ${group}.merged.vcf.gz
        bcftools index --tbi ${group}.merged.vcf.gz
    else
        bcftools merge \\
            --force-samples \\
            --merge none \\
            -O z \\
            -o ${group}.merged.vcf.gz \\
            ${vcf_list}
        bcftools index --tbi ${group}.merged.vcf.gz
    fi
    """
}


// ===========================================================================
// WORKFLOW
// ===========================================================================

workflow {

    // -------------------------------------------------------------------------
    // Parse samples.csv
    // Columns: sample, fastq_1, fastq_2, strandedness, TK, Group, baseline
    //   TK       : patient ID linking longitudinal samples (e.g. "TK_patient1")
    //   Group    : sensitive / resistant
    //   baseline : "true" = earliest timepoint per patient (TK12-equivalent)
    //              "false" = later timepoints (TK13, TK14)
    //              ""      = single-timepoint samples (no progression analysis)
    // -------------------------------------------------------------------------
    def ch_reads = channel.fromPath(params.samplesheet)
        | splitCsv(header: true)
        | map { row ->
            def meta = [
                sample   : row.sample,
                group    : row.Group,
                tk       : row.TK,
                baseline : row.baseline ?: ""
            ]
            tuple(meta, file(row.fastq_1), file(row.fastq_2))
        }

    // -------------------------------------------------------------------------
    // Reference channels
    // -------------------------------------------------------------------------
    def ch_fasta             = channel.value(file(params.genome))
    def ch_gtf               = channel.value(file(params.gtf))
    def ch_ref_proteome      = channel.value(file(params.reference_proteome))
    def ch_vep_cache         = channel.value(file(params.vep_cache))
    def ch_arriba_blacklist  = channel.value(file(params.arriba_blacklist))
    def ch_arriba_known      = channel.value(file(params.arriba_known_fusions))
    def ch_arriba_domains    = channel.value(file(params.arriba_protein_domains))

    // -------------------------------------------------------------------------
    // Build STAR index once
    // -------------------------------------------------------------------------
    def ch_star_index = params.star_index
        ? channel.value(file(params.star_index))
        : STAR_GENOMEGENERATE(ch_fasta, ch_gtf).collect().map { dirs -> dirs[0] }

    // -------------------------------------------------------------------------
    // Reference genome indexing (runs once)
    // -------------------------------------------------------------------------
    def ch_fai = SAMTOOLS_FAIDX(ch_fasta)
    def ch_ref = GATK_DICT(ch_fai)

    // -------------------------------------------------------------------------
    // Read processing: trim -> align -> sort -> markdup -> remapq -> split-N
    // -------------------------------------------------------------------------
    def ch_trimmed  = params.skip_trimming ? ch_reads : TRIM_GALORE(ch_reads).reads
    def ch_aligned  = STAR_ALIGN(ch_trimmed, ch_star_index)
    def ch_sorted   = SAMTOOLS_SORT_INDEX(ch_aligned.bam)
    def ch_markdup  = MARKDUPLICATES(ch_sorted).bam
    def ch_remapped = REMAP_MAPQ(ch_markdup)
    def ch_split    = SPLIT_N_CIGAR(ch_remapped, ch_ref)

    // -------------------------------------------------------------------------
    // BQSR (optional; skip if no dbSNP available)
    // -------------------------------------------------------------------------
    def ch_for_calling
    if (!params.skip_bqsr && params.dbsnp) {
        def ch_dbsnp     = channel.value(file(params.dbsnp))
        def ch_dbsnp_tbi = channel.value(file("${params.dbsnp}.tbi"))
        def ch_recal     = BQSR_RECAL(ch_split, ch_ref, ch_dbsnp, ch_dbsnp_tbi)
        ch_for_calling   = BQSR_APPLY(ch_recal, ch_ref)
    } else {
        ch_for_calling = ch_split
    }

    // -------------------------------------------------------------------------
    // Variant calling: HaplotypeCaller -> GenotypeGVCFs -> filter -> PASS
    // -------------------------------------------------------------------------
    def ch_gvcf     = HAPLOTYPE_CALLER(ch_for_calling, ch_ref)
    def ch_raw_vcf  = GENOTYPE_GVCFS(ch_gvcf, ch_ref)
    def ch_filtered = VARIANT_FILTRATION(ch_raw_vcf, ch_ref)
    def ch_pass     = SELECT_PASS(ch_filtered, ch_ref)

    // -------------------------------------------------------------------------
    // VEP annotation -> per-sample variant protein FASTA
    // -------------------------------------------------------------------------
    def ch_vep      = VEP_ANNOTATE(ch_pass, ch_ref, ch_vep_cache)
    PYPGATK_FASTA(ch_vep, ch_gtf, ch_ref_proteome)

    // -------------------------------------------------------------------------
    // Merge PASS VCFs per group (sensitive / resistant)
    // -------------------------------------------------------------------------
    ch_pass
        | map { meta, vcf, tbi -> tuple(meta.group, vcf, tbi) }
        | groupTuple(by: 0)
        | MERGE_PER_GROUP

    // -------------------------------------------------------------------------
    // Fusion calling with Arriba
    // Join BAM channel with chimeric SAM channel by sample name
    // -------------------------------------------------------------------------
    def ch_chimeric = STAR_ALIGN.out.chimeric

    def ch_arriba_input = ch_for_calling
        | map { meta, bam, bai -> tuple(meta.sample, meta, bam, bai) }
        | join(
            ch_chimeric | map { meta, sam -> tuple(meta.sample, sam) },
            by: 0
          )
        | map { sample_id, meta, bam, bai, chimeric_sam ->
            tuple(meta, bam, bai, chimeric_sam)
          }

    ARRIBA(ch_arriba_input, ch_ref, ch_gtf, ch_arriba_blacklist, ch_arriba_known, ch_arriba_domains)

    // -------------------------------------------------------------------------
    // TK progression analysis: pseudo-paired subtraction
    //
    // Logic:
    //   1. Split ch_pass into baseline channel (meta.baseline == "true")
    //      and progression channel (meta.baseline == "false").
    //   2. Key both by meta.tk (patient ID) to pair correctly.
    //   3. combine() pairs each progression sample with its patient's baseline.
    //   4. TK_PROGRESSION_SUB subtracts baseline variants.
    //   5. VEP + PROGRESSION_FASTA on the subtracted VCF.
    //
    // Samples with meta.baseline == "" (single timepoints: TK9, TK10, TK16, TK18)
    // are handled by standard PYPGATK_FASTA above; they have no baseline to subtract.
    // -------------------------------------------------------------------------
    def ch_pass_split = ch_pass
        | branch { meta, vcf, tbi ->
            baseline    : meta.baseline == "true"
            progression : meta.baseline == "false"
            single      : true        // single-timepoint, no subtraction needed
        }

    // Key baseline VCFs by TK patient group
    def ch_baselines = ch_pass_split.baseline
        | map { meta, vcf, tbi -> tuple(meta.tk, vcf, tbi) }

    // Pair each progression sample with its patient's baseline
    def ch_progression_input = ch_pass_split.progression
        | map { meta, vcf, tbi -> tuple(meta.tk, meta, vcf, tbi) }
        | combine(ch_baselines, by: 0)
        | map { tk, meta, prog_vcf, prog_tbi, base_vcf, base_tbi ->
            tuple(meta, prog_vcf, prog_tbi, base_vcf, base_tbi)
        }

    // Subtract baseline, then annotate and generate FASTA for progression variants
    def ch_prog_vcf = TK_PROGRESSION_SUB(ch_progression_input)

    // VEP annotate progression-specific VCF before FASTA generation
    def ch_prog_vep = VEP_ANNOTATE(ch_prog_vcf, ch_ref, ch_vep_cache)
    PROGRESSION_FASTA(ch_prog_vep, ch_gtf, ch_ref_proteome)
}
