#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplesheet = "${projectDir}/samples.csv"
params.outdir = "${projectDir}/results"
params.read_length = 150
params.sra_dir = null
params.fusion_flank_aa = 50
params.splice_min_coverage = 2.5
params.splice_min_junction_reads = 3
params.splice_min_isoform_fraction = 0.05
params.splice_min_protein_aa = 60
params.splice_class_codes = 'j,u'
params.run_proteogenomic_validation = false
params.maxquant_txt = "${projectDir}/txtMQMBR"
params.maxquant_mqpar = null
params.maxquant_canonical_fasta = null
params.maxquant_contaminants = null
params.maxquant_raw_map = null
params.ensembl_pep = null
params.read_validation_padding = 150
params.rna_variant_min_depth = 10
params.rna_variant_min_alt_reads = 3
params.rna_variant_min_alt_fraction = 0.05
params.rna_fusion_min_split_reads = 1
params.rna_fusion_min_total_support = 2
params.haplotype_scatter_count = 24
params.hc_calling_confidence = 20
params.hc_dont_use_soft_clipped_bases = true
params.hc_pcr_indel_model = 'CONSERVATIVE'
params.snp_filter_qd = 2.0
params.snp_filter_fs = 60.0
params.snp_filter_sor = 3.0
params.snp_filter_mq = 40.0
params.snp_filter_mq_rank_sum = -12.5
params.snp_filter_read_pos_rank_sum = -8.0
params.indel_filter_qd = 2.0
params.indel_filter_fs = 200.0
params.indel_filter_sor = 10.0
params.indel_filter_read_pos_rank_sum = -20.0
params.finding_review_mapq = 20
params.finding_review_baseq = 20
params.finding_review_reference_reads = 20
params.finding_classes = 'rna_variant,progression_variant,fusion,splice_junction'
params.finding_primary_class_order = 'rna_variant,progression_variant'
params.finding_priority_mode = 'all'
params.finding_priority_genes = ''
params.finding_priority_impacts = ''
params.finding_priority_consequences = ''
params.generate_priority_igv_reports = true
params.igv_report_classes = 'rna_variant,progression_variant,fusion,splice_junction'
params.igv_report_max_reads = 100
params.igv_report_max_file_size_mb = 64
params.igv_report_gene_filter = ''
params.igv_report_sample_filter = ''
params.igv_report_timeout_seconds = 600
params.igv_report_title_prefix = 'PGTK finding'
params.run_external_vcf_comparison = false
params.external_vcf_dir = "${projectDir}/sarek"
params.external_vcf_suffix = '.haplotypecaller.vcf.gz'
params.reference_downloads = null
params.container_cache = null
params.pysam_image = null
params.go_obo = "${projectDir}/reference_downloads/go-basic.obo"
params.go_gaf = "${projectDir}/reference_downloads/goa_human.gaf.gz"
params.go_min_size = 10
params.go_max_size = 500
params.go_fdr_threshold = 0.1
params.go_namespaces = 'all'
params.run_expression_go = true
params.gene_count_feature_type = 'exon'
params.gene_count_id_attribute = 'gene_id'
params.gene_count_symbol_attribute = 'gene_name'
params.gene_count_biotypes = 'all'
params.gene_count_strandedness = 0
params.gene_count_min_mapq = 10
params.gene_count_min_overlap = 1
params.gene_count_count_read_pairs = true
params.gene_count_require_both_ends = false
params.gene_count_exclude_chimeric = true
params.gene_count_primary_only = true
params.gene_count_allow_multi_overlap = false
params.gene_count_count_multimapping = false
params.expression_pseudocount = 0.5
params.expression_cpm_threshold = 1.0
params.expression_tpm_threshold = 0.0
params.expression_rank_metric = 'log2_tpm_fold_change'
params.expression_rank_min_nonzero_scores = 1
params.go_expression_background = 'genome'
params.go_variant_background = 'genome'
params.go_variant_biotypes = 'protein_coding'


def strictBooleanParam(value, String name) {
    if (value instanceof Boolean) return value
    if (value == null) error "${name} requires true or false"
    def normalized = value.toString().trim().toLowerCase()
    if (normalized in ['true','1','yes','y','on']) return true
    if (normalized in ['false','0','no','n','off']) return false
    error "${name} must be true or false; received '${value}'"
}

def booleanText(value, String name) {
    return strictBooleanParam(value, name) ? 'true' : 'false'
}


def resolveExternalVcf(String directory, String srr, String suffix) {
    def root = new File(directory)
    if (!root.isDirectory()) error "External caller folder not found: ${directory}"
    def expected = "${srr}${suffix}"
    def preferred = [
        new File(root, expected),
        new File(new File(root, srr), expected)
    ].findAll { it.isFile() }.unique { it.canonicalPath }
    if (preferred.size() == 1) {
        if (preferred[0].length() == 0L) error "External caller VCF is empty: ${preferred[0]}"
        return file(preferred[0].toString(), checkIfExists:true)
    }
    if (preferred.size() > 1) error "Ambiguous external caller VCF for ${srr}: ${preferred*.toString().join(', ')}"
    def matches = []
    root.eachFileRecurse { candidate ->
        if (candidate.isFile() && candidate.name == expected) matches << candidate
    }
    matches = matches.unique { it.canonicalPath }
    if (matches.size() == 0) error "External caller VCF not found. Expected ${root}/${expected} or ${root}/${srr}/${expected}"
    if (matches.size() > 1) error "Ambiguous external caller VCF for ${srr}; found ${matches.size()} matches below ${directory}: ${matches*.toString().join(', ')}"
    if (matches[0].length() == 0L) error "External caller VCF is empty: ${matches[0]}"
    return file(matches[0].toString(), checkIfExists:true)
}

def resolveContaminants(String configured, String maxquantTxt) {
    if (configured) return file(configured, checkIfExists:true)
    def roots = [new File(maxquantTxt), new File(maxquantTxt).parentFile, new File(projectDir.toString())].findAll { it?.isDirectory() }
    def matches=[]
    roots.each { root -> root.eachFileRecurse { f -> if (f.isFile() && f.name.toLowerCase() in ['contaminants.fasta','contaminants.fa']) matches << f } }
    def homeCandidate = new File(System.getProperty('user.home'), 'scripts/MaxQuant_v2.8.1.0/bin/conf/contaminants.fasta')
    if (homeCandidate.isFile()) matches << homeCandidate
    matches=matches.unique { it.canonicalPath }
    if (matches.size()!=1) error "Set --maxquant_contaminants or place exactly one contaminants.fasta under txtMQMBR, its parent, or projectDir; found ${matches.size()}"
    return file(matches[0].toString(), checkIfExists:true)
}

def resolveMaxQuantMqpar(maxquant_txt, override_path) {
    if (override_path) return file(override_path, checkIfExists:true)
    def txt = new File(maxquant_txt.toString())
    def candidates = [new File(txt, 'mqpar.xml'), new File(txt.parentFile, 'mqpar.xml')]
    def found = candidates.find { it?.isFile() && it.length() > 0 }
    if (!found) error "Cannot find mqpar.xml inside --maxquant_txt or its parent; use --maxquant_mqpar to override"
    return file(found.toString(), checkIfExists:true)
}

def maxQuantFastaPaths(mqpar_file) {
    def matcher = mqpar_file.text =~ /(?is)<fastaFilePath>\s*([^<]+?)\s*<\/fastaFilePath>/
    def paths = matcher.collect { matched ->
        matched[1].trim()
    }
    if (!paths) error "mqpar.xml contains no fastaFilePath entries"
    return paths
}

def fastaBasename(path_value) {
    return path_value.toString().replace('\\', '/').tokenize('/').last()
}

def resolveExistingFasta(configured_path, mqpar_file, maxquant_txt) {
    def basename = fastaBasename(configured_path)
    def txt = new File(maxquant_txt.toString())
    def home = System.getenv('HOME') ?: ''
    def candidates = [
        new File(configured_path.toString()),
        new File(mqpar_file.parent.toString(), basename),
        new File(txt, basename),
        new File(txt.parentFile, basename),
        new File(projectDir.toString(), basename),
        home ? new File(home, "FastaDB/${basename}") : null,
    ].findAll { it != null }
    def found = candidates.find { it.isFile() && it.length() > 0 }
    if (!found) error "FASTA '${basename}' from mqpar.xml was not found in the recorded path, beside mqpar.xml, in MQTXT, its parent, projectDir, or ~/FastaDB"
    return file(found.toString(), checkIfExists:true)
}

def resolveMaxQuantCanonicalFastas(mqpar_file, maxquant_txt, override_value, sample_fasta_basenames) {
    if (override_value) {
        def values = override_value instanceof List ? override_value : override_value.toString().split(',').collect { it.trim() }.findAll { it }
        return values.collect { file(it, checkIfExists:true) }
    }
    def sample_names = sample_fasta_basenames.collect { it.toLowerCase() } as Set
    def canonical_entries = maxQuantFastaPaths(mqpar_file).findAll { configured ->
        def basename = fastaBasename(configured)
        !sample_names.contains(basename.toLowerCase()) && !basename.toLowerCase().contains('contaminant')
    }
    if (!canonical_entries) error "No canonical FASTA remains after matching mqpar.xml entries to pipeline-generated sample FASTAs"
    return canonical_entries.collect { resolveExistingFasta(it, mqpar_file, maxquant_txt) }
}

process DOWNLOAD_REFERENCES {
    tag 'GRCh38_Ensembl111'
    cpus 4; memory '16 GB'; time '12h'; disk '150 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    input:
    path genome_archive
    path gtf_archive
    path cdna_archive
    path proteome_archive
    path vep_archive
    path arriba_archive
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

    mkdir -p refs/vep_cache refs/arriba_unpack

    gzip -dc ${genome_archive} > refs/genome.fa
    gzip -dc ${gtf_archive} > refs/genes.gtf
    gzip -dc ${cdna_archive} > refs/cdna.fa
    gzip -dc ${proteome_archive} > refs/human_reviewed_isoforms.fasta

    tar -xzf ${vep_archive} -C refs/vep_cache
    tar -xzf ${arriba_archive} -C refs/arriba_unpack

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
    cpus 12; memory '32 GB'; time '36h'; disk '750 GB'
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
    cpus 1; memory '4 GB'; time '12h'; disk '1000 GB'
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
    cpus 4; memory '12 GB'; time '16h'; disk '500 GB'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc/fastqc_raw", mode:'copy'
    input: tuple val(meta), path(r1), path(r2)
    output:
    path "${meta.sample}.raw_fastqc", emit: qc
    script:
    """
    mkdir ${meta.sample}.raw_fastqc
    ln -s ${r1} ${meta.sample}.raw.R1.fastq.gz
    ln -s ${r2} ${meta.sample}.raw.R2.fastq.gz
    fastqc --threads ${task.cpus} --outdir ${meta.sample}.raw_fastqc ${meta.sample}.raw.R1.fastq.gz ${meta.sample}.raw.R2.fastq.gz
    """
}

process TRIM_GALORE {
    tag "${meta.sample}"
    cpus 8; memory '24 GB'; time '36h'; disk '750 GB'
    container 'quay.io/biocontainers/trim-galore:0.6.10--hdfd78af_0'
    publishDir "${params.outdir}/qc/trim_galore", mode:'copy', pattern:'*_trimming_report.txt'
    input: tuple val(meta), path(r1), path(r2)
    output:
    tuple val(meta), path("${meta.sample}_R1.trimmed.fastq.gz"), path("${meta.sample}_R2.trimmed.fastq.gz"), emit: reads
    path '*_trimming_report.txt', emit: reports
    script:
    def trimCores = Math.max(1, task.cpus.intdiv(2))
    """
    trim_galore --paired --quality 20 --length 36 --cores ${trimCores} --gzip --basename ${meta.sample} ${r1} ${r2}
    mv ${meta.sample}_val_1.fq.gz ${meta.sample}_R1.trimmed.fastq.gz
    mv ${meta.sample}_val_2.fq.gz ${meta.sample}_R2.trimmed.fastq.gz
    """
}

process FASTQC_TRIMMED {
    tag "${meta.sample}:trimmed"
    cpus 4; memory '12 GB'; time '16h'; disk '500 GB'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/qc/fastqc_trimmed", mode:'copy'
    input: tuple val(meta), path(r1), path(r2)
    output:
    path "${meta.sample}.trimmed_fastqc", emit: qc
    script:
    """
    mkdir ${meta.sample}.trimmed_fastqc
    ln -s ${r1} ${meta.sample}.trimmed.R1.fastq.gz
    ln -s ${r2} ${meta.sample}.trimmed.R2.fastq.gz
    fastqc --threads ${task.cpus} --outdir ${meta.sample}.trimmed_fastqc ${meta.sample}.trimmed.R1.fastq.gz ${meta.sample}.trimmed.R2.fastq.gz
    """
}

process STAR_INDEX {
    tag 'GRCh38_Ensembl111'
    cpus 16; memory '64 GB'; time '16h'; disk '400 GB'
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
    cpus 20; memory '128 GB'; time '48h'; disk '1000 GB'
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
    cpus 12; memory '48 GB'; time '36h'; disk '1000 GB'
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
    cpus 4; memory '8 GB'; time '8h'; disk '50 GB'
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
    tag 'GRCh38'; cpus 4; memory '16 GB'; time '4h'; disk '30 GB'
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
    cpus 2; memory '48 GB'; time '8h'; disk '250 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: tuple val(meta), path(bam), path(bai)
    output:
    tuple val(meta), path("${meta.sample}.markdup.bam"), path("${meta.sample}.markdup.bam.bai"), emit: bam
    path "${meta.sample}.metrics.txt", emit: metrics
    script:
    def javaHeapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.80) as int)
    def javaGcThreads = Math.max(1, Math.min(task.cpus as int, 8))
    """
    set -euo pipefail
    mkdir -p gatk_tmp
    trap 'rm -rf gatk_tmp' EXIT
    gatk --java-options "-Xms1g -Xmx${javaHeapGb}g -XX:ParallelGCThreads=${javaGcThreads} -Djava.io.tmpdir=\${PWD}/gatk_tmp" \
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
    tag "${meta.sample}"
    cpus 2; memory '24 GB'; time '16h'; disk '400 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: tuple val(meta), path(bam), path(bai); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.split.bam"), path("${meta.sample}.split.bam.bai")
    script:
    def javaHeapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.80) as int)
    """
    set -euo pipefail
    mkdir -p gatk_tmp
    trap 'rm -rf gatk_tmp' EXIT
    gatk --java-options "-Xms1g -Xmx${javaHeapGb}g -Djava.io.tmpdir=\${PWD}/gatk_tmp" SplitNCigarReads \
        -R ${genome} -I ${bam} -O ${meta.sample}.split.bam \
        --create-output-bam-index true
    if [[ -s ${meta.sample}.split.bai && ! -e ${meta.sample}.split.bam.bai ]]; then
        mv ${meta.sample}.split.bai ${meta.sample}.split.bam.bai
    fi
    test -s ${meta.sample}.split.bam.bai
    """
}

process PREPARE_HAPLOTYPE_INTERVALS {
    tag "GRCh38:${params.haplotype_scatter_count}_shards"
    cpus 2; memory '8 GB'; time '4h'; disk '50 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: tuple path(genome), path(fai), path(dict)
    output: path 'hc_intervals/*.interval_list', emit: intervals
    script:
    def javaHeapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.80) as int)
    """
    set -euo pipefail
    mkdir hc_intervals
    gatk --java-options "-Xms1g -Xmx${javaHeapGb}g" SplitIntervals \
        -R ${genome} \
        -O hc_intervals \
        --scatter-count ${params.haplotype_scatter_count} \
        --subdivision-mode INTERVAL_SUBDIVISION
    test "\$(find hc_intervals -name '*.interval_list' -type f | wc -l)" -eq ${params.haplotype_scatter_count}
    """
}

process HAPLOTYPE_CALLER {
    tag "${meta.sample}:${interval.baseName}"
    cpus 8; memory '20 GB'; time '24h'; disk '150 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input:
    tuple val(meta), path(bam), path(bai), path(interval)
    tuple path(genome), path(fai), path(dict)
    output:
    tuple val(meta), val(interval.baseName), path("${meta.sample}.${interval.baseName}.g.vcf.gz"), path("${meta.sample}.${interval.baseName}.g.vcf.gz.tbi")
    script:
    def javaHeapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.80) as int)
    def softClipArg = strictBooleanParam(params.hc_dont_use_soft_clipped_bases, '--hc_dont_use_soft_clipped_bases') ? '--dont-use-soft-clipped-bases true' : '--dont-use-soft-clipped-bases false'
    def pcrIndelArg = params.hc_pcr_indel_model ? "--pcr-indel-model ${params.hc_pcr_indel_model}" : ''
    """
    set -euo pipefail
    mkdir -p gatk_tmp
    trap 'rm -rf gatk_tmp' EXIT
    gatk --java-options "-Xms1g -Xmx${javaHeapGb}g -Djava.io.tmpdir=\${PWD}/gatk_tmp" HaplotypeCaller \
        -R ${genome} \
        -I ${bam} \
        -L ${interval} \
        -O ${meta.sample}.${interval.baseName}.g.vcf.gz \
        -ERC GVCF \
        ${softClipArg} \
        ${pcrIndelArg} \
        --standard-min-confidence-threshold-for-calling ${params.hc_calling_confidence} \
        --native-pair-hmm-threads ${task.cpus}
    """
}

process VALIDATE_HAPLOTYPE_SHARDS {
    tag "${meta.sample}:validate_haplotype_shards"
    cpus 1; memory '2 GB'; time '2h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/qc/haplotype_shards", mode:'copy'
    input:
    tuple val(meta), path(gvcfs), path(tbis)
    path validator
    output:
    tuple val(meta), path(gvcfs), path(tbis), path("${meta.sample}.haplotype_shards.tsv"), path("${meta.sample}.haplotype_shards.summary.txt")
    script:
    """
    python ${validator} --sample ${meta.sample} --expected ${params.haplotype_scatter_count} --gvcf ${gvcfs} --tbi ${tbis} --output-prefix ${meta.sample}.haplotype_shards
    """
}

process GATHER_HAPLOTYPE_GVCF {
    tag "${meta.sample}:gather_gvcf"
    cpus 2; memory '8 GB'; time '8h'; disk '200 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/gvcf", mode:'copy'
    input:
    tuple val(meta), path(gvcfs), path(tbis), path(shard_table), path(shard_summary)
    output:
    tuple val(meta), path("${meta.sample}.g.vcf.gz"), path("${meta.sample}.g.vcf.gz.tbi")
    script:
    def ordered = gvcfs.sort { it.name }
    def inputs = ordered.collect { "-I ${it}" }.join(' ')
    def javaHeapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.80) as int)
    """
    set -euo pipefail
    gatk --java-options "-Xms1g -Xmx${javaHeapGb}g" GatherVcfs \
        ${inputs} \
        -O ${meta.sample}.g.vcf.gz \
        --CREATE_INDEX false
    gatk --java-options "-Xms1g -Xmx${javaHeapGb}g" IndexFeatureFile \
        -I ${meta.sample}.g.vcf.gz
    test -s ${meta.sample}.g.vcf.gz
    test -s ${meta.sample}.g.vcf.gz.tbi
    """
}

process GENOTYPE_VARIANTS {
    tag "${meta.sample}:genotype"
    cpus 2; memory '8 GB'; time '8h'; disk '80 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_raw", mode:'copy'
    input: tuple val(meta), path(gvcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.raw.vcf.gz"), path("${meta.sample}.raw.vcf.gz.tbi")
    script:
    def javaHeapGb = Math.max(1, Math.floor(task.memory.toGiga() * 0.80) as int)
    """
    set -euo pipefail
    gatk --java-options "-Xms1g -Xmx${javaHeapGb}g" GenotypeGVCFs \
        -R ${genome} -V ${gvcf} -O ${meta.sample}.raw.vcf.gz
    test -s ${meta.sample}.raw.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.raw.vcf.gz
    test -s ${meta.sample}.raw.vcf.gz
    test -s ${meta.sample}.raw.vcf.gz.tbi
    """
}

process NORMALIZE_VARIANTS {
    tag "${meta.sample}:normalize"
    cpus 2; memory '8 GB'; time '8h'; disk '80 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_normalized", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.normalized.vcf.gz"), path("${meta.sample}.normalized.vcf.gz.tbi")
    script:
    """
    set -euo pipefail
    gatk LeftAlignAndTrimVariants -R ${genome} -V ${vcf} \
        -O ${meta.sample}.normalized.vcf.gz --split-multi-allelics true
    test -s ${meta.sample}.normalized.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.normalized.vcf.gz
    """
}

process SELECT_SNPS {
    tag "${meta.sample}:snps"
    cpus 1; memory '4 GB'; time '4h'; disk '40 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_snp/raw", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.snp.raw.vcf.gz"), path("${meta.sample}.snp.raw.vcf.gz.tbi")
    script:
    """
    set -euo pipefail
    gatk SelectVariants -R ${genome} -V ${vcf} --select-type-to-include SNP -O ${meta.sample}.snp.raw.vcf.gz
    test -s ${meta.sample}.snp.raw.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.snp.raw.vcf.gz
    """
}

process SELECT_INDELS {
    tag "${meta.sample}:indels"
    cpus 1; memory '4 GB'; time '4h'; disk '40 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_indel/raw", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.indel.raw.vcf.gz"), path("${meta.sample}.indel.raw.vcf.gz.tbi")
    script:
    """
    set -euo pipefail
    gatk SelectVariants -R ${genome} -V ${vcf} --select-type-to-include INDEL -O ${meta.sample}.indel.raw.vcf.gz
    test -s ${meta.sample}.indel.raw.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.indel.raw.vcf.gz
    """
}

process FILTER_SNPS {
    tag "${meta.sample}:filter_snps"
    cpus 1; memory '4 GB'; time '4h'; disk '40 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_snp/filtered", mode:'copy', pattern:'*.snp.filtered.vcf.gz*'
    publishDir "${params.outdir}/vcf_snp/pass", mode:'copy', pattern:'*.snp.pass.vcf.gz*'
    input: tuple val(meta), path(vcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.snp.filtered.vcf.gz"), path("${meta.sample}.snp.filtered.vcf.gz.tbi"), path("${meta.sample}.snp.pass.vcf.gz"), path("${meta.sample}.snp.pass.vcf.gz.tbi")
    script:
    """
    set -euo pipefail
    gatk VariantFiltration -R ${genome} -V ${vcf} \
        --filter-name SNP_QD2 --filter-expression 'QD < ${params.snp_filter_qd}' \
        --filter-name SNP_FS60 --filter-expression 'FS > ${params.snp_filter_fs}' \
        --filter-name SNP_SOR3 --filter-expression 'SOR > ${params.snp_filter_sor}' \
        --filter-name SNP_MQ40 --filter-expression 'MQ < ${params.snp_filter_mq}' \
        --filter-name SNP_MQRankSum-12.5 --filter-expression 'MQRankSum < ${params.snp_filter_mq_rank_sum}' \
        --filter-name SNP_ReadPosRankSum-8 --filter-expression 'ReadPosRankSum < ${params.snp_filter_read_pos_rank_sum}' \
        -O ${meta.sample}.snp.filtered.vcf.gz
    test -s ${meta.sample}.snp.filtered.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.snp.filtered.vcf.gz
    gatk SelectVariants -R ${genome} -V ${meta.sample}.snp.filtered.vcf.gz --exclude-filtered -O ${meta.sample}.snp.pass.vcf.gz
    test -s ${meta.sample}.snp.pass.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.snp.pass.vcf.gz
    """
}

process FILTER_INDELS {
    tag "${meta.sample}:filter_indels"
    cpus 1; memory '4 GB'; time '4h'; disk '40 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_indel/filtered", mode:'copy', pattern:'*.indel.filtered.vcf.gz*'
    publishDir "${params.outdir}/vcf_indel/pass", mode:'copy', pattern:'*.indel.pass.vcf.gz*'
    input: tuple val(meta), path(vcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.indel.filtered.vcf.gz"), path("${meta.sample}.indel.filtered.vcf.gz.tbi"), path("${meta.sample}.indel.pass.vcf.gz"), path("${meta.sample}.indel.pass.vcf.gz.tbi")
    script:
    """
    set -euo pipefail
    gatk VariantFiltration -R ${genome} -V ${vcf} \
        --filter-name INDEL_QD2 --filter-expression 'QD < ${params.indel_filter_qd}' \
        --filter-name INDEL_FS200 --filter-expression 'FS > ${params.indel_filter_fs}' \
        --filter-name INDEL_SOR10 --filter-expression 'SOR > ${params.indel_filter_sor}' \
        --filter-name INDEL_ReadPosRankSum-20 --filter-expression 'ReadPosRankSum < ${params.indel_filter_read_pos_rank_sum}' \
        -O ${meta.sample}.indel.filtered.vcf.gz
    test -s ${meta.sample}.indel.filtered.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.indel.filtered.vcf.gz
    gatk SelectVariants -R ${genome} -V ${meta.sample}.indel.filtered.vcf.gz --exclude-filtered -O ${meta.sample}.indel.pass.vcf.gz
    test -s ${meta.sample}.indel.pass.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.indel.pass.vcf.gz
    """
}

process MERGE_FILTERED_VARIANTS {
    tag "${meta.sample}:merge_types"
    cpus 1; memory '4 GB'; time '4h'; disk '80 GB'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_filtered", mode:'copy', pattern:'*.filtered.vcf.gz*'
    publishDir "${params.outdir}/vcf_pass", mode:'copy', pattern:'*.pass.vcf.gz*'
    input:
    tuple val(meta), path(snp_filtered), path(snp_filtered_tbi), path(snp_pass), path(snp_pass_tbi), path(indel_filtered), path(indel_filtered_tbi), path(indel_pass), path(indel_pass_tbi)
    output:
    tuple val(meta), path("${meta.sample}.filtered.vcf.gz"), path("${meta.sample}.filtered.vcf.gz.tbi"), emit: filtered
    tuple val(meta), path("${meta.sample}.pass.vcf.gz"), path("${meta.sample}.pass.vcf.gz.tbi"), emit: pass
    script:
    """
    set -euo pipefail
    gatk MergeVcfs -I ${snp_filtered} -I ${indel_filtered} -O ${meta.sample}.filtered.vcf.gz
    test -s ${meta.sample}.filtered.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.filtered.vcf.gz
    gatk MergeVcfs -I ${snp_pass} -I ${indel_pass} -O ${meta.sample}.pass.vcf.gz
    test -s ${meta.sample}.pass.vcf.gz.tbi || gatk IndexFeatureFile -I ${meta.sample}.pass.vcf.gz
    """
}

process VARIANT_STAGE_QC {
    tag "${meta.sample}:variant_stage_qc"
    cpus 1; memory '2 GB'; time '2h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/qc/variant_stages", mode:'copy'
    input:
    tuple val(meta), path(raw_vcf), path(raw_tbi), path(pass_vcf), path(pass_tbi), path(rna_vcf), path(rna_tbi)
    path genome
    path summarizer
    output:
    path "${meta.sample}.variant_stages.tsv", emit: table
    path "${meta.sample}.variant_stages.report.md", emit: report
    path "${meta.sample}.variant_stages.provenance.tsv", emit: provenance
    script:
    """
    python ${summarizer} --sample ${meta.sample} --raw ${raw_vcf} --pass-vcf ${pass_vcf} --rna ${rna_vcf} --genome ${genome} --calling-confidence ${params.hc_calling_confidence} --soft-clipped-setting ${booleanText(params.hc_dont_use_soft_clipped_bases, '--hc_dont_use_soft_clipped_bases')} --pcr-indel-model ${params.hc_pcr_indel_model} --output-prefix ${meta.sample}.variant_stages
    """
}

process COMPARE_EXTERNAL_VCF {
    tag "${meta.sample}:${meta.comparison_stage}:external_vcf_comparison"
    cpus 1; memory '4 GB'; time '4h'; disk '30 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/comparison/external_vcf", mode:'copy'
    input:
    tuple val(meta), path(pgtk_vcf), path(pgtk_tbi), path(external_vcf)
    path script_file
    output:
    path "${meta.sample}.${meta.comparison_stage}.external_comparison.summary.tsv", emit: summary
    path "${meta.sample}.${meta.comparison_stage}.external_comparison.shared.tsv", emit: shared
    path "${meta.sample}.${meta.comparison_stage}.external_comparison.pgtk_only.tsv", emit: pgtk_only
    path "${meta.sample}.${meta.comparison_stage}.external_comparison.external_only.tsv", emit: external_only
    path "${meta.sample}.${meta.comparison_stage}.external_comparison.report.md", emit: report
    script:
    """
    python ${script_file} --sample ${meta.sample} --stage ${meta.comparison_stage} --pgtk ${pgtk_vcf} --external ${external_vcf} --output-prefix ${meta.sample}.${meta.comparison_stage}.external_comparison
    """
}

process BCFTOOLS_STATS {
    tag "${meta.sample}"
    cpus 2; memory '4 GB'; time '8h'; disk '50 GB'
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
    tag "${meta.sample}"; cpus 12; memory '8 GB'; time '24h'; disk '100 GB'
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

process VALIDATE_RNA_VARIANTS {
    tag "${meta.sample}"
    cpus 1; memory '4 GB'; time '16h'; disk '100 GB'
    container "${params.pysam_image}"
    publishDir "${params.outdir}/rna_validation/variants", mode:'copy'
    input:
    tuple val(meta), path(vcf), path(tbi)
    path genome
    path validator
    output:
    tuple val(meta), path("${meta.sample}.rna.validated.vcf.gz"), path("${meta.sample}.rna.validated.vcf.gz.tbi"), emit: validated
    tuple val(meta), path("${meta.sample}.rna.audit.tsv"), emit: audit
    tuple val(meta), path("${meta.sample}.rna.rejected.tsv"), emit: rejected
    script:
    """
    set -euo pipefail
    python3 ${validator} variant \\
        --input ${vcf} --genome ${genome} --sample ${meta.sample} \\
        --min-depth ${params.rna_variant_min_depth} \\
        --min-alt-reads ${params.rna_variant_min_alt_reads} \\
        --min-alt-fraction ${params.rna_variant_min_alt_fraction} \\
        --output-prefix ${meta.sample}.rna
    """
}

process VALIDATE_RNA_FUSIONS {
    tag "${meta.sample}"
    cpus 1; memory '4 GB'; time '8h'; disk '50 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/rna_validation/fusions", mode:'copy'
    input:
    tuple val(meta), path(fusions)
    path validator
    output:
    tuple val(meta), path("${meta.sample}.fusion.validated.tsv"), emit: validated
    tuple val(meta), path("${meta.sample}.fusion.audit.tsv"), emit: audit
    tuple val(meta), path("${meta.sample}.fusion.rejected.tsv"), emit: rejected
    script:
    """
    set -euo pipefail
    python ${validator} fusion \\
        --input ${fusions} --sample ${meta.sample} \\
        --min-split-reads ${params.rna_fusion_min_split_reads} \\
        --min-total-support ${params.rna_fusion_min_total_support} \\
        --output-prefix ${meta.sample}.fusion
    """
}

process VALIDATE_RNA_SPLICE_TRANSCRIPTS {
    tag "${meta.sample}"
    cpus 2; memory '8 GB'; time '24h'; disk '150 GB'
    container "${params.pysam_image}"
    publishDir "${params.outdir}/rna_validation/splicing", mode:'copy'
    input:
    tuple val(meta), path(novel_gtf), path(bam), path(bai)
    path validator
    output:
    tuple val(meta), path("${meta.sample}.splice.validated.gtf"), emit: validated
    tuple val(meta), path("${meta.sample}.splice.audit.tsv"), emit: audit
    tuple val(meta), path("${meta.sample}.splice.rejected.tsv"), emit: rejected
    script:
    """
    set -euo pipefail
    python3 ${validator} splice \\
        --input ${novel_gtf} --bam ${bam} --sample ${meta.sample} \\
        --min-junction-reads ${params.splice_min_junction_reads} \\
        --output-prefix ${meta.sample}.splice
    """
}

process VALIDATE_VARIANT_CODONS {
    tag "${meta.sample}:genome_read_codon_validation"
    cpus 2; memory '8 GB'; time '8h'; disk '100 GB'
    container "${params.pysam_image}"
    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai)
    path genome
    path validation_script
    output:
    tuple val(meta), path("${meta.sample}.variant_codon_validation.all.tsv"), path("${meta.sample}.variant_codon_validation.validated.tsv"), path("${meta.sample}.variant_codon_validation.partial.tsv"), path("${meta.sample}.variant_codon_validation.failed.tsv"), path("${meta.sample}.variant_codon_validation.category_summary.tsv"), path("${meta.sample}.variant_codon_validation.summary.txt"), path("${meta.sample}.variant_codon_validation.report.md")
    script:
    """
    set -euo pipefail
    python3 ${validation_script} \
        --vcf ${vcf} \
        --bam ${meta.sample}=${bam} \
        --genome ${genome} \
        --threads ${task.cpus} \
        --min-base-quality 20 \
        --min-mapping-quality 20 \
        --min-alt-reads ${params.rna_variant_min_alt_reads} \
        --min-alt-fraction ${params.rna_variant_min_alt_fraction} \
        --output-prefix ${meta.sample}.variant_codon_validation
    """
}

process VALIDATE_VARIANT_READ_PROVENANCE {
    tag "${meta.sample}:variant_supporting_read_provenance"
    cpus 2; memory '12 GB'; time '8h'; disk '100 GB'
    container "${params.pysam_image}"
    input:
    tuple val(meta), path(vcf), path(tbi), path(bam), path(bai)
    path samplesheet
    path validation_script
    output:
    tuple val(meta), path("${meta.sample}.variant_read_provenance.supporting_reads.tsv"), path("${meta.sample}.variant_read_provenance.summary.txt"), path("${meta.sample}.variant_read_provenance.report.md")
    script:
    """
    set -euo pipefail
    python3 ${validation_script} \
        --vcf ${vcf} \
        --bam ${meta.sample}=${bam} \
        --samples ${samplesheet} \
        --threads ${task.cpus} \
        --min-base-quality 20 \
        --min-mapping-quality 20 \
        --output-prefix ${meta.sample}.variant_read_provenance
    """
}

process MERGE_VARIANT_CODON_VALIDATION {
    tag 'merge_variant_codon_validation'
    cpus 1; memory '4 GB'; time '6h'; disk '100 GB'
    publishDir "${params.outdir}/rna_validation/variant_codons", mode:'copy'
    input:
    path inputs
    path merger
    output:
    path 'variant_codon_validation.all.tsv', emit: all
    path 'variant_codon_validation.validated.tsv', emit: validated
    path 'variant_codon_validation.partial.tsv', emit: partial
    path 'variant_codon_validation.failed.tsv', emit: failed
    path 'variant_codon_validation.category_summary.tsv', emit: category_summary
    path 'variant_codon_validation.summary.txt', emit: summary
    path 'variant_codon_validation.report.md', emit: report
    script:
    """
    python ${merger} --mode codon --inputs ${inputs} --output-prefix variant_codon_validation
    """
}

process MERGE_VARIANT_READ_PROVENANCE {
    tag 'merge_variant_read_provenance'
    cpus 1; memory '4 GB'; time '6h'; disk '200 GB'
    publishDir "${params.outdir}/rna_validation/variant_read_provenance", mode:'copy'
    input:
    path inputs
    path merger
    output:
    path 'variant_read_provenance.supporting_reads.tsv', emit: reads
    path 'variant_read_provenance.summary.txt', emit: summary
    path 'variant_read_provenance.report.md', emit: report
    script:
    """
    python ${merger} --mode provenance --inputs ${inputs} --output-prefix variant_read_provenance
    """
}

process ANALYZE_CODON_MISMATCHES {
    tag 'codon_translation_mismatch_investigation'
    cpus 1; memory '2 GB'; time '4h'; disk '50 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/rna_validation/variant_codons", mode:'copy'
    input:
    path codon_table
    path analysis_script
    output:
    path 'codon_mismatch_analysis.detailed.tsv', emit: detailed
    path 'codon_mismatch_analysis.summary.tsv', emit: summary
    path 'codon_mismatch_analysis.manual_review.tsv', emit: manual_review
    path 'codon_mismatch_analysis.report.md', emit: report
    script:
    """
    set -euo pipefail
    python ${analysis_script} \
        --input ${codon_table} \
        --output-prefix codon_mismatch_analysis
    """
}

process PYPGATK_FASTA {
    tag "${meta.sample}"; cpus 1; memory '8 GB'; time '8h'; disk '150 GB'
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
    tag "${meta.sample}"; cpus 8; memory '32 GB'; time '12h'; disk '60 GB'
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
    cpus 2; memory '8 GB'; time '12h'; disk '100 GB'
    container "${params.container_cache}/pvactools-7.1.1.img"
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
    cpus 8; memory '32 GB'; time '36h'; disk '400 GB'
    container "${params.container_cache}/stringtie-3.0.3.img"
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
    cpus 2; memory '8 GB'; time '16h'; disk '150 GB'
    container "${params.container_cache}/gffcompare-0.12.10.img"
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
        for file in ./*; do test -f "\$file" || continue; bytes=\$(wc -c < "\$file"); printf '%s %s bytes\n' "\${file#./}" "\$bytes"; done | sort >&2
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
    cpus 8; memory '32 GB'; time '48h'; disk '300 GB'
    container "${params.container_cache}/transdecoder-6.0.0.img"
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
    cpus 1; memory '4 GB'; time '12h'; disk '100 GB'
    publishDir "${params.outdir}/combined_fasta", mode:'copy'
    input:
    tuple val(meta), path(variant_fasta), path(fusion_fasta), path(splice_fasta)
    output: tuple val(meta), path("${meta.sample}.exploratory_proteogenomics.fasta")
    script:
    """
    set -euo pipefail
    cat ${variant_fasta} ${fusion_fasta} ${splice_fasta} > ${meta.sample}.custom_events.raw.fasta
    awk '
        /^>/ {
            if (header != "" && sequence != "" && !seen[sequence]++) { print header; print sequence }
            header = \$0
            sequence = ""
            next
        }
        { sequence = sequence \$0 }
        END {
            if (header != "" && sequence != "" && !seen[sequence]++) { print header; print sequence }
        }
    ' ${meta.sample}.custom_events.raw.fasta > ${meta.sample}.exploratory_proteogenomics.fasta
    test -s ${meta.sample}.exploratory_proteogenomics.fasta
    if grep -q '^>CANONICAL|' ${meta.sample}.exploratory_proteogenomics.fasta; then
        echo 'ERROR: canonical sequence detected in custom event FASTA' >&2
        exit 1
    fi
    """
}

process PROGRESSION_SUBTRACT {
    tag "${meta.sample}"; cpus 2; memory '8 GB'; time '4h'; disk '50 GB'
    container 'quay.io/biocontainers/bcftools:1.21--h8b25389_0'
    publishDir "${params.outdir}/progression_vcf", mode:'copy'
    input: tuple val(meta), path(pv), path(pt), path(bv), path(bt)
    output:
    tuple val(meta),
        path("${meta.sample}.nonbaseline_only.vep.vcf.gz"), path("${meta.sample}.nonbaseline_only.vep.vcf.gz.tbi"),
        path("${meta.sample}.baseline_only.vep.vcf.gz"), path("${meta.sample}.baseline_only.vep.vcf.gz.tbi"),
        path("${meta.sample}.shared_with_baseline.vep.vcf.gz"), path("${meta.sample}.shared_with_baseline.vep.vcf.gz.tbi"),
        path("${meta.sample}.subtraction.summary.tsv")
    script:
    """
    set -euo pipefail
    mkdir isec
    bcftools isec -p isec -O z ${pv} ${bv}
    cp isec/0000.vcf.gz ${meta.sample}.nonbaseline_only.vep.vcf.gz
    cp isec/0001.vcf.gz ${meta.sample}.baseline_only.vep.vcf.gz
    cp isec/0002.vcf.gz ${meta.sample}.shared_with_baseline.vep.vcf.gz
    for f in ${meta.sample}.nonbaseline_only.vep.vcf.gz ${meta.sample}.baseline_only.vep.vcf.gz ${meta.sample}.shared_with_baseline.vep.vcf.gz; do bcftools index --tbi -f "\$f"; done
    printf 'Sample\tCategory\tRecords\n' > ${meta.sample}.subtraction.summary.tsv
    for category in nonbaseline_only baseline_only shared_with_baseline; do
        count=\$(bcftools view -H ${meta.sample}.\${category}.vep.vcf.gz | wc -l)
        printf '${meta.sample}\t%s\t%s\n' "\$category" "\$count" >> ${meta.sample}.subtraction.summary.tsv
    done
    """
}

process COUNT_GENES_PER_SAMPLE {
    tag "${meta.sample}:gene_counts"
    cpus 4; memory '4 GB'; time '4h'; disk '100 GB'
    container "${params.container_cache}/subread-2.0.8.img"
    publishDir "${params.outdir}/expression/per_sample", mode:'copy'
    input:
    tuple val(meta), path(bam), path(bai)
    path gtf
    output:
    tuple val(meta), path("${meta.sample}.gene_counts.tsv"), path("${meta.sample}.gene_counts.tsv.summary")
    script:
    def paired = strictBooleanParam(params.gene_count_count_read_pairs, '--gene_count_count_read_pairs') ? '-p --countReadPairs' : ''
    def bothEnds = strictBooleanParam(params.gene_count_require_both_ends, '--gene_count_require_both_ends') ? '-B' : ''
    def chimeric = strictBooleanParam(params.gene_count_exclude_chimeric, '--gene_count_exclude_chimeric') ? '-C' : ''
    def primary = strictBooleanParam(params.gene_count_primary_only, '--gene_count_primary_only') ? '--primary' : ''
    def overlap = strictBooleanParam(params.gene_count_allow_multi_overlap, '--gene_count_allow_multi_overlap') ? '-O' : ''
    def multi = strictBooleanParam(params.gene_count_count_multimapping, '--gene_count_count_multimapping') ? '-M' : ''
    """
    set -euo pipefail
    featureCounts \
        -T ${task.cpus} \
        -a ${gtf} \
        -o ${meta.sample}.gene_counts.tsv \
        -t ${params.gene_count_feature_type} \
        -g ${params.gene_count_id_attribute} \
        -s ${params.gene_count_strandedness} \
        -Q ${params.gene_count_min_mapq} \
        --minOverlap ${params.gene_count_min_overlap} \
        ${paired} ${bothEnds} ${chimeric} ${primary} ${overlap} ${multi} \
        ${bam}
    test -s ${meta.sample}.gene_counts.tsv
    test -s ${meta.sample}.gene_counts.tsv.summary
    """
}

process MERGE_GENE_EXPRESSION {
    tag 'merge_gene_expression'
    cpus 1; memory '8 GB'; time '4h'; disk '30 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/expression", mode:'copy'
    input:
    path count_tables
    path gtf
    path analysis_script
    output:
    path 'gene_expression.gene_expression.tsv', emit: matrix
    path 'gene_expression.summary.tsv', emit: summary
    script:
    """
    set -euo pipefail
    python3 ${analysis_script} merge-counts \
        --counts ${count_tables} \
        --gtf ${gtf} \
        --feature-type ${params.gene_count_feature_type} \
        --id-attribute ${params.gene_count_id_attribute} \
        --symbol-attribute ${params.gene_count_symbol_attribute} \
        --biotypes ${params.gene_count_biotypes} \
        --output-prefix gene_expression
    """
}

process ANALYZE_EXPRESSION_SAMPLE_GO {
    tag "${meta.sample}:expressed_GO"
    cpus 1; memory '2 GB'; time '1h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/expression/go/per_sample", mode:'copy'
    input:
    val(meta)
    path expression_matrix
    val all_samples
    path go_mapping
    path analysis_script
    output:
    tuple val(meta), path("${meta.sample}.expression_go.expression_ora.tsv"), path("${meta.sample}.expression_go.summary.tsv")
    script:
    """
    set -euo pipefail
    python3 ${analysis_script} sample-ora \
        --matrix ${expression_matrix} \
        --sample ${meta.sample} \
        --subject ${meta.tk} \
        --group ${meta.group} \
        --all-samples '${all_samples}' \
        --go-mapping ${go_mapping} \
        --background ${params.go_expression_background} \
        --cpm-threshold ${params.expression_cpm_threshold} \
        --tpm-threshold ${params.expression_tpm_threshold} \
        --go-min-size ${params.go_min_size} \
        --go-max-size ${params.go_max_size} \
        --fdr-threshold ${params.go_fdr_threshold} \
        --namespaces ${params.go_namespaces} \
        --output-prefix ${meta.sample}.expression_go
    """
}

process ANALYZE_EXPRESSION_RANKED_GO {
    tag "${meta.tk}:${meta.sample}_vs_${baseline_sample}:ranked_GO"
    cpus 1; memory '2 GB'; time '1h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/expression/go/ranked", mode:'copy'
    input:
    tuple val(meta), val(baseline_sample)
    path expression_matrix
    val all_samples
    path go_mapping
    path analysis_script
    output:
    tuple val(meta), val(baseline_sample), path("${meta.sample}_vs_${baseline_sample}.expression_go.ranked_go.tsv"), path("${meta.sample}_vs_${baseline_sample}.expression_go.summary.tsv")
    script:
    """
    set -euo pipefail
    python3 ${analysis_script} ranked-go \
        --matrix ${expression_matrix} \
        --sample ${meta.sample} \
        --baseline-sample ${baseline_sample} \
        --subject ${meta.tk} \
        --group ${meta.group} \
        --all-samples '${all_samples}' \
        --go-mapping ${go_mapping} \
        --background ${params.go_expression_background} \
        --pseudocount ${params.expression_pseudocount} \
        --rank-metric ${params.expression_rank_metric} \
        --min-nonzero-scores ${params.expression_rank_min_nonzero_scores} \
        --go-min-size ${params.go_min_size} \
        --go-max-size ${params.go_max_size} \
        --fdr-threshold ${params.go_fdr_threshold} \
        --namespaces ${params.go_namespaces} \
        --output-prefix ${meta.sample}_vs_${baseline_sample}.expression_go
    """
}

process MERGE_EXPRESSION_GO {
    tag 'merge_expression_GO'
    cpus 1; memory '2 GB'; time '1h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/expression/go", mode:'copy'
    input:
    path ora_tables
    path ranked_tables
    path summary_tables
    path samplesheet
    path analysis_script
    output:
    path 'expression_go.expression_ora.tsv', emit: ora
    path 'expression_go.ranked_go.tsv', emit: ranked
    path 'expression_go.summary.tsv', emit: summary
    script:
    """
    set -euo pipefail
    python3 ${analysis_script} merge-expression-go \
        --samples ${samplesheet} \
        --ora ${ora_tables} \
        --ranked ${ranked_tables} \
        --summary ${summary_tables} \
        --output-prefix expression_go
    """
}

process PREPARE_EXPRESSION_MULTIQC_CONTENT {
    tag 'expression_GO_multiqc_content'
    cpus 1; memory '2 GB'; time '2h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    input:
    path expression_ora
    path ranked_go
    path expression_summary
    path variant_set_go
    path variant_set_summary
    path content_builder
    output:
    path 'expression_multiqc_content'
    script:
    """
    set -euo pipefail
    python3 ${content_builder} --output-dir expression_multiqc_content --expression-ora ${expression_ora} --ranked-go ${ranked_go} --expression-summary ${expression_summary} --variant-set-go ${variant_set_go} --variant-set-summary ${variant_set_summary}
    """
}

process ANALYZE_PROGRESSION_VARIANT_SETS {
    tag 'progression_common_and_exclusive_GO'
    cpus 1; memory '8 GB'; time '8h'; disk '30 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/progression_biology/sets", mode:'copy'
    input:
    path sample_gene_tables
    path samplesheet
    path go_mapping
    path gtf
    path analysis_script
    output:
    path 'progression_variant_sets.variant_set_go.tsv', emit: enrichment
    path 'progression_variant_sets.summary.tsv', emit: summary
    script:
    """
    set -euo pipefail
    python3 ${analysis_script} variant-sets \
        --genes ${sample_gene_tables} \
        --background ${params.go_variant_background} \
        --samples ${samplesheet} \
        --go-mapping ${go_mapping} \
        --gtf ${gtf} \
        --id-attribute ${params.gene_count_id_attribute} \
        --symbol-attribute ${params.gene_count_symbol_attribute} \
        --biotypes ${params.go_variant_biotypes} \
        --go-min-size ${params.go_min_size} \
        --go-max-size ${params.go_max_size} \
        --fdr-threshold ${params.go_fdr_threshold} \
        --namespaces ${params.go_namespaces} \
        --output-prefix progression_variant_sets
    """
}

process PREPARE_GO_ANNOTATIONS {
    tag 'GO_annotations'
    cpus 1; memory '4 GB'; time '4h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/progression_biology/go_reference", mode:'copy'
    input: path go_obo; path go_gaf; path preparation_script
    output:
    path 'go_annotations.mapping.tsv', emit: mapping
    path 'go_annotations.metadata.tsv', emit: metadata
    script:
    """
    set -euo pipefail
    python3 ${preparation_script} --obo ${go_obo} --gaf ${go_gaf} --output-prefix go_annotations
    """
}
process ANALYZE_PROGRESSION_SAMPLE {
    tag "${meta.sample}:GO_progression"
    cpus 1; memory '4 GB'; time '4h'; disk '30 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/progression_biology/samples", mode:'copy'
    input:
    tuple val(meta), path(progression_vcf), path(progression_tbi), path(background_vcf), path(background_tbi)
    path go_mapping
    path analysis_script
    output:
    tuple val(meta), path("${meta.sample}.progression_biology.alleles.tsv"), path("${meta.sample}.progression_biology.genes.tsv"), path("${meta.sample}.progression_biology.go_enrichment.tsv"), path("${meta.sample}.progression_biology.candidates.tsv"), path("${meta.sample}.progression_biology.summary.tsv")
    script:
    """
    set -euo pipefail
    python3 ${analysis_script} --sample ${meta.sample} --subject ${meta.tk} --group ${meta.group} --progression-vcf ${progression_vcf} --background-vcf ${background_vcf} --go-mapping ${go_mapping} --go-min-size ${params.go_min_size} --go-max-size ${params.go_max_size} --fdr-threshold ${params.go_fdr_threshold} --output-prefix ${meta.sample}.progression_biology
    """
}
process COMPARE_PROGRESSION_PAIR {
    tag "${pair_meta.subject}:${pair_meta.sample_a}_vs_${pair_meta.sample_b}:GO"
    cpus 1; memory '2 GB'; time '2h'; disk '20 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/progression_biology/pairs", mode:'copy'
    input:
    tuple val(pair_meta), path(alleles_a), path(genes_a), path(go_a), path(alleles_b), path(genes_b), path(go_b)
    path comparison_script
    output:
    tuple val(pair_meta), path("${pair_meta.pair_id}.progression_pair.alleles.tsv"), path("${pair_meta.pair_id}.progression_pair.genes.tsv"), path("${pair_meta.pair_id}.progression_pair.go_contrasts.tsv"), path("${pair_meta.pair_id}.progression_pair.summary.tsv")
    script:
    """
    set -euo pipefail
    python3 ${comparison_script} --subject ${pair_meta.subject} --sample-a ${pair_meta.sample_a} --sample-b ${pair_meta.sample_b} --alleles-a ${alleles_a} --alleles-b ${alleles_b} --genes-a ${genes_a} --genes-b ${genes_b} --go-a ${go_a} --go-b ${go_b} --output-prefix ${pair_meta.pair_id}.progression_pair
    """
}
process MERGE_PROGRESSION_BIOLOGY {
    tag 'merge_GO_progression_biology'
    cpus 1; memory '4 GB'; time '4h'; disk '50 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/progression_biology", mode:'copy'
    input:
    path sample_alleles; path sample_genes; path sample_go; path sample_candidates; path sample_summaries
    path pair_alleles; path pair_genes; path pair_go; path pair_summaries
    path go_metadata; path merge_script
    output:
    path 'progression_biology.progression_alleles.tsv', emit: alleles
    path 'progression_biology.progression_genes.tsv', emit: genes
    path 'progression_biology.go_enrichment.tsv', emit: enrichment
    path 'progression_biology.pairwise_allele_contrasts.tsv', emit: pairwise_alleles
    path 'progression_biology.pairwise_gene_contrasts.tsv', emit: pairwise_genes
    path 'progression_biology.pairwise_go_contrasts.tsv', emit: pairwise_categories
    path 'progression_biology.pairwise_summary.tsv', emit: pairwise_summary
    path 'progression_biology.candidate_priority.tsv', emit: candidates
    path 'progression_biology.multiqc_summary.tsv', emit: multiqc_summary
    path 'progression_biology.go_metadata.tsv', emit: go_metadata_output
    path 'progression_biology.report.md', emit: report
    script:
    """
    set -euo pipefail
    python3 ${merge_script} --sample-alleles ${sample_alleles} --sample-genes ${sample_genes} --sample-go ${sample_go} --sample-candidates ${sample_candidates} --sample-summary ${sample_summaries} --pair-alleles ${pair_alleles} --pair-genes ${pair_genes} --pair-go ${pair_go} --pair-summary ${pair_summaries} --go-metadata ${go_metadata} --output-prefix progression_biology
    """
}
process BUILD_IGV_EVIDENCE_BUNDLE {
    tag 'all_rna_progression_igv'
    cpus 2; memory '16 GB'; time '24h'; disk '300 GB'
    container "${params.pysam_image}"
    publishDir "${params.outdir}/igv/all_evidence", mode:'copy'
    input:
    path rna_vcfs
    path progression_vcfs
    path fusion_tables
    path splice_tables
    path sorted_bams
    path genome
    path bundle_script
    output:
    path 'pgtk_igv.events.tsv', emit: events
    path 'pgtk_igv.events.bed', emit: bed
    path 'pgtk_igv.events.bedpe', emit: bedpe
    path 'pgtk_igv.sample_manifest.tsv', emit: manifest
    path 'pgtk_igv.igv.batch.txt', emit: batch
    path 'pgtk_igv.igv.session.xml', emit: session
    path 'pgtk_igv.summary.txt', emit: summary
    path 'pgtk_igv.*.events.bam*', emit: bams
    script:
    def bamArgs=sorted_bams.findAll { it.name.endsWith('.bam') }.collect { bam -> "--bam ${bam.baseName.tokenize('.')[0]}=${bam}" }.join(' ')
    """
    python3 ${bundle_script} --genome ${genome} --rna-vcf ${rna_vcfs} --progression-vcf ${progression_vcfs} --fusion-table ${fusion_tables} --splice-table ${splice_tables} ${bamArgs} --padding ${params.read_validation_padding} --output-prefix pgtk_igv
    """
}

process BUILD_FINDING_IGV_REVIEWS {
    tag 'consolidated_strict_finding_igv_bundle'
    cpus 2; memory '12 GB'; time '48h'; disk '500 GB'
    container "${params.pysam_image}"
    publishDir "${params.outdir}/igv/findings", mode:'copy'
    input:
    path events
    path event_bams
    path genome
    path review_script
    output:
    path 'finding_reviews', emit: reviews
    script:
    def bamFiles = event_bams.findAll { it.name ==~ /^pgtk_igv\.[^.]+\.events\.bam$/ }
    if (!bamFiles) error 'BUILD_FINDING_IGV_REVIEWS received no event BAM files'
    def bamArgs = bamFiles.collect { bam ->
        def matcher = bam.name =~ /^pgtk_igv\.([^.]+)\.events\.bam$/
        if (!matcher.matches()) error "Cannot derive sample identifier from event BAM: ${bam.name}"
        "--bam '${matcher.group(1)}=${bam}'"
    }.join(' ')
    """
    set -euo pipefail
    mkdir -p finding_reviews
    python3 ${review_script} \\
        --events ${events} \\
        ${bamArgs} \\
        --genome ${genome} \\
        --output-dir finding_reviews \\
        --padding ${params.read_validation_padding} \\
        --mapq ${params.finding_review_mapq} \\
        --baseq ${params.finding_review_baseq} \\
        --reference-display-reads ${params.finding_review_reference_reads} \\
        --alt-display-reads ${params.igv_report_max_reads} \\
        --finding-classes '${params.igv_report_classes}' \\
        --primary-class-order '${params.finding_primary_class_order}' \\
        --priority-mode '${params.finding_priority_mode}' \\
        --priority-genes '${params.finding_priority_genes}' \\
        --priority-impacts '${params.finding_priority_impacts}' \\
        --priority-consequences '${params.finding_priority_consequences}' \\
        --priority-limit 0 \\
        --gene-filter '${params.igv_report_gene_filter}' \\
        --sample-filter '${params.igv_report_sample_filter}'
    test -s finding_reviews/findings_manifest.tsv
    test -s finding_reviews/bam_manifest.tsv
    test -s finding_reviews/support_labels.bed
    test -s finding_reviews/priority_findings.bed
    test -s finding_reviews/review.igv.batch.txt
    test -s finding_reviews/igv.session.xml
    finding_count=\$(awk 'END { print NR - 1 }' finding_reviews/findings_manifest.tsv)
    test "\$finding_count" -gt 0
    bam_rows=\$(awk 'END { print NR - 1 }' finding_reviews/bam_manifest.tsv)
    sample_count=\$(( bam_rows / 4 ))
    entry_count=\$(find finding_reviews -mindepth 1 | wc -l)
    maximum_entries=\$(( 20 + 8 * sample_count ))
    if test "\$entry_count" -gt "\$maximum_entries"; then
        echo "ERROR: consolidated finding bundle has \$entry_count entries; expected at most \$maximum_entries" >&2
        exit 1
    fi
    printf 'Generated one consolidated strict IGV bundle for %s findings using %s filesystem entries\n' "\$finding_count" "\$entry_count"
    """
}

process ANALYZE_VARIANT_LANDSCAPE {
    tag 'variant_types_and_nonsynonymous_GO'
    cpus 1; memory '8 GB'; time '2h'; disk '30 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/variant_landscape", mode:'copy', pattern:'variant_landscape.*'
    input:
    path raw_vcfs
    path normalized_vcfs
    path filtered_vcfs
    path pass_vcfs
    path vep_vcfs
    path rna_vcfs
    path nonbaseline_vcfs
    path baseline_only_vcfs
    path shared_vcfs
    path go_mapping
    path analysis_script
    output:
    path 'variant_landscape.summary.tsv', emit: summary
    path 'variant_landscape.nonsynonymous_genes.tsv', emit: genes
    path 'variant_landscape.go_significant.tsv', emit: go_significant
    path 'variant_landscape.go_top.tsv', emit: go_top
    path 'variant_landscape.go_summary.tsv', emit: go_summary
    path 'variant_landscape.report.md', emit: report
    path 'variant_landscape.multiqc', emit: multiqc
    script:
    def args = []
    raw_vcfs.each { args << "--vcf raw_genotyped=${it}" }
    normalized_vcfs.each { args << "--vcf normalized=${it}" }
    filtered_vcfs.each { args << "--vcf hard_filtered_all=${it}" }
    pass_vcfs.each { args << "--vcf hard_filter_pass=${it}" }
    vep_vcfs.each { args << "--vcf vep_pass=${it}" }
    rna_vcfs.each { args << "--vcf rna_validated=${it}" }
    nonbaseline_vcfs.each { args << "--vcf progression_nonbaseline_only=${it}" }
    baseline_only_vcfs.each { args << "--vcf progression_baseline_only=${it}" }
    shared_vcfs.each { args << "--vcf progression_shared_with_baseline=${it}" }
    """
    set -euo pipefail
    python3 ${analysis_script} ${args.join(' ')} \
        --go-mapping ${go_mapping} \
        --go-min-size ${params.go_min_size} \
        --go-max-size ${params.go_max_size} \
        --fdr-threshold ${params.go_fdr_threshold} \
        --go-top 100 \
        --max-significant-rows 100000 \
        --max-output-mb 100 \
        --output-prefix variant_landscape
    """
}

process BUILD_FINDING_EXPLORER {
    tag 'complete_database_free_finding_explorer'
    cpus 1; memory '4 GB'; time '2h'; disk '20 GB'
    container "${params.container_cache}/quay.io-biocontainers-igv-reports-1.16.0--pyh7e72e81_0.img"
    publishDir "${params.outdir}/igv/findings", mode:'copy'
    input:
    path finding_reviews
    path genome
    path explorer_script
    path server_launcher
    output:
    path 'finding_explorer', emit: explorer
    script:
    """
    set -euo pipefail
    mkdir -p finding_explorer
    python3 ${explorer_script} --manifest ${finding_reviews}/findings_manifest.tsv --excluded-reads ${finding_reviews}/excluded_reads.tsv --bam-manifest ${finding_reviews}/bam_manifest.tsv --genome ${genome} --flanking ${params.read_validation_padding} --output-dir finding_explorer
    cp ${server_launcher} finding_explorer/serve_explorer.sh
    total=\$(awk -F': ' '\$1=="Findings" {print \$2}' finding_explorer/coverage_summary.txt)
    records=\$(awk -F': ' '\$1=="Embedded compact records" {print \$2}' finding_explorer/coverage_summary.txt)
    discarded=\$(awk -F': ' '\$1=="Findings discarded" {print \$2}' finding_explorer/coverage_summary.txt)
    databases=\$(awk -F': ' '\$1=="Database files" {print \$2}' finding_explorer/coverage_summary.txt)
    test "\$total" -gt 0
    test "\$total" -eq "\$records"
    test "\$discarded" -eq 0
    test "\$databases" -eq 0
    test -z "\$(find finding_explorer -type f -name '*.sqlite' -print -quit)"
    test -z "\$(find finding_explorer -type f -name '*.sqlite3' -print -quit)"
    test -z "\$(find finding_explorer -type f -name '*.db' -print -quit)"
    printf 'Generated database-free explorer for %s findings; discarded 0\n' "\$total"
    """
}

process BUILD_COMPARATIVE_ADVANTAGE_REPORT {
    tag 'comparative_biological_advantage'
    cpus 1; memory '4 GB'; time '8h'; disk '100 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/comparative_advantage", mode:'copy'
    input:
    path samplesheet
    path raw_vcfs
    path pass_vcfs
    path rna_vcfs
    path progression_vcfs
    path fusion_tables
    path splice_details
    path variant_fastas
    path fusion_fastas
    path splice_fastas
    path combined_fastas
    path stage_qc_tables
    path external_comparison_tables
    path codon_summary
    path provenance_summary
    path report_script
    output:
    path 'comparative_advantage.report.md', emit: report
    path 'comparative_advantage.variant_stage_inventory.tsv', emit: variant_inventory
    path 'comparative_advantage.fasta_inventory.tsv', emit: fasta_inventory
    path 'comparative_advantage.rna_event_inventory.tsv', emit: rna_inventory
    path 'comparative_advantage.external_caller_comparison.tsv', emit: external_comparison
    path 'comparative_advantage.multiqc_summary.tsv', emit: multiqc_summary
    script:
    """
    python ${report_script} --samples ${samplesheet} --raw-vcf ${raw_vcfs} --pass-vcf ${pass_vcfs} --rna-vcf ${rna_vcfs} --progression-vcf ${progression_vcfs} --fusion-table ${fusion_tables} --splice-detail ${splice_details} --variant-fasta ${variant_fastas} --fusion-fasta ${fusion_fastas} --splice-fasta ${splice_fastas} --combined-fasta ${combined_fastas} --stage-qc ${stage_qc_tables} --external-comparison ${external_comparison_tables} --codon-summary ${codon_summary} --provenance-summary ${provenance_summary} --output-prefix comparative_advantage
    """
}

process PREPARE_COMPARATIVE_MULTIQC_CONTENT {
    tag 'complete_core_multiqc_content'
    cpus 1; memory '4 GB'; time '4h'; disk '30 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    input:
    path samplesheet
    path report
    path variant_inventory
    path fasta_inventory
    path rna_inventory
    path external_comparison
    path summary
    path progression_report
    path progression_summary
    path progression_enrichment
    path progression_pairwise_categories
    path complete_report
    path rna_failure_report
    path rna_variant_explanations
    path report_builder
    val maxquant_enabled
    output: path 'comparative_multiqc_content'
    script:
    """
    set -euo pipefail
    python3 ${report_builder} --output-dir comparative_multiqc_content --samples ${samplesheet} --variant-inventory ${variant_inventory} --fasta-inventory ${fasta_inventory} --rna-inventory ${rna_inventory} --external-comparison ${external_comparison} --summary ${summary} --progression-summary ${progression_summary} --progression-enrichment ${progression_enrichment} --progression-pairwise ${progression_pairwise_categories} --complete-report ${complete_report} --rna-failure-report ${rna_failure_report} --rna-variant-explanations ${rna_variant_explanations} --comparative-report ${report} --progression-report ${progression_report} --maxquant-enabled ${maxquant_enabled}
    """
}

process MULTIQC_QC_DATA {
    tag 'qc_data_pass'
    cpus 4; memory '16 GB'; time '12h'; disk '150 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    input: path qc_files
    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data
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

process BUILD_COMPLETE_FINDINGS_REPORT {
    tag 'all_findings'
    cpus 1; memory '2 GB'; time '12h'; disk '100 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/reports", mode:'copy'
    input:
    path samplesheet
    path multiqc_data
    path variant_audits
    path fusion_audits
    path splice_audits
    path all_vep_vcfs
    path all_arriba
    path all_arriba_discarded
    path all_assembled
    path all_annotated_gtfs
    path all_novel_gtfs
    path variant_fastas
    path fusion_fastas
    path splice_fastas
    path progression_vcfs
    path reporter
    output:
    path 'complete_findings.report.md', emit: report
    path 'complete_findings.rna_event_audit.tsv', emit: audit
    path 'complete_findings.rna_validation_failures.tsv', emit: failures
    path 'complete_findings.rna_validation_failures.md', emit: failure_report
    path 'complete_findings.software_inventory.tsv', emit: inventory
    path 'complete_findings.multiqc_general_stats.tsv', emit: qc_stats
    path 'complete_findings.rna_variant_validation_explanations.md', emit: variant_explanations
    script:
    """
    set -euo pipefail
    python ${reporter} \\
      --samples ${samplesheet} \\
      --multiqc-data ${multiqc_data} \\
      --variant-audit ${variant_audits} \\
      --fusion-audit ${fusion_audits} \\
      --splice-audit ${splice_audits} \\
      --vep-vcf ${all_vep_vcfs} \\
      --arriba ${all_arriba} \\
      --arriba-discarded ${all_arriba_discarded} \\
      --assembled-gtf ${all_assembled} \\
      --annotated-gtf ${all_annotated_gtfs} \\
      --novel-gtf ${all_novel_gtfs} \\
      --variant-fasta ${variant_fastas} \\
      --fusion-fasta ${fusion_fastas} \\
      --splice-fasta ${splice_fastas} \\
      --progression-vcf ${progression_vcfs} \\
      --prefix complete_findings
    """
}

process VALIDATE_MAXQUANT_INPUTS {
    tag 'maxquant_inputs'
    cpus 1; memory '2 GB'; time '4h'; disk '30 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/proteogenomics_validation", mode:'copy'
    input:
    path peptides
    path evidence
    path msms
    path protein_groups
    path mqpar
    path canonical_fastas
    path contaminants_fasta
    output:
    path 'maxquant_inputs.validated.txt'
    script:
    """
    set -euo pipefail
    python - ${peptides} ${evidence} ${msms} ${protein_groups} ${mqpar} ${contaminants_fasta} ${canonical_fastas} <<'PY_MQ'
import csv
import pathlib
import sys
import xml.etree.ElementTree as ET

peptides, evidence, msms, protein_groups, mqpar, contaminants = map(pathlib.Path, sys.argv[1:7])
canonicals = [pathlib.Path(value) for value in sys.argv[7:]]
if not canonicals:
    raise SystemExit('no canonical FASTA files resolved')
for path in (peptides, evidence, msms, protein_groups, mqpar, contaminants, *canonicals):
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f'missing or empty input: {path}')
required = {
    peptides: {'Sequence', 'PEP', 'Score', 'MS/MS Count', 'Evidence IDs', 'MS/MS IDs', 'Protein group IDs'},
    evidence: {'id', 'Sequence', 'Raw file', 'Experiment', 'PEP', 'Score', 'MS/MS IDs'},
    msms: {'id', 'Sequence', 'Raw file', 'Scan number', 'Score', 'PEP'},
    protein_groups: {'id', 'Protein IDs'},
}
for table, expected in required.items():
    with table.open(encoding='utf-8', errors='replace', newline='') as handle:
        fields = set(csv.DictReader(handle, delimiter='\t').fieldnames or [])
    missing = expected - fields
    if missing:
        raise SystemExit(f'{table.name} missing columns: {sorted(missing)}')
root = ET.parse(mqpar).getroot()
version = (root.findtext('maxQuantVersion') or '').strip()
fastas = [(node.text or '').strip() for node in root.findall('./fastaFiles/FastaFileInfo/fastaFilePath')]
if not version or not fastas:
    raise SystemExit('mqpar.xml lacks MaxQuant version or searched FASTA paths')
if (root.findtext('includeContaminants') or '').strip().lower() != 'true':
    raise SystemExit('mqpar.xml does not report includeContaminants=True')
with open('maxquant_inputs.validated.txt', 'w', encoding='utf-8') as handle:
    print(f'MaxQuant version: {version}', file=handle)
    print(f'Searched FASTA files: {len(fastas)}', file=handle)
    print(f'Canonical FASTA files supplied: {len(canonicals)}', file=handle)
    for canonical in canonicals:
        print(f'Canonical FASTA: {canonical}', file=handle)
    print(f'Contaminant FASTA supplied: {contaminants}', file=handle)
    print(f'Match Between Runs: {(root.findtext("matchBetweenRuns") or "unknown").strip()}', file=handle)
    print(f'Match unidentified features: {(root.findtext("matchUnidentifiedFeatures") or "unknown").strip()}', file=handle)
    print(f'MBR FDR enabled: {(root.findtext("matchBetweenRunsFdr") or "unknown").strip()}', file=handle)
    print(f'Matching time window: {(root.findtext("matchingTimeWindow") or "unknown").strip()}', file=handle)
    print(f'Alignment time window: {(root.findtext("alignmentTimeWindow") or "unknown").strip()}', file=handle)
PY_MQ
    """
}

process MAP_MAXQUANT_PEPTIDES {
    tag 'peptide_fasta_mapping'
    cpus 2; memory '16 GB'; time '24h'; disk '200 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/proteogenomics_validation", mode:'copy'
    input:
    path validation_stamp
    path peptides
    path canonical_fastas
    path contaminants_fasta
    path combined_fastas
    path mapper_script
    output:
    path 'peptide_fasta_mapping.mapping.tsv', emit: mapping
    path 'peptide_fasta_mapping.candidates.tsv', emit: candidates
    path 'peptide_fasta_mapping.summary.txt', emit: summary
    script:
    """
    set -euo pipefail
    python ${mapper_script} \\
        --peptides ${peptides} \\
        --fasta ${canonical_fastas} ${contaminants_fasta} ${combined_fastas} \\
        --output-prefix peptide_fasta_mapping
    """
}

process ANNOTATE_MAXQUANT_VARIANTS {
    tag 'variant_peptide_annotation'
    cpus 2; memory '16 GB'; time '24h'; disk '200 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/proteogenomics_validation", mode:'copy'
    input:
    path validation_stamp
    path candidates
    path vep_vcfs
    path combined_fastas
    path ensembl_pep
    path annotation_script
    output:
    path 'variant_peptide_annotation.detailed.tsv', emit: detailed
    path 'variant_peptide_annotation.prioritized.tsv', emit: prioritized
    path 'variant_peptide_annotation.validated.tsv', emit: validated
    path 'variant_peptide_annotation.rejected.tsv', emit: rejected
    path 'variant_peptide_annotation.summary.txt', emit: summary
    path 'variant_peptide_annotation.unresolved.tsv', emit: unresolved
    script:
    """
    set -euo pipefail
    python ${annotation_script} \\
        --candidates ${candidates} \\
        --vep-vcf ${vep_vcfs} \\
        --variant-fasta ${combined_fastas} \\
        --ensembl-pep ${ensembl_pep} \\
        --il-equivalent \\
        --output-prefix variant_peptide_annotation
    """
}

process ANALYZE_MAXQUANT_JUNCTIONS {
    tag 'junction_peptide_analysis'
    cpus 2; memory '16 GB'; time '24h'; disk '200 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/proteogenomics_validation", mode:'copy'
    input:
    path validation_stamp
    path peptides
    path canonical_fastas
    path contaminants_fasta
    path fusion_fastas
    path splice_fastas
    path arriba_tables
    path junction_script
    output:
    path 'junction_peptide_analysis.all_mappings.tsv', emit: all_mappings
    path 'junction_peptide_analysis.fusion_candidates.tsv', emit: fusion_candidates
    path 'junction_peptide_analysis.splice_candidates.tsv', emit: splice_candidates
    path 'junction_peptide_analysis.inferred_junctions.tsv', emit: inferred
    path 'junction_peptide_analysis.summary.txt', emit: summary
    script:
    """
    set -euo pipefail
    python ${junction_script} \\
        --peptides ${peptides} \\
        --canonical-fasta ${canonical_fastas} ${contaminants_fasta} \\
        --fusion-fasta ${fusion_fastas} \\
        --splice-fasta ${splice_fastas} \\
        --arriba ${arriba_tables} \\
        --output-prefix junction_peptide_analysis
    """
}

process VALIDATE_MAXQUANT_SPLICE_JUNCTIONS {
    tag 'validated_splice_junctions'
    cpus 2; memory '12 GB'; time '24h'; disk '200 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/proteogenomics_validation", mode:'copy'
    input:
    path validation_stamp
    path candidates
    path splice_fastas
    path transcript_gtfs
    path reference_gtf
    path validation_script
    output:
    path 'validated_splice_junctions.detailed.tsv', emit: detailed
    path 'validated_splice_junctions.junction_spanning.tsv', emit: spanning
    path 'validated_splice_junctions.prioritized_novel_junctions.tsv', emit: prioritized
    path 'validated_splice_junctions.unresolved.tsv', emit: unresolved
    path 'validated_splice_junctions.summary.txt', emit: summary
    script:
    """
    set -euo pipefail
    python3 ${validation_script} \\
        --candidates ${candidates} \\
        --splice-fasta ${splice_fastas} \\
        --transcript-gtf ${transcript_gtfs} \\
        --reference-gtf ${reference_gtf} \\
        --output-prefix validated_splice_junctions
    """
}

process BUILD_PROTEOGENOMICS_EVIDENCE_REPORT {
    tag 'proteogenomics_evidence_report'
    cpus 1; memory '4 GB'; time '12h'; disk '200 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/proteogenomics_validation", mode:'copy'
    input:
    path validation_stamp
    path samplesheet
    path mqpar
    path evidence
    path msms
    path protein_groups
    path vep_vcfs
    path variant_annotation
    path peptide_mapping
    path splice_validation
    path searched_fastas
    val raw_map_mode
    path raw_file_map
    path report_script
    output:
    path 'proteogenomics_evidence.variants.tsv', emit: variants
    path 'proteogenomics_evidence.junctions.tsv', emit: junctions
    path 'proteogenomics_evidence.direct_msms_variants.tsv', emit: direct_variants
    path 'proteogenomics_evidence.mbr_only_variants.tsv', emit: mbr_variants
    path 'proteogenomics_evidence.direct_msms_junctions.tsv', emit: direct_junctions
    path 'proteogenomics_evidence.mbr_only_junctions.tsv', emit: mbr_junctions
    path 'proteogenomics_evidence.raw_file_mapping.tsv', emit: raw_mapping
    path 'proteogenomics_evidence.audit.tsv', emit: audit
    path 'proteogenomics_evidence.rejected_associations.tsv', emit: rejected
    path 'proteogenomics_evidence.validation_failures.md', emit: failure_report
    path 'proteogenomics_evidence.unique_variants.tsv', emit: unique_variants
    path 'proteogenomics_evidence.unique_junctions.tsv', emit: unique_junctions
    path 'proteogenomics_evidence.unique_direct_msms_variants.tsv', emit: unique_direct_variants
    path 'proteogenomics_evidence.unique_mbr_only_variants.tsv', emit: unique_mbr_variants
    path 'proteogenomics_evidence.sample_matched_direct_msms_variants.tsv', emit: sample_matched_direct_variants
    path 'proteogenomics_evidence.cross_sample_direct_msms_variants.tsv', emit: cross_sample_direct_variants
    path 'proteogenomics_evidence.sample_matched_mbr_only_variants.tsv', emit: sample_matched_mbr_variants
    path 'proteogenomics_evidence.cross_sample_mbr_only_variants.tsv', emit: cross_sample_mbr_variants
    path 'proteogenomics_evidence.sample_matched_direct_msms_junctions.tsv', emit: sample_matched_direct_junctions
    path 'proteogenomics_evidence.cross_sample_direct_msms_junctions.tsv', emit: cross_sample_direct_junctions
    path 'proteogenomics_evidence.sample_matched_mbr_only_junctions.tsv', emit: sample_matched_mbr_junctions
    path 'proteogenomics_evidence.cross_sample_mbr_only_junctions.tsv', emit: cross_sample_mbr_junctions
    path 'proteogenomics_evidence.evidence_classification.md', emit: classification_report
    path 'proteogenomics_evidence.report.md', emit: report
    path 'proteogenomics_evidence.summary.txt', emit: summary
    script:
    def rawMapArg = raw_map_mode == 'explicit' ? "--raw-file-map ${raw_file_map}" : ""
    """
    set -euo pipefail
    python ${report_script} \\
        --samples ${samplesheet} \\
        --mqpar ${mqpar} \\
        --evidence ${evidence} \\
        --msms ${msms} \\
        --protein-groups ${protein_groups} \\
        --vep-vcf ${vep_vcfs} \\
        --variant-annotation ${variant_annotation} \\
        --peptide-mapping ${peptide_mapping} \\
        --splice-validation ${splice_validation} \\
        --searched-fasta ${searched_fastas} \\
        ${rawMapArg} \\
        --output-prefix proteogenomics_evidence
    """
}


process BUILD_INTEGRATED_VARIANT_EVIDENCE {
    tag 'strict_integrated_variant_evidence'
    cpus 1; memory '4 GB'; time '8h'; disk '100 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/proteogenomics_validation", mode:'copy'
    input:
    path variants
    path codon_validation
    path codon_mismatch_analysis
    path integration_script
    output:
    path 'integrated_variant_evidence.all.tsv', emit: all
    path 'integrated_variant_evidence.strict.tsv', emit: strict
    path 'integrated_variant_evidence.excluded.tsv', emit: excluded
    path 'integrated_variant_evidence.report.md', emit: report
    script:
    """
    set -euo pipefail
    python ${integration_script} \
        --variants ${variants} \
        --codon-validation ${codon_validation} \
        --codon-mismatch-analysis ${codon_mismatch_analysis} \
        --output-prefix integrated_variant_evidence
    """
}

process VALIDATE_PROTEOGENOMIC_READS {
    tag 'proteogenomic_read_validation'
    cpus 1; memory '32 GB'; time '48h'; disk '600 GB'
    container "${params.pysam_image}"
    publishDir "${params.outdir}/proteogenomics_validation/read_validation", mode:'copy'
    input:
    path variants
    path junctions
    path splice_detail
    path arriba_tables
    path sorted_bams
    path reference_gtf
    path genome
    path validation_script
    output:
    path 'proteogenomic_read_validation.events.tsv'
    path 'proteogenomic_read_validation.reads.tsv'
    path 'proteogenomic_read_validation.fusions.tsv'
    path 'proteogenomic_read_validation.summary.txt', emit: summary
    path 'proteogenomic_read_validation.report.md', emit: report
    path 'proteogenomic_read_validation.variants.bed'
    path 'proteogenomic_read_validation.fusions.bedpe'
    path 'proteogenomic_read_validation.junctions.bed'
    path 'proteogenomic_read_validation.extraction_regions.bed'
    path 'proteogenomic_read_validation.igv.batch.txt'
    path 'proteogenomic_read_validation.*.bam*'
    script:
    def bamArgs = sorted_bams.findAll { file -> file.name.endsWith('.bam') }.collect { bam -> "--bam ${bam.baseName.tokenize('.')[0]}=${bam}" }.join(' ')
    """
    set -euo pipefail
    python3 ${validation_script} \
        --variants ${variants} \
        --junctions ${junctions} \
        --splice-detail ${splice_detail} \
        --arriba ${arriba_tables} \
        ${bamArgs} \
        --gtf ${reference_gtf} \
        --genome ${genome} \
        --padding ${params.read_validation_padding} \
        --output-prefix proteogenomic_read_validation
    """
}

process PREPARE_FINAL_MULTIQC_CONTENT {
    tag 'compact_integrated_report_content'
    cpus 1; memory '2 GB'; time '2h'; disk '10 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    input:
    path variant_inventory
    path rna_inventory
    path progression_summary
    path external_comparison
    path expression_ora
    path ranked_go
    path variant_set_go
    path proteogenomics_summary
    path integrated_report
    path evidence_classification
    path read_summary
    path codon_summary
    path provenance_summary
    path compact_builder
    val maxquant_enabled
    output:
    path 'multiqc_custom_content', emit: content
    script:
    """
    set -euo pipefail
    python3 ${compact_builder} \
      --output-dir multiqc_custom_content \
      --variant-inventory ${variant_inventory} \
      --rna-inventory ${rna_inventory} \
      --progression-summary ${progression_summary} \
      --external-comparison ${external_comparison} \
      --expression-ora ${expression_ora} \
      --ranked-go ${ranked_go} \
      --variant-set-go ${variant_set_go} \
      --proteogenomics-summary ${proteogenomics_summary} \
      --integrated-report ${integrated_report} \
      --evidence-classification ${evidence_classification} \
      --read-summary ${read_summary} \
      --codon-summary ${codon_summary} \
      --provenance-summary ${provenance_summary} \
      --maxquant-enabled ${maxquant_enabled}
    """
}

process MULTIQC_FINAL {
    tag 'integrated_final_report'
    cpus 4; memory '16 GB'; time '12h'; disk '150 GB'
    container 'quay.io/biocontainers/multiqc:1.35--pyhdfd78af_1'
    publishDir "${params.outdir}/multiqc", mode:'copy'
    input:
    path qc_files
    path custom_content
    path multiqc_config
    output:
    path 'multiqc_report.html', emit: report
    path 'multiqc_report_data', emit: data
    script:
    """
    set -euo pipefail
    multiqc . \\
        --force \\
        --title 'PGTK complete RNA-seq and proteogenomics report' \\
        --filename multiqc_report.html \\
        --config ${multiqc_config} \\
        --data-dir \\
        --data-format tsv
    rm -f multiqc_report_data/llms-full.txt
    test "\$(stat -c %s multiqc_report.html)" -le 25000000 || { echo "ERROR: final MultiQC HTML exceeds 25 MB" >&2; exit 1; }
    """
}

workflow {
    if (!params.sra_dir) error '--sra_dir is required'
    if (!params.reference_downloads) error '--reference_downloads is required'
    if (!params.container_cache) error '--container_cache is required'
    if (!params.pysam_image) error '--pysam_image is required'
    if (!params.ensembl_pep) error '--ensembl_pep is required'
    if (!(params.finding_priority_mode in ['all','filter'])) error '--finding_priority_mode must be all or filter'
    if ((params.igv_report_timeout_seconds as int) <= 0) error '--igv_report_timeout_seconds must be a positive integer'
    if ((params.igv_report_max_reads as int) <= 0) error '--igv_report_max_reads must be positive'
    if ((params.igv_report_max_file_size_mb as int) < 0) error '--igv_report_max_file_size_mb must be zero or positive'
    if ((params.finding_review_mapq as int) < 0) error '--finding_review_mapq must be non-negative'
    if ((params.finding_review_baseq as int) < 0) error '--finding_review_baseq must be non-negative'
    if ((params.finding_review_reference_reads as int) < 0) error '--finding_review_reference_reads must be non-negative'
    if ((params.read_validation_padding as int) < 0) error '--read_validation_padding must be non-negative'
    samples = channel.fromPath(params.samplesheet, checkIfExists:true).splitCsv(header:true).map { row ->
        if (!row.sample || !row.srr) error 'samples.csv requires sample and srr; TK, Group and baseline are optional'
        def sample = row.sample.trim()
        def srr = row.srr.trim()
        def baselineValue = (row.baseline ?: 'false').trim().toLowerCase()
        if (!(baselineValue in ['true','false'])) error "baseline must be true or false for ${sample}"
        def meta=[sample:sample, srr:srr, tk:(row.TK ?: sample).trim(), group:(row.Group ?: sample).trim(), baseline:baselineValue]
        def sraFile = file("${params.sra_dir}/${srr}/${srr}.sra", checkIfExists: true)
        tuple(meta, srr, sraFile)
    }
    reference_genome_archive = file("${params.reference_downloads}/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz", checkIfExists:true)
    reference_gtf_archive = file("${params.reference_downloads}/Homo_sapiens.GRCh38.111.gtf.gz", checkIfExists:true)
    reference_cdna_archive = file("${params.reference_downloads}/Homo_sapiens.GRCh38.cdna.all.fa.gz", checkIfExists:true)
    reference_proteome_archive = file("${params.reference_downloads}/human_reviewed_isoforms.fasta.gz", checkIfExists:true)
    reference_vep_archive = file("${params.reference_downloads}/homo_sapiens_vep_111_GRCh38.tar.gz", checkIfExists:true)
    reference_arriba_archive = file("${params.reference_downloads}/arriba_v2.4.0.tar.gz", checkIfExists:true)
    refs=DOWNLOAD_REFERENCES(
        reference_genome_archive,
        reference_gtf_archive,
        reference_cdna_archive,
        reference_proteome_archive,
        reference_vep_archive,
        reference_arriba_archive
    )
    ref=REF_INDEX(refs.genome)
    staridx=STAR_INDEX(refs.genome,refs.gtf)
    downloaded=SRA_TO_FASTQ(samples)
    reads=downloaded.map { m,r1,r2 -> tuple(m.sample,m,r1,r2) }.groupTuple(by:0).map { id,ms,r1s,r2s -> tuple(ms[0],r1s,r2s) } | CAT_FASTQ
    raw_qc=FASTQC_RAW(reads)
    trim_result=TRIM_GALORE(reads)
    trimmed=trim_result.reads
    trimmed_qc=FASTQC_TRIMMED(trimmed)
    star_result=STAR_ALIGN(trimmed,staridx)
    arriba_result=ARRIBA(star_result.bam,refs.genome,refs.gtf,refs.blacklist,refs.known,refs.domains)
    sortedbam=SORT_INDEX_BAM(star_result.bam)
    expression_analysis_script=file("${projectDir}/expression_go_analysis.py", checkIfExists:true)
    gene_count_parts=COUNT_GENES_PER_SAMPLE(sortedbam,refs.gtf)
    gene_expression=MERGE_GENE_EXPRESSION(gene_count_parts.map { m,c,su -> c }.collect(),refs.gtf,expression_analysis_script)
    assembled=STRINGTIE_ASSEMBLY(sortedbam,refs.gtf)
    novel_result=GFFCOMPARE_NOVEL(assembled,refs.gtf)
    rna_validator = file("${projectDir}/validate_rna_events.py", checkIfExists:true)
    validated_fusions=VALIDATE_RNA_FUSIONS(arriba_result.accepted,rna_validator)
    fusion_fasta=FUSION_FASTA(validated_fusions.validated)
    novel_keyed=novel_result.novel.map { m,g -> tuple(m.sample,m,g) }
    bam_keyed=sortedbam.map { m,b,bai -> tuple(m.sample,b,bai) }
    splice_validation_inputs=novel_keyed.join(bam_keyed).map { sample,m,g,b,bai -> tuple(m,g,b,bai) }
    validated_splice=VALIDATE_RNA_SPLICE_TRANSCRIPTS(splice_validation_inputs,rna_validator)
    splice_fasta=SPLICE_PROTEIN_FASTA(validated_splice.validated,refs.genome)
    flagstat=SAMTOOLS_FLAGSTAT(sortedbam)
    md_result=MARK_DUPLICATES(sortedbam)
    hc_intervals=PREPARE_HAPLOTYPE_INTERVALS(ref)
    split=SPLIT_N_CIGAR(md_result.bam,ref)
    hc_inputs=split.combine(hc_intervals.intervals.flatten()).map { m,b,bai,interval -> tuple(m,b,bai,interval) }
    hc_parts=HAPLOTYPE_CALLER(hc_inputs,ref)
    hc_grouped=hc_parts.map { m,shard,g,t -> tuple(m.sample,m,shard,g,t) }
        .groupTuple(by:0)
        .map { sample,metas,shards,gvcfs,tbis ->
            def ordered=[shards,gvcfs,tbis].transpose().sort { left,right -> left[0] <=> right[0] }
            tuple(metas[0],ordered.collect { it[1] },ordered.collect { it[2] })
        }
    shard_validator=file("${projectDir}/validate_haplotype_shards.py", checkIfExists:true)
    validated_hc_shards=VALIDATE_HAPLOTYPE_SHARDS(hc_grouped,shard_validator)
    gvcf=GATHER_HAPLOTYPE_GVCF(validated_hc_shards)
    raw_variants=GENOTYPE_VARIANTS(gvcf,ref)
    normalized_variants=NORMALIZE_VARIANTS(raw_variants,ref)
    selected_snps=SELECT_SNPS(normalized_variants,ref)
    selected_indels=SELECT_INDELS(normalized_variants,ref)
    snp_filtered=FILTER_SNPS(selected_snps,ref)
    indel_filtered=FILTER_INDELS(selected_indels,ref)
    separated_filter_inputs=snp_filtered.map { m,sf,sfi,sp,spi -> tuple(m.sample,m,sf,sfi,sp,spi) }
        .join(indel_filtered.map { m,inf,infi,inp,inpi -> tuple(m.sample,inf,infi,inp,inpi) })
        .map { sample,m,sf,sfi,sp,spi,inf,infi,inp,inpi -> tuple(m,sf,sfi,sp,spi,inf,infi,inp,inpi) }
    genotype_result=MERGE_FILTERED_VARIANTS(separated_filter_inputs)
    pass=genotype_result.pass
    variant_stats=BCFTOOLS_STATS(pass)
    annotated=VEP_ANNOTATE(pass,ref,refs.vep_cache)
    validated_variants=VALIDATE_RNA_VARIANTS(annotated,refs.genome,rna_validator)
    stage_qc_inputs=raw_variants.map { m,v,tbi -> tuple(m.sample,m,v,tbi) }
        .join(genotype_result.pass.map { m,v,tbi -> tuple(m.sample,v,tbi) })
        .join(validated_variants.validated.map { m,v,tbi -> tuple(m.sample,v,tbi) })
        .map { sample,m,raw,raw_tbi,pass_vcf,pass_tbi,rna,rna_tbi -> tuple(m,raw,raw_tbi,pass_vcf,pass_tbi,rna,rna_tbi) }
    stage_summarizer=file("${projectDir}/summarize_variant_stages.py", checkIfExists:true)
    variant_stage_qc=VARIANT_STAGE_QC(stage_qc_inputs,refs.genome,stage_summarizer)
    external_comparison_tables = channel.value([file("${projectDir}/external_comparison.none.tsv", checkIfExists:true)])
    if (strictBooleanParam(params.run_external_vcf_comparison, '--run_external_vcf_comparison')) {
        if (!params.external_vcf_dir) error '--external_vcf_dir is required when --run_external_vcf_comparison true'
        comparison_script=file("${projectDir}/compare_external_vcf.py", checkIfExists:true)
        external_raw_inputs=raw_variants.map { m,v,tbi ->
            def comparisonMeta = m + [comparison_stage:'raw']
            tuple(comparisonMeta,v,tbi,resolveExternalVcf(params.external_vcf_dir.toString(),m.srr,params.external_vcf_suffix.toString()))
        }
        external_pass_inputs=genotype_result.pass.map { m,v,tbi ->
            def comparisonMeta = m + [comparison_stage:'pass']
            tuple(comparisonMeta,v,tbi,resolveExternalVcf(params.external_vcf_dir.toString(),m.srr,params.external_vcf_suffix.toString()))
        }
        external_rna_inputs=validated_variants.validated.map { m,v,tbi ->
            def comparisonMeta = m + [comparison_stage:'rna_validated']
            tuple(comparisonMeta,v,tbi,resolveExternalVcf(params.external_vcf_dir.toString(),m.srr,params.external_vcf_suffix.toString()))
        }
        external_inputs=external_raw_inputs.mix(external_pass_inputs,external_rna_inputs)
        external_comparison_result=COMPARE_EXTERNAL_VCF(external_inputs,comparison_script)
        external_comparison_tables = external_comparison_result.summary.collect()
    }
    variant_read_provenance_script = file("${projectDir}/validate_variant_read_provenance.py", checkIfExists:true)
    variant_codon_script = file("${projectDir}/validate_variant_codons.py", checkIfExists:true)
    variant_validation_merger = file("${projectDir}/merge_variant_validation.py", checkIfExists:true)
    validated_variant_keyed = validated_variants.validated.map { m,v,t -> tuple(m.sample,m,v,t) }
    sorted_bam_keyed_for_variant_validation = sortedbam.map { m,b,bai -> tuple(m.sample,b,bai) }
    variant_validation_inputs = validated_variant_keyed.join(sorted_bam_keyed_for_variant_validation).map { sample,m,v,t,b,bai -> tuple(m,v,t,b,bai) }
    variant_codon_parts = VALIDATE_VARIANT_CODONS(variant_validation_inputs, refs.genome, variant_codon_script)
    variant_provenance_parts = VALIDATE_VARIANT_READ_PROVENANCE(variant_validation_inputs, file(params.samplesheet, checkIfExists:true), variant_read_provenance_script)
    variant_codon_validation = MERGE_VARIANT_CODON_VALIDATION(variant_codon_parts.map { m,a,v,p,f,c,s,r -> [a,v,p,f,c,s,r] }.flatten().collect(), variant_validation_merger)
    codon_mismatch_script = file("${projectDir}/analyze_codon_mismatches.py", checkIfExists:true)
    codon_mismatch_analysis = ANALYZE_CODON_MISMATCHES(variant_codon_validation.all, codon_mismatch_script)
    variant_read_provenance = MERGE_VARIANT_READ_PROVENANCE(variant_provenance_parts.map { m,t,s,r -> [t,s,r] }.flatten().collect(), variant_validation_merger)
    variant_fasta=PYPGATK_FASTA(validated_variants.validated,refs.gtf,refs.cdna)

    variant_keyed=variant_fasta.map { m,f -> tuple(m.sample,m,f) }
    fusion_keyed=fusion_fasta.map { m,f -> tuple(m.sample,m,f) }
    splice_keyed=splice_fasta.map { m,f -> tuple(m.sample,m,f) }
    combined_inputs=variant_keyed
        .join(fusion_keyed)
        .join(splice_keyed)
        .map { sample,m1,vf,m2,ff,m3,sf -> tuple(m1,vf,ff,sf) }
    combined_fasta=COMBINE_PROTEIN_FASTA(combined_inputs)
    groups=validated_variants.validated.branch { m,v,t -> baseline:m.baseline=='true'; progression:m.baseline=='false'; other:true }
    bases=groups.baseline.map { m,v,t -> tuple(m.tk,v,t) }
    pairs=groups.progression.map { m,v,t -> tuple(m.tk,m,v,t) }.combine(bases,by:0).map { k,m,pv,pt,bv,bt -> tuple(m,pv,pt,bv,bt) }
    prog=PROGRESSION_SUBTRACT(pairs)
    if (!params.go_obo) error '--go_obo is required for unbiased GO enrichment'
    if (!params.go_gaf) error '--go_gaf is required for unbiased GO enrichment'
    go_preparation_script=file("${projectDir}/prepare_go_annotations.py", checkIfExists:true)
    progression_biology_script=file("${projectDir}/analyze_progression_biology.py", checkIfExists:true)
    progression_pair_script=file("${projectDir}/compare_progression_pair.py", checkIfExists:true)
    progression_merge_script=file("${projectDir}/merge_progression_biology.py", checkIfExists:true)
    go_reference=PREPARE_GO_ANNOTATIONS(file(params.go_obo,checkIfExists:true),file(params.go_gaf,checkIfExists:true),go_preparation_script)
    variant_landscape_script=file("${projectDir}/analyze_variant_landscape.py",checkIfExists:true)
    variant_landscape=ANALYZE_VARIANT_LANDSCAPE(
        raw_variants.map { m,v,t -> v }.collect(), normalized_variants.map { m,v,t -> v }.collect(),
        genotype_result.filtered.map { m,v,t -> v }.collect(), genotype_result.pass.map { m,v,t -> v }.collect(),
        annotated.map { m,v,t -> v }.collect(), validated_variants.validated.map { m,v,t -> v }.collect(),
        prog.map { m,nv,nt,bv,bt,sv,st,su -> nv }.collect(), prog.map { m,nv,nt,bv,bt,sv,st,su -> bv }.collect(),
        prog.map { m,nv,nt,bv,bt,sv,st,su -> sv }.collect(), go_reference.mapping, variant_landscape_script
    )
    progression_sample_inputs=prog.map { m,nv,nt,bv,bt,sv,st,su -> tuple(m.sample,m,nv,nt) }.join(validated_variants.validated.map { m,v,t -> tuple(m.sample,v,t) }).map { sample,m,nv,nt,rv,rt -> tuple(m,nv,nt,rv,rt) }
    progression_sample_results=ANALYZE_PROGRESSION_SAMPLE(progression_sample_inputs,go_reference.mapping,progression_biology_script)
    progression_variant_sets=ANALYZE_PROGRESSION_VARIANT_SETS(
        progression_sample_results.map { m,a,g,e,c,su -> g }.collect(),
        file(params.samplesheet,checkIfExists:true), go_reference.mapping, refs.gtf, expression_analysis_script
    )
    expression_multiqc_content = channel.empty()
    if (strictBooleanParam(params.run_expression_go, '--run_expression_go')) {
        expression_metadata = channel.fromPath(params.samplesheet, checkIfExists:true)
            .splitCsv(header:true)
            .map { row ->
                def sample = row.sample.trim()
                def baselineValue = (row.baseline ?: 'false').trim().toLowerCase()
                [sample:sample, tk:(row.TK ?: sample).trim(), group:(row.Group ?: sample).trim(), baseline:baselineValue]
            }
        expression_all_samples = channel.fromPath(params.samplesheet, checkIfExists:true)
            .splitCsv(header:true)
            .map { row -> row.sample.trim() }
            .collect()
            .map { names -> names.join(',') }
        expression_ora_inputs = expression_metadata
        expression_subject_groups = channel.fromPath(params.samplesheet, checkIfExists:true)
            .splitCsv(header:true)
            .map { row ->
                def sample = row.sample.trim()
                def baselineValue = (row.baseline ?: 'false').trim().toLowerCase()
                def meta = [sample:sample, tk:(row.TK ?: sample).trim(), group:(row.Group ?: sample).trim(), baseline:baselineValue]
                tuple(meta.tk, meta)
            }
            .groupTuple(by:0)
        expression_rank_inputs = expression_subject_groups.flatMap { subject, members ->
            def baselines = members.findAll { it.baseline == 'true' }
            def progressions = members.findAll { it.baseline == 'false' }
            if (baselines.size() != 1) return []
            progressions.collect { meta -> tuple(meta, baselines[0].sample) }
        }
        expression_sample_go = ANALYZE_EXPRESSION_SAMPLE_GO(
            expression_ora_inputs, gene_expression.matrix, expression_all_samples,
            go_reference.mapping, expression_analysis_script
        )
        expression_ranked_go = ANALYZE_EXPRESSION_RANKED_GO(
            expression_rank_inputs, gene_expression.matrix, expression_all_samples,
            go_reference.mapping, expression_analysis_script
        )
        expression_go = MERGE_EXPRESSION_GO(
            expression_sample_go.map { m,ora,summary -> ora }.collect(),
            expression_ranked_go.map { m,b,ranked,summary -> ranked }.collect(),
            expression_sample_go.map { m,ora,summary -> summary }
                .mix(expression_ranked_go.map { m,b,ranked,summary -> summary })
                .collect(),
            file(params.samplesheet,checkIfExists:true), expression_analysis_script
        )
        expression_multiqc_builder=file("${projectDir}/build_expression_multiqc_content.py",checkIfExists:true)
        expression_multiqc_content = PREPARE_EXPRESSION_MULTIQC_CONTENT(
            expression_go.ora, expression_go.ranked, expression_go.summary,
            progression_variant_sets.enrichment, progression_variant_sets.summary,
            expression_multiqc_builder
        )
    }
    progression_pair_left=progression_sample_results.map { m,a,g,e,c,su -> tuple(m.tk,m,a,g,e) }
    progression_pair_right=progression_sample_results.map { m,a,g,e,c,su -> tuple(m.tk,m,a,g,e) }
    progression_pair_inputs=progression_pair_left.combine(progression_pair_right,by:0).filter { subject,ma,aa,ga,ea,mb,ab,gb,eb -> ma.sample < mb.sample }.map { subject,ma,aa,ga,ea,mb,ab,gb,eb -> def pm=[subject:subject,sample_a:ma.sample,sample_b:mb.sample,pair_id:"${subject}.${ma.sample}_vs_${mb.sample}"]; tuple(pm,aa,ga,ea,ab,gb,eb) }
    progression_pair_results=COMPARE_PROGRESSION_PAIR(progression_pair_inputs,progression_pair_script)
    progression_biology=MERGE_PROGRESSION_BIOLOGY(
        progression_sample_results.map { m,a,g,e,c,su -> a }.collect(), progression_sample_results.map { m,a,g,e,c,su -> g }.collect(), progression_sample_results.map { m,a,g,e,c,su -> e }.collect(), progression_sample_results.map { m,a,g,e,c,su -> c }.collect(), progression_sample_results.map { m,a,g,e,c,su -> su }.collect(),
        progression_pair_results.map { m,a,g,e,su -> a }.collect(), progression_pair_results.map { m,a,g,e,su -> g }.collect(), progression_pair_results.map { m,a,g,e,su -> e }.collect(), progression_pair_results.map { m,a,g,e,su -> su }.collect(),
        go_reference.metadata, progression_merge_script
    )
    qc_files = raw_qc.qc
        .mix(trimmed_qc.qc, trim_result.reports, star_result.logs,
             flagstat, md_result.metrics, variant_stats)
        .collect()
    multiqc_result=MULTIQC_QC_DATA(qc_files)
    igv_bundle_script=file("${projectDir}/build_igv_evidence_bundle.py",checkIfExists:true)
    igv_bundle=BUILD_IGV_EVIDENCE_BUNDLE(
        validated_variants.validated.map { m,v,t -> v }.collect(),
        prog.map { m,nv,nt,bv,bt,sv,st,su -> nv }.collect(),
        validated_fusions.validated.map { m,f -> f }.collect(),
        validated_splice.audit.map { m,f -> f }.collect(),
        sortedbam.map { m,b,bai -> [b,bai] }.flatten().collect(),
        refs.genome,
        igv_bundle_script
    )
    finding_review_script=file("${projectDir}/build_finding_igv_reviews.py",checkIfExists:true)
    finding_reviews=BUILD_FINDING_IGV_REVIEWS(
        igv_bundle.events,
        igv_bundle.bams.collect(),
        refs.genome,
        finding_review_script
    )
    if (strictBooleanParam(params.generate_priority_igv_reports, '--generate_priority_igv_reports')) {
        finding_explorer_script=file("${projectDir}/build_finding_explorer.py",checkIfExists:true)
        finding_explorer_launcher=file("${projectDir}/serve_finding_explorer.sh",checkIfExists:true)
        finding_explorer=BUILD_FINDING_EXPLORER(finding_reviews.reviews, refs.genome, finding_explorer_script, finding_explorer_launcher)
    }
    comparative_report_script=file("${projectDir}/build_comparative_advantage_report.py", checkIfExists:true)
    comparative_report=BUILD_COMPARATIVE_ADVANTAGE_REPORT(
        file(params.samplesheet,checkIfExists:true),
        raw_variants.map { m,v,t -> v }.collect(),
        genotype_result.pass.map { m,v,t -> v }.collect(),
        validated_variants.validated.map { m,v,t -> v }.collect(),
        prog.map { m,nv,nt,bv,bt,sv,st,su -> nv }.collect(),
        validated_fusions.validated.map { m,f -> f }.collect(),
        validated_splice.audit.map { m,f -> f }.collect(),
        variant_fasta.map { m,f -> f }.collect(),
        fusion_fasta.map { m,f -> f }.collect(),
        splice_fasta.map { m,f -> f }.collect(),
        combined_fasta.map { m,f -> f }.collect(),
        variant_stage_qc.table.collect(),
        external_comparison_tables,
        variant_codon_validation.summary,
        variant_read_provenance.summary,
        comparative_report_script
    )


    complete_reporter=file("${projectDir}/build_complete_report.py", checkIfExists:true)
    report_samplesheet=file(params.samplesheet, checkIfExists:true)
    complete_findings=BUILD_COMPLETE_FINDINGS_REPORT(
        report_samplesheet,
        multiqc_result.data,
        validated_variants.audit.map { m,f -> f }.collect(),
        validated_fusions.audit.map { m,f -> f }.collect(),
        validated_splice.audit.map { m,f -> f }.collect(),
        annotated.map { m,v,t -> v }.collect(),
        arriba_result.accepted.map { m,f -> f }.collect(),
        arriba_result.discarded.collect(),
        assembled.map { m,g -> g }.collect(),
        novel_result.annotated.map { m,g -> g }.collect(),
        novel_result.novel.map { m,g -> g }.collect(),
        variant_fasta.map { m,f -> f }.collect(),
        fusion_fasta.map { m,f -> f }.collect(),
        splice_fasta.map { m,f -> f }.collect(),
        prog.map { m,nv,nt,bv,bt,sv,st,su -> nv }.collect(),
        complete_reporter
    )

    multiqc_report_builder=file("${projectDir}/build_pgtk_multiqc_content.py",checkIfExists:true)
    comparative_multiqc=PREPARE_COMPARATIVE_MULTIQC_CONTENT(
        file(params.samplesheet,checkIfExists:true), comparative_report.report,
        comparative_report.variant_inventory, comparative_report.fasta_inventory,
        comparative_report.rna_inventory, comparative_report.external_comparison,
        comparative_report.multiqc_summary, progression_biology.report,
        progression_biology.multiqc_summary, progression_biology.enrichment,
        progression_biology.pairwise_categories, complete_findings.report,
        complete_findings.failure_report, complete_findings.variant_explanations,
        multiqc_report_builder, booleanText(params.run_proteogenomic_validation, '--run_proteogenomic_validation')
    )

    multiqc_config_file=file("${projectDir}/multiqc_config.yaml", checkIfExists:true)
    if (strictBooleanParam(params.run_proteogenomic_validation, '--run_proteogenomic_validation')) {
        if (!new File(params.maxquant_txt.toString()).isDirectory()) error "MaxQuant folder not found: ${params.maxquant_txt}"
        mq_peptides = file("${params.maxquant_txt}/peptides.txt", checkIfExists:true)
        mq_evidence = file("${params.maxquant_txt}/evidence.txt", checkIfExists:true)
        mq_msms = file("${params.maxquant_txt}/msms.txt", checkIfExists:true)
        mq_protein_groups = file("${params.maxquant_txt}/proteinGroups.txt", checkIfExists:true)
        mq_mqpar = resolveMaxQuantMqpar(params.maxquant_txt, params.maxquant_mqpar)
        mq_contaminants = resolveContaminants(params.maxquant_contaminants?.toString(), params.maxquant_txt.toString())
        raw_map_mode = params.maxquant_raw_map ? 'explicit' : 'default'
        mq_raw_map = params.maxquant_raw_map ? file(params.maxquant_raw_map, checkIfExists:true) : file("${projectDir}/maxquant_raw_file_map.none.tsv", checkIfExists:true)
        ensembl_pep = file(params.ensembl_pep, checkIfExists:true)

        mapper_script = file("${projectDir}/map_peptides_to_fasta.py", checkIfExists:true)
        annotation_script = file("${projectDir}/annotate_variant_peptides.py", checkIfExists:true)
        junction_script = file("${projectDir}/analyze_chimeric_splice_peptides.py", checkIfExists:true)
        splice_validation_script = file("${projectDir}/validate_splice_junction_peptides.py", checkIfExists:true)
        report_script = file("${projectDir}/proteogenomics_evidence_report.py", checkIfExists:true)
        read_validation_script = file("${projectDir}/validate_proteogenomic_reads.py", checkIfExists:true)
        report_samplesheet = file(params.samplesheet, checkIfExists:true)

        combined_fastas_for_validation = combined_fasta.map { m,f -> f }.collect()
        sample_fasta_basenames = new File(params.samplesheet.toString()).readLines().drop(1).findAll { it.trim() }.collect { row ->
            def sample = row.split(',', -1)[0].trim()
            "${sample}.exploratory_proteogenomics.fasta"
        }
        mq_canonicals = resolveMaxQuantCanonicalFastas(mq_mqpar, params.maxquant_txt, params.maxquant_canonical_fasta, sample_fasta_basenames)
        validation_stamp = VALIDATE_MAXQUANT_INPUTS(mq_peptides, mq_evidence, mq_msms, mq_protein_groups, mq_mqpar, mq_canonicals, mq_contaminants)
        fusion_fastas_for_validation = fusion_fasta.map { m,f -> f }.collect()
        splice_fastas_for_validation = splice_fasta.map { m,f -> f }.collect()
        assembled_gtfs_for_validation = validated_splice.validated.map { m,g -> g }.collect()
        arriba_tables_for_validation = validated_fusions.validated.map { m,f -> f }.collect()
        vep_vcfs_for_validation = validated_variants.validated.map { m,v,t -> v }.collect()
        searched_fastas_for_report = combined_fastas_for_validation.map { fastas -> fastas + mq_canonicals }

        peptide_mapping = MAP_MAXQUANT_PEPTIDES(validation_stamp, mq_peptides, mq_canonicals, mq_contaminants, combined_fastas_for_validation, mapper_script)
        variant_annotation = ANNOTATE_MAXQUANT_VARIANTS(validation_stamp, peptide_mapping.mapping, vep_vcfs_for_validation, combined_fastas_for_validation, ensembl_pep, annotation_script)
        junction_analysis = ANALYZE_MAXQUANT_JUNCTIONS(validation_stamp, mq_peptides, mq_canonicals, mq_contaminants, fusion_fastas_for_validation, splice_fastas_for_validation, arriba_tables_for_validation, junction_script)
        splice_validation = VALIDATE_MAXQUANT_SPLICE_JUNCTIONS(validation_stamp, junction_analysis.splice_candidates, splice_fastas_for_validation, assembled_gtfs_for_validation, refs.gtf, splice_validation_script)
        evidence_report = BUILD_PROTEOGENOMICS_EVIDENCE_REPORT(validation_stamp, report_samplesheet, mq_mqpar, mq_evidence, mq_msms, mq_protein_groups, vep_vcfs_for_validation, variant_annotation.detailed, peptide_mapping.mapping, splice_validation.detailed, searched_fastas_for_report, raw_map_mode, mq_raw_map, report_script)
        integrated_variant_script = file("${projectDir}/build_integrated_variant_evidence.py", checkIfExists:true)
        integrated_variant_evidence = BUILD_INTEGRATED_VARIANT_EVIDENCE(evidence_report.variants, variant_codon_validation.all, codon_mismatch_analysis.detailed, integrated_variant_script)
        sorted_bams_for_validation = sortedbam.map { m,b,bai -> [b,bai] }.flatten().collect()
        read_validation = VALIDATE_PROTEOGENOMIC_READS(evidence_report.variants, evidence_report.junctions, splice_validation.detailed, arriba_tables_for_validation, sorted_bams_for_validation, refs.gtf, refs.genome, read_validation_script)
        validation_semantics_doc = file("${projectDir}/PIPELINE_VALIDATION_SEMANTICS.md", checkIfExists:true)
        final_multiqc_content = PREPARE_FINAL_MULTIQC_CONTENT(
            comparative_report.variant_inventory,
            comparative_report.rna_inventory,
            progression_biology.multiqc_summary,
            comparative_report.external_comparison,
            expression_go.ora,
            expression_go.ranked,
            progression_variant_sets.enrichment,
            evidence_report.summary,
            integrated_variant_evidence.report,
            evidence_report.classification_report,
            read_validation.summary,
            variant_codon_validation.summary,
            variant_read_provenance.summary,
            file("${projectDir}/build_compact_multiqc_content.py", checkIfExists:true),
            booleanText(params.run_proteogenomic_validation, '--run_proteogenomic_validation')
        )
        final_multiqc_inputs = final_multiqc_content.content.mix(variant_landscape.multiqc).collect()
        MULTIQC_FINAL(qc_files, final_multiqc_inputs, multiqc_config_file)
    } else {
        final_multiqc_inputs = comparative_multiqc.mix(expression_multiqc_content, variant_landscape.multiqc).collect()
        MULTIQC_FINAL(qc_files, final_multiqc_inputs, multiqc_config_file)
    }
}
