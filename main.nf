#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.samplesheet = "${projectDir}/samples.csv"
params.outdir = "${projectDir}/results"
params.read_length = 150
params.skip_trimming = false
params.sra_dir = "${projectDir}/sra_cache"
params.genome_url = 'https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.gtf_url = 'https://ftp.ensembl.org/pub/release-111/gtf/homo_sapiens/Homo_sapiens.GRCh38.111.gtf.gz'
params.vep_cache_url = 'https://ftp.ensembl.org/pub/release-111/variation/indexed_vep_cache/homo_sapiens_vep_111_GRCh38.tar.gz'
params.proteome_url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29'
params.arriba_url = 'https://github.com/suhrig/arriba/releases/download/v2.4.0/arriba_v2.4.0.tar.gz'

process DOWNLOAD_REFERENCES {
    tag 'GRCh38_Ensembl111'
    cpus 2; memory '8 GB'; time '12h'; disk '100 GB'; queue 'normal'
    publishDir "${params.outdir}/references", mode: 'copy'
    output:
    path 'refs/genome.fa', emit: genome
    path 'refs/genes.gtf', emit: gtf
    path 'refs/human_reviewed_isoforms.fasta', emit: proteome
    path 'refs/vep_cache', emit: vep_cache
    path 'refs/arriba_blacklist.tsv.gz', emit: blacklist
    path 'refs/arriba_known_fusions.tsv.gz', emit: known
    path 'refs/arriba_protein_domains.gff3', emit: domains
    script:
    """
    set -euo pipefail
    command -v curl >/dev/null; command -v gzip >/dev/null; command -v tar >/dev/null
    mkdir -p refs/vep_cache refs/arriba_unpack
    curl -fL --retry 5 '${params.genome_url}' | gzip -dc > refs/genome.fa
    curl -fL --retry 5 '${params.gtf_url}' | gzip -dc > refs/genes.gtf
    curl -fL --retry 5 '${params.proteome_url}' | gzip -dc > refs/human_reviewed_isoforms.fasta
    curl -fL --retry 5 '${params.vep_cache_url}' -o vep.tar.gz
    tar -xzf vep.tar.gz -C refs/vep_cache
    curl -fL --retry 5 '${params.arriba_url}' -o arriba.tar.gz
    tar -xzf arriba.tar.gz -C refs/arriba_unpack
    cp \$(find refs/arriba_unpack -type f -name 'blacklist_hg38_GRCh38*.tsv.gz' | head -1) refs/arriba_blacklist.tsv.gz
    cp \$(find refs/arriba_unpack -type f -name 'known_fusions_hg38_GRCh38*.tsv.gz' | head -1) refs/arriba_known_fusions.tsv.gz
    cp \$(find refs/arriba_unpack -type f -name 'protein_domains_hg38_GRCh38*.gff3' | head -1) refs/arriba_protein_domains.gff3
    test -s refs/genome.fa; test -s refs/genes.gtf; test -s refs/human_reviewed_isoforms.fasta
    test -d refs/vep_cache/homo_sapiens/111_GRCh38
    """
}

process SRA_TO_FASTQ {
    tag "${meta.sample}:${srr}"
    cpus 8; memory '16 GB'; time '24h'; disk '150 GB'; queue 'normal'
    container 'quay.io/biocontainers/sra-tools:3.2.1--h4304569_0'
    input: tuple val(meta), val(srr), path(sra_file)
    output: tuple val(meta), path("${srr}_1.fastq.gz"), path("${srr}_2.fastq.gz")
    script:
    """
    set -euo pipefail
    fasterq-dump --split-files --threads ${task.cpus} ${sra_file}
    gzip ${srr}_1.fastq ${srr}_2.fastq
    """
}

process CAT_FASTQ {
    tag "${meta.sample}"
    cpus 1; memory '2 GB'; time '4h'; disk '200 GB'; queue 'normal'
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

process TRIM_GALORE {
    tag "${meta.sample}"
    cpus 4; memory '8 GB'; time '12h'; disk '150 GB'; queue 'normal'
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

process STAR_INDEX {
    tag 'GRCh38_Ensembl111'
    cpus 32; memory '32 GB'; time '24h'; disk '120 GB'; queue 'bigmem'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'
    publishDir "${params.outdir}/references/star_index", mode:'copy'
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
    cpus 32; memory '32 GB'; time '24h'; disk '200 GB'; queue 'bigmem'
    container 'quay.io/biocontainers/star:2.7.11b--h43eeafb_1'
    input: tuple val(meta), path(r1), path(r2); path index
    output: tuple val(meta), path("${meta.sample}.Aligned.out.bam")
    script:
    """
    STAR \
        --genomeDir ${index} \
        --readFilesIn ${r1} ${r2} \
        --readFilesCommand zcat \
        --runThreadN ${task.cpus} \
        --twopassMode Basic \
        --outFileNamePrefix ${meta.sample}. \
        --outStd Log \
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
    cpus 8; memory '16 GB'; time '24h'; disk '200 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/bam/star", mode:'copy', pattern:'*.Aligned.sortedByCoord.out.bam*'
    input: tuple val(meta), path(bam)
    output: tuple val(meta), path("${meta.sample}.Aligned.sortedByCoord.out.bam"), path("${meta.sample}.Aligned.sortedByCoord.out.bam.bai")
    script:
    """
    samtools sort -@ ${task.cpus} -o ${meta.sample}.Aligned.sortedByCoord.out.bam ${bam}
    samtools index -@ ${task.cpus} ${meta.sample}.Aligned.sortedByCoord.out.bam
    """
}

process REF_INDEX {
    tag 'GRCh38'; cpus 2; memory '8 GB'; time '4h'; disk '30 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: path genome
    output: tuple path('genome.fa'), path('genome.fa.fai'), path('genome.dict')
    script:
    """
    cp ${genome} genome.fa
    samtools faidx genome.fa
    gatk CreateSequenceDictionary -R genome.fa -O genome.dict
    """
}

process MARK_DUPLICATES {
    tag "${meta.sample}"; cpus 4; memory '16 GB'; time '24h'; disk '200 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: tuple val(meta), path(bam), path(bai)
    output: tuple val(meta), path("${meta.sample}.markdup.bam"), path("${meta.sample}.markdup.bam.bai")
    script:
    """
    gatk MarkDuplicates -I ${bam} -O ${meta.sample}.markdup.bam -M ${meta.sample}.metrics.txt --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
    mv ${meta.sample}.markdup.bai ${meta.sample}.markdup.bam.bai
    """
}

process SPLIT_N_CIGAR {
    tag "${meta.sample}"; cpus 4; memory '16 GB'; time '24h'; disk '200 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    input: tuple val(meta), path(bam), path(bai); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.split.bam"), path("${meta.sample}.split.bam.bai")
    script:
    """
    gatk SplitNCigarReads -R ${genome} -I ${bam} -O ${meta.sample}.split.bam --create-output-bam-index true
    """
}

process HAPLOTYPE_CALLER {
    tag "${meta.sample}"; cpus 10; memory '16 GB'; time '48h'; disk '120 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/gvcf", mode:'copy'
    input: tuple val(meta), path(bam), path(bai); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.g.vcf.gz"), path("${meta.sample}.g.vcf.gz.tbi")
    script:
    """
    gatk HaplotypeCaller -R ${genome} -I ${bam} -O ${meta.sample}.g.vcf.gz -ERC GVCF --dont-use-soft-clipped-bases true --standard-min-confidence-threshold-for-calling 20 --native-pair-hmm-threads ${task.cpus}
    """
}

process GENOTYPE_FILTER {
    tag "${meta.sample}"; cpus 2; memory '12 GB'; time '12h'; disk '40 GB'; queue 'normal'
    container 'quay.io/biocontainers/gatk4:4.6.1.0--py310hdfd78af_0'
    publishDir "${params.outdir}/vcf_pass", mode:'copy', pattern:'*.pass.vcf.gz*'
    input: tuple val(meta), path(gvcf), path(tbi); tuple path(genome), path(fai), path(dict)
    output: tuple val(meta), path("${meta.sample}.pass.vcf.gz"), path("${meta.sample}.pass.vcf.gz.tbi")
    script:
    """
    gatk GenotypeGVCFs -R ${genome} -V ${gvcf} -O ${meta.sample}.raw.vcf.gz
    gatk VariantFiltration -R ${genome} -V ${meta.sample}.raw.vcf.gz --window 35 --cluster 3 --filter-expression 'QD < 2.0' --filter-name QD2 --filter-expression 'FS > 30.0' --filter-name FS30 --filter-expression 'MQ < 40.0' --filter-name MQ40 --filter-expression 'MQRankSum < -12.5' --filter-name MQRankSum-12.5 --filter-expression 'ReadPosRankSum < -8.0' --filter-name ReadPos-8 -O ${meta.sample}.filtered.vcf.gz
    gatk SelectVariants -R ${genome} -V ${meta.sample}.filtered.vcf.gz --exclude-filtered -O ${meta.sample}.pass.vcf.gz
    """
}

process VEP_ANNOTATE {
    tag "${meta.sample}"; cpus 10; memory '16 GB'; time '24h'; disk '60 GB'; queue 'normal'
    container 'quay.io/biocontainers/ensembl-vep:111.0--pl5321h2a3209d_0'
    publishDir "${params.outdir}/vep", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); tuple path(genome), path(fai), path(dict); path cache
    output: tuple val(meta), path("${meta.sample}.vep.vcf.gz"), path("${meta.sample}.vep.vcf.gz.tbi")
    script:
    """
    vep --input_file ${vcf} --output_file ${meta.sample}.vep.vcf --format vcf --vcf --cache --offline --cache_version 111 --dir_cache ${cache} --species homo_sapiens --assembly GRCh38 --fasta ${genome} --pick --canonical --protein --symbol --numbers --biotype --total_length --hgvs --fork ${task.cpus} --force_overwrite
    bgzip ${meta.sample}.vep.vcf && tabix -p vcf ${meta.sample}.vep.vcf.gz
    """
}

process PYPGATK_FASTA {
    tag "${meta.sample}"; cpus 2; memory '12 GB'; time '12h'; disk '40 GB'; queue 'normal'
    container 'quay.io/biocontainers/pypgatk:0.0.24--pyhdfd78af_0'
    publishDir "${params.outdir}/variant_fasta", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); path gtf; path proteome
    output: tuple val(meta), path("${meta.sample}.variant_proteins.fasta")
    script:
    """
    pypgatk vcf-to-proteindb --input-vcf ${vcf} --gene-annotations-gtf ${gtf} --protein-db-fasta ${proteome} --af-field AF --annotation-field-name CSQ --consequence-filter missense_variant,frameshift_variant,stop_gained,stop_lost,start_lost,splice_donor_variant,splice_acceptor_variant,inframe_insertion,inframe_deletion --output-proteindb ${meta.sample}.variant_proteins.fasta
    sed -i 's/^>/>${meta.sample}|/' ${meta.sample}.variant_proteins.fasta
    """
}

process ARRIBA {
    tag "${meta.sample}"; cpus 2; memory '12 GB'; time '12h'; disk '60 GB'; queue 'normal'
    container 'quay.io/biocontainers/arriba:2.4.0--h0033a41_2'
    publishDir "${params.outdir}/fusions", mode:'copy'
    input: tuple val(meta), path(bam); path genome; path gtf; path blacklist; path known; path domains
    output: tuple val(meta), path("${meta.sample}.fusions.tsv"); path "${meta.sample}.fusions.discarded.tsv"
    script:
    """
    arriba -x ${bam} -a ${genome} -g ${gtf} -b ${blacklist} -k ${known} -p ${domains} -o ${meta.sample}.fusions.tsv -O ${meta.sample}.fusions.discarded.tsv
    """
}

process PROGRESSION_SUBTRACT {
    tag "${meta.sample}"; cpus 1; memory '4 GB'; time '4h'; disk '20 GB'; queue 'normal'
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
    tag "${meta.sample}"; cpus 2; memory '12 GB'; time '12h'; disk '40 GB'; queue 'normal'
    container 'quay.io/biocontainers/pypgatk:0.0.24--pyhdfd78af_0'
    publishDir "${params.outdir}/progression_fasta", mode:'copy'
    input: tuple val(meta), path(vcf), path(tbi); path gtf; path proteome
    output: path "${meta.sample}.progression_proteins.fasta"
    script:
    """
    pypgatk vcf-to-proteindb --input-vcf ${vcf} --gene-annotations-gtf ${gtf} --protein-db-fasta ${proteome} --af-field AF --annotation-field-name CSQ --consequence-filter missense_variant,frameshift_variant,stop_gained,stop_lost,start_lost,splice_donor_variant,splice_acceptor_variant,inframe_insertion,inframe_deletion --output-proteindb ${meta.sample}.progression_proteins.fasta
    sed -i 's/^>/>${meta.sample}|PROGRESSION|/' ${meta.sample}.progression_proteins.fasta
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
    trimmed=params.skip_trimming ? reads : TRIM_GALORE(reads).reads
    starbam=STAR_ALIGN(trimmed,staridx)
    ARRIBA(starbam,refs.genome,refs.gtf,refs.blacklist,refs.known,refs.domains)
    sortedbam=SORT_INDEX_BAM(starbam)
    md=MARK_DUPLICATES(sortedbam)
    split=SPLIT_N_CIGAR(md,ref)
    gvcf=HAPLOTYPE_CALLER(split,ref)
    pass=GENOTYPE_FILTER(gvcf,ref)
    annotated=VEP_ANNOTATE(pass,ref,refs.vep_cache)
    PYPGATK_FASTA(annotated,refs.gtf,refs.proteome)
    groups=annotated.branch { m,v,t -> baseline:m.baseline=='true'; progression:m.baseline=='false'; other:true }
    bases=groups.baseline.map { m,v,t -> tuple(m.tk,v,t) }
    pairs=groups.progression.map { m,v,t -> tuple(m.tk,m,v,t) }.combine(bases,by:0).map { k,m,pv,pt,bv,bt -> tuple(m,pv,pt,bv,bt) }
    prog=PROGRESSION_SUBTRACT(pairs)
    PROGRESSION_FASTA(prog,refs.gtf,refs.proteome)
}
