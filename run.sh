curl -I https://s3.nird.sigma2.no
java --version
#openjdk 21.0.8 2025-07-15
#git clone https://github.com/animesh/PGTK
#vim ./PGTK/samples.csv
#cd PGTK
#wget -i /cluster/projects/nn9036k/TK/sra/fastq/ena_fastq_urls.txt
#ls -1 TK/*_1.fq.gz | awk -F '_' '{print $1$2$5}' | sed 's/TK\///g' > S1
#ls -1 TK/*_1.fq.gz | awk -F '_' '{print $1}' | sed 's/TK\///g' > S2
#printf 'lane_%s\n' {1..7} > S3
#ls -1 $PWD/TK/*1.fq.gz  > S4
#ls -1 $PWD/TK/*2.fq.gz  > S5
#echo "patient,sample,lane,fastq_1,fastq_2" > samples.csv
#paste -d ',' S? >> samples.csv
#nextflow  self-update
/cluster/home/ash022/scripts/nextflow -v
#nextflow version 26.04.6.12646
sudo apt purge apptainer apptainer-suid
sudo apt autoremove
# Add the key
curl -s https://packages.microsoft.com/keys/microsoft.asc | gpg --dearmor | sudo tee /usr/share/keyrings/apptainer.gpg > /dev/null
# Add the repository to your sources list
echo "deb [signed-by=/usr/share/keyrings/apptainer.gpg] https://packages.microsoft.com/debian/$(lsb_release -cs)/prod $(lsb_release -cs) main" | sudo tee /etc/apt/sources.list.d/apptainer.list
# Update and install
sudo apt update
sudo apt install -y apptainer
# lane info
awk -F',' 'NR>1 {f=$3; print "\n=== " $2 " (" f ") ==="; system("zcat " f " | awk \"NR%4==1 {split(\$2, p, \\\":\\\"); if(p[4] ~ /^[0-9]+$/) valid[p[4]]++; else invalid++} END {for(l in valid) print \\\"  Lane \\\" l \\\": \\\" valid[l] \\\" reads\\\"; if(invalid) print \\\"  INVALID/MISSING: \\\" invalid \\\" reads\\\"}\"")}' PGTK/samples.csv
=== SRR31089070 (/mnt/z/Download/TK/fastq/SRR31089070_1.fastq.gz) ===
  Lane 4: 70036409 reads

=== SRR31089071 (/mnt/z/Download/TK/fastq/SRR31089071_1.fastq.gz) ===
  Lane 1: 22312937 reads
  Lane 4: 52924095 reads

=== SRR31089072 (/mnt/z/Download/TK/fastq/SRR31089072_1.fastq.gz) ===
  Lane 7: 79039664 reads

=== SRR31089073 (/mnt/z/Download/TK/fastq/SRR31089073_1.fastq.gz) ===
  Lane 3: 27617267 reads
  Lane 2: 28335788 reads
  Lane 7: 31129441 reads

=== SRR31089074 (/mnt/z/Download/TK/fastq/SRR31089074_1.fastq.gz) ===
  Lane 1: 43538252 reads
  Lane 4: 20478700 reads

=== SRR31089075 (/mnt/z/Download/TK/fastq/SRR31089075_1.fastq.gz) ===
  Lane 4: 88147083 reads

=== SRR31089076 (/mnt/z/Download/TK/fastq/SRR31089076_1.fastq.gz) ===
  Lane 3: 189222463 reads
  Lane 6: 86114371 reads
  Lane 5: 99913894 reads
  Lane 2: 192408164 reads
# multiple lanes per sample, so need to rename laneNA to lane0 for Sarek to work
sed  's/NA/lane0/' PGTK/samples_laneNA.csv > PGTK/sample_lane0.csv

nextflow clean -f -q

nextflow run https://github.com/nf-core/sarek --revision 3.9.0 --input PGTK/sample_lane0.csv --outdir resultsTKvep -profile singularity -c <(echo "process { withName: 'NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM' { cpus = 12; ext.args = { '-t 12' }; ext.args2 = { '--threads 12' } } }")  --tools haplotypecaller,vep --max_memory '64.GB' --max_cpus 12 -w /mnt/z/Download/

nextflow clean -f -q

nextflow run https://github.com/nf-core/sarek \
    --revision 3.9.0 \
    --input PGTK/sample_lane0.csv \
    --outdir resultsTKvep \
    -profile singularity \
    -c <(echo "process { withName: 'NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM' { cpus = 12 } }") \
    --tools haplotypecaller,vep \
    --max_memory '64.GB' \
    --max_cpus 12 \
    -w /mnt/z/Download/


 N E X T F L O W   ~  version 26.04.6

Launching `https://github.com/nf-core/sarek` [backstabbing_baekeland] revision: b97952e5ba [master]


------------------------------------------------------
                                        ,--./,-.
        ___     __   __   __   ___     /,-._.--~'
  |\ | |__  __ /  ` /  \ |__) |__         }  {
  | \| |       \__, \__/ |  \ |___     \`-._,-`-,
                                        `._,._,'
      ____
    .´ _  `.
   /  |\`-_ \      __        __   ___     
  |   | \  `-|    |__`  /\  |__) |__  |__/
   \ |   \  /     .__| /¯¯\ |  \ |___ |  \
    `|____\´

  nf-core/sarek 3.9.0
------------------------------------------------------

Input/output options
  input                  : PGTK/sample_lane0.csv
  outdir                 : resultsTKvep

Main options
  intervals              : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed
  tools                  : haplotypecaller,vep

Variant Calling
  cf_chrom_len           : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/Length/Homo_sapiens_assembly38.len
  pon                    : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz
  pon_tbi                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz.tbi

Reference genome options
  ascat_genome           : hg38
  ascat_alleles          : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_alleles_hg38.zip
  ascat_loci             : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/G1000_loci_hg38.zip
  ascat_loci_gc          : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/GC_G1000_hg38.zip
  ascat_loci_rt          : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/ASCAT/RT_G1000_hg38.zip
  bwa                    : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/BWAIndex/
  bwamem2                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/BWAmem2Index/
  chr_dir                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/Chromosomes
  dbsnp                  : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz
  dbsnp_tbi              : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/dbsnp_146.hg38.vcf.gz.tbi
  dbsnp_vqsr             : --resource:dbsnp,known=false,training=true,truth=false,prior=2.0 dbsnp_146.hg38.vcf.gz
  dict                   : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.dict
  dragmap                : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/dragmap/
  fasta                  : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
  fasta_fai              : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta.fai
  germline_resource      : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz
  germline_resource_tbi  : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/af-only-gnomad.hg38.vcf.gz.tbi
  known_indels           : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz
  known_indels_tbi       : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/{Mills_and_1000G_gold_standard.indels.hg38,beta/Homo_sapiens_assembly38.known_indels}.vcf.gz.tbi
  known_indels_vqsr      : --resource:gatk,known=false,training=true,truth=true,prior=10.0 Homo_sapiens_assembly38.known_indels.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=10.0 Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
  known_snps             : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz
  known_snps_tbi         : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000G_omni2.5.hg38.vcf.gz.tbi
  known_snps_vqsr        : --resource:1000G,known=false,training=true,truth=true,prior=10.0 1000G_omni2.5.hg38.vcf.gz
  mappability            : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/Control-FREEC/out100m2_hg38.gem
  msisensor2_models      : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/MSIsensor2/models_hg38//
  msisensorpro_scan      : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/MSIsensorPro/Homo_sapiens_assembly38.msisensor_scan.list
  ngscheckmate_bed       : s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/NGSCheckMate/SNP_GRCh38_hg38_wChr.bed
  sentieon_dnascope_model: s3://ngi-igenomes/igenomes//Homo_sapiens/GATK/GRCh38/Annotation/Sentieon/SentieonDNAscopeModel1.1.model
  snpeff_db              : GRCh38.99
  vep_cache_version      : 115
  vep_genome             : GRCh38
  vep_species            : homo_sapiens

Generic options
  trace_report_suffix    : 2026-07-25_21-41-16

Core Nextflow options
  revision               : master
  runName                : backstabbing_baekeland
  containerEngine        : singularity
  launchDir              : /root
  workDir                : /mnt/z/Download
  projectDir             : /root/.nextflow/assets/.repos/nf-core/sarek/clones/b97952e5bac68d5deb93d4a3349a45f146be9830
  userName               : root
  profile                : singularity
  configFiles            : /root/.nextflow/assets/.repos/nf-core/sarek/clones/b97952e5bac68d5deb93d4a3349a45f146be9830/nextflow.config, /dev/fd/63

!! Only displaying parameters that differ from the pipeline defaults !!
------------------------------------------------------

* The pipeline
    https://doi.org/10.12688/f1000research.16665.2
    https://doi.org/10.1093/nargab/lqae031
    https://doi.org/10.5281/zenodo.3476425

* The nf-core framework
    https://doi.org/10.1038/s41587-020-0439-x

* Software dependencies
    https://github.com/nf-core/sarek/blob/master/CITATIONS.md

executor >  local (1)
[-        ] NFC…ALS:CREATE_INTERVALS_BED -
executor >  local (1)
[-        ] NFC…ALS:CREATE_INTERVALS_BED -
[-        ] NFC…GZIPTABIX_INTERVAL_SPLIT -
[-        ] NFC…PTABIX_INTERVAL_COMBINED -
executor >  local (1)
[-        ] NFC…ALS:CREATE_INTERVALS_BED | 0 of 1
[-        ] NFC…GZIPTABIX_INTERVAL_SPLIT -
[-        ] NFC…PTABIX_INTERVAL_COMBINED | 0 of 1
[-        ] NFC…NG_DECOMPRESS_TO_FQ_PAIR -
[-        ] NFC…RING_DECOMPRESS_TO_R1_FQ -
[-        ] NFC…RING_DECOMPRESS_TO_R2_FQ -
[-        ] NFC…UT:SAMTOOLS_VIEW_MAP_MAP -
[-        ] NFCORE_SAREK:SAREK:FASTQC    | 0 of 3
[8e/344995] NFC…ASTP (SRR31089072-lane0) | 0 of 3
Plus 37 more processes waiting for tasks…
Staging foreign file: s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/Annotation/intervals/wgs_calling_regions_noseconds.hg38.bed ...
...

[-        ] NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                  -
[ac/2f866d] NFCORE_SAREK:SAREK:FASTQC (SRR31089073-lane0)                                 [100%] 3 of 3 ✔
[bd/848566] NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTP (SRR31089072-lane0)            [100%] 3 of 3 ✔
[4c/e4ff70] NFC…E_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM (SRR31089072) [100%] 36 of 36 ✔
[50/59b898] NFC…STQ_PREPROCESS_GATK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (SRR31089072) [100%] 3 of 3 ✔
[56/7f72d2] NFC…BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS (SRR31089072) [100%] 3 of 3 ✔
[ea/c57e54] NFC…_GATK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH (SRR31089072) [100%] 3 of 3 ✔
[ba/6fc41e] NFC…PREPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (SRR31089072) [100%] 63 of 63 ✔
[f3/d0d53d] NFC…REPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_GATHERBQSRREPORTS (SRR31089072) [100%] 3 of 3 ✔
[70/6b35d8] NFC…K:SAREK:FASTQ_PREPROCESS_GATK:BAM_APPLYBQSR:GATK4_APPLYBQSR (SRR31089072) [100%] 63 of 63 ✔
[b4/ebcf44] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:MERGE_CRAM (SRR31089072) [100%] 3 of 3 ✔
[4e/a87cb8] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:INDEX_CRAM (SRR31089072) [100%] 3 of 3 ✔
[34/0ba974] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:SAMTOOLS_STATS (SRR31089073)   [100%] 3 of 3 ✔
[05/d71257] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:MOSDEPTH (SRR31089072)         [100%] 3 of 3 ✔
[93/45af9a] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:GATK4_HAPLOTYPECALLER (SRR31089072) [100%] 63 of 63 ✔
[ad/271104] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:MERGE_HAPLOTYPECALLER (SRR31089072) [100%] 3 of 3 ✔
[48/82e444] NFC…NG_GERMLINE_ALL:VCF_VARIANT_FILTERING_GATK:CNNSCOREVARIANTS (SRR31089072) [100%] 3 of 3 ✔
[9d/c31259] NFC…RMLINE_ALL:VCF_VARIANT_FILTERING_GATK:FILTERVARIANTTRANCHES (SRR31089072) [100%] 3 of 3 ✔
[60/383a35] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:BCFTOOLS_STATS (SRR31089072)      [100%] 3 of 3 ✔
[bf/b96405] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_COUNT (SRR31089072) [100%] 3 of 3 ✔
[4a/b6231e] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_QUAL (SRR31089072)  [100%] 3 of 3 ✔
executor >  local (303)
[6e/f45ce4] NFC…E_INTERVALS:CREATE_INTERVALS_BED (wgs_calling_regions_noseconds.hg38.bed) [100%] 1 of 1 ✔
[6f/8df7a7] NFC…PARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_SPLIT (chr13_86252980-111703855) [100%] 21 of 21 ✔
[41/bfd37a] NFC…S:TABIX_BGZIPTABIX_INTERVAL_COMBINED (wgs_calling_regions_noseconds.hg38) [100%] 1 of 1 ✔
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_FQ_PAIR                               -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R1_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R2_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                  -
[ac/2f866d] NFCORE_SAREK:SAREK:FASTQC (SRR31089073-lane0)                                 [100%] 3 of 3 ✔
[bd/848566] NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTP (SRR31089072-lane0)            [100%] 3 of 3 ✔
[4c/e4ff70] NFC…E_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM (SRR31089072) [100%] 36 of 36 ✔
[50/59b898] NFC…STQ_PREPROCESS_GATK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (SRR31089072) [100%] 3 of 3 ✔
[56/7f72d2] NFC…BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS (SRR31089072) [100%] 3 of 3 ✔
[ea/c57e54] NFC…_GATK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH (SRR31089072) [100%] 3 of 3 ✔
[ba/6fc41e] NFC…PREPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (SRR31089072) [100%] 63 of 63 ✔
[f3/d0d53d] NFC…REPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_GATHERBQSRREPORTS (SRR31089072) [100%] 3 of 3 ✔
[70/6b35d8] NFC…K:SAREK:FASTQ_PREPROCESS_GATK:BAM_APPLYBQSR:GATK4_APPLYBQSR (SRR31089072) [100%] 63 of 63 ✔
[b4/ebcf44] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:MERGE_CRAM (SRR31089072) [100%] 3 of 3 ✔
[4e/a87cb8] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:INDEX_CRAM (SRR31089072) [100%] 3 of 3 ✔
[34/0ba974] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:SAMTOOLS_STATS (SRR31089073)   [100%] 3 of 3 ✔
[05/d71257] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:MOSDEPTH (SRR31089072)         [100%] 3 of 3 ✔
[93/45af9a] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:GATK4_HAPLOTYPECALLER (SRR31089072) [100%] 63 of 63 ✔
[ad/271104] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:MERGE_HAPLOTYPECALLER (SRR31089072) [100%] 3 of 3 ✔
[48/82e444] NFC…NG_GERMLINE_ALL:VCF_VARIANT_FILTERING_GATK:CNNSCOREVARIANTS (SRR31089072) [100%] 3 of 3 ✔
[9d/c31259] NFC…RMLINE_ALL:VCF_VARIANT_FILTERING_GATK:FILTERVARIANTTRANCHES (SRR31089072) [100%] 3 of 3 ✔
[60/383a35] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:BCFTOOLS_STATS (SRR31089072)      [100%] 3 of 3 ✔
[bf/b96405] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_COUNT (SRR31089072) [100%] 3 of 3 ✔
[4a/b6231e] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_QUAL (SRR31089072)  [100%] 3 of 3 ✔
executor >  local (303)
[6e/f45ce4] NFC…E_INTERVALS:CREATE_INTERVALS_BED (wgs_calling_regions_noseconds.hg38.bed) [100%] 1 of 1 ✔
[6f/8df7a7] NFC…PARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_SPLIT (chr13_86252980-111703855) [100%] 21 of 21 ✔
[41/bfd37a] NFC…S:TABIX_BGZIPTABIX_INTERVAL_COMBINED (wgs_calling_regions_noseconds.hg38) [100%] 1 of 1 ✔
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_FQ_PAIR                               -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R1_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R2_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                  -
[ac/2f866d] NFCORE_SAREK:SAREK:FASTQC (SRR31089073-lane0)                                 [100%] 3 of 3 ✔
[bd/848566] NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTP (SRR31089072-lane0)            [100%] 3 of 3 ✔
[4c/e4ff70] NFC…E_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM (SRR31089072) [100%] 36 of 36 ✔
[50/59b898] NFC…STQ_PREPROCESS_GATK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (SRR31089072) [100%] 3 of 3 ✔
[56/7f72d2] NFC…BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS (SRR31089072) [100%] 3 of 3 ✔
[ea/c57e54] NFC…_GATK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH (SRR31089072) [100%] 3 of 3 ✔
[ba/6fc41e] NFC…PREPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (SRR31089072) [100%] 63 of 63 ✔
[f3/d0d53d] NFC…REPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_GATHERBQSRREPORTS (SRR31089072) [100%] 3 of 3 ✔
[70/6b35d8] NFC…K:SAREK:FASTQ_PREPROCESS_GATK:BAM_APPLYBQSR:GATK4_APPLYBQSR (SRR31089072) [100%] 63 of 63 ✔
[b4/ebcf44] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:MERGE_CRAM (SRR31089072) [100%] 3 of 3 ✔
[4e/a87cb8] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:INDEX_CRAM (SRR31089072) [100%] 3 of 3 ✔
[34/0ba974] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:SAMTOOLS_STATS (SRR31089073)   [100%] 3 of 3 ✔
[05/d71257] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:MOSDEPTH (SRR31089072)         [100%] 3 of 3 ✔
[93/45af9a] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:GATK4_HAPLOTYPECALLER (SRR31089072) [100%] 63 of 63 ✔
[ad/271104] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:MERGE_HAPLOTYPECALLER (SRR31089072) [100%] 3 of 3 ✔
[48/82e444] NFC…NG_GERMLINE_ALL:VCF_VARIANT_FILTERING_GATK:CNNSCOREVARIANTS (SRR31089072) [100%] 3 of 3 ✔
[9d/c31259] NFC…RMLINE_ALL:VCF_VARIANT_FILTERING_GATK:FILTERVARIANTTRANCHES (SRR31089072) [100%] 3 of 3 ✔
[60/383a35] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:BCFTOOLS_STATS (SRR31089072)      [100%] 3 of 3 ✔
[bf/b96405] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_COUNT (SRR31089072) [100%] 3 of 3 ✔
[4a/b6231e] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_QUAL (SRR31089072)  [100%] 3 of 3 ✔
executor >  local (303)
[6e/f45ce4] NFC…E_INTERVALS:CREATE_INTERVALS_BED (wgs_calling_regions_noseconds.hg38.bed) [100%] 1 of 1 ✔
[6f/8df7a7] NFC…PARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_SPLIT (chr13_86252980-111703855) [100%] 21 of 21 ✔
[41/bfd37a] NFC…S:TABIX_BGZIPTABIX_INTERVAL_COMBINED (wgs_calling_regions_noseconds.hg38) [100%] 1 of 1 ✔
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_FQ_PAIR                               -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R1_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R2_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                  -
[ac/2f866d] NFCORE_SAREK:SAREK:FASTQC (SRR31089073-lane0)                                 [100%] 3 of 3 ✔
[bd/848566] NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTP (SRR31089072-lane0)            [100%] 3 of 3 ✔
[4c/e4ff70] NFC…E_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM (SRR31089072) [100%] 36 of 36 ✔
[50/59b898] NFC…STQ_PREPROCESS_GATK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (SRR31089072) [100%] 3 of 3 ✔
[56/7f72d2] NFC…BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS (SRR31089072) [100%] 3 of 3 ✔
[ea/c57e54] NFC…_GATK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH (SRR31089072) [100%] 3 of 3 ✔
[ba/6fc41e] NFC…PREPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (SRR31089072) [100%] 63 of 63 ✔
[f3/d0d53d] NFC…REPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_GATHERBQSRREPORTS (SRR31089072) [100%] 3 of 3 ✔
[70/6b35d8] NFC…K:SAREK:FASTQ_PREPROCESS_GATK:BAM_APPLYBQSR:GATK4_APPLYBQSR (SRR31089072) [100%] 63 of 63 ✔
[b4/ebcf44] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:MERGE_CRAM (SRR31089072) [100%] 3 of 3 ✔
[4e/a87cb8] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:INDEX_CRAM (SRR31089072) [100%] 3 of 3 ✔
[34/0ba974] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:SAMTOOLS_STATS (SRR31089073)   [100%] 3 of 3 ✔
[05/d71257] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:MOSDEPTH (SRR31089072)         [100%] 3 of 3 ✔
[93/45af9a] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:GATK4_HAPLOTYPECALLER (SRR31089072) [100%] 63 of 63 ✔
[ad/271104] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:MERGE_HAPLOTYPECALLER (SRR31089072) [100%] 3 of 3 ✔
[48/82e444] NFC…NG_GERMLINE_ALL:VCF_VARIANT_FILTERING_GATK:CNNSCOREVARIANTS (SRR31089072) [100%] 3 of 3 ✔
[9d/c31259] NFC…RMLINE_ALL:VCF_VARIANT_FILTERING_GATK:FILTERVARIANTTRANCHES (SRR31089072) [100%] 3 of 3 ✔
executor >  local (303)
[6e/f45ce4] NFC…E_INTERVALS:CREATE_INTERVALS_BED (wgs_calling_regions_noseconds.hg38.bed) [100%] 1 of 1 ✔
[6f/8df7a7] NFC…PARE_INTERVALS:TABIX_BGZIPTABIX_INTERVAL_SPLIT (chr13_86252980-111703855) [100%] 21 of 21 ✔
[41/bfd37a] NFC…S:TABIX_BGZIPTABIX_INTERVAL_COMBINED (wgs_calling_regions_noseconds.hg38) [100%] 1 of 1 ✔
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_FQ_PAIR                               -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R1_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:SPRING_DECOMPRESS_TO_R2_FQ                                 -
[-        ] NFCORE_SAREK:SAREK:CONVERT_FASTQ_INPUT:SAMTOOLS_VIEW_MAP_MAP                  -
[ac/2f866d] NFCORE_SAREK:SAREK:FASTQC (SRR31089073-lane0)                                 [100%] 3 of 3 ✔
[bd/848566] NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTP (SRR31089072-lane0)            [100%] 3 of 3 ✔
[4c/e4ff70] NFC…E_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM (SRR31089072) [100%] 36 of 36 ✔
[50/59b898] NFC…STQ_PREPROCESS_GATK:BAM_MARKDUPLICATES:GATK4_MARKDUPLICATES (SRR31089072) [100%] 3 of 3 ✔
[56/7f72d2] NFC…BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:SAMTOOLS_STATS (SRR31089072) [100%] 3 of 3 ✔
[ea/c57e54] NFC…_GATK:BAM_MARKDUPLICATES:CRAM_QC_MOSDEPTH_SAMTOOLS:MOSDEPTH (SRR31089072) [100%] 3 of 3 ✔
[ba/6fc41e] NFC…PREPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_BASERECALIBRATOR (SRR31089072) [100%] 63 of 63 ✔
[f3/d0d53d] NFC…REPROCESS_GATK:BAM_BASERECALIBRATOR:GATK4_GATHERBQSRREPORTS (SRR31089072) [100%] 3 of 3 ✔
[70/6b35d8] NFC…K:SAREK:FASTQ_PREPROCESS_GATK:BAM_APPLYBQSR:GATK4_APPLYBQSR (SRR31089072) [100%] 63 of 63 ✔
[b4/ebcf44] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:MERGE_CRAM (SRR31089072) [100%] 3 of 3 ✔
[4e/a87cb8] NFC…ESS_GATK:BAM_APPLYBQSR:CRAM_MERGE_INDEX_SAMTOOLS:INDEX_CRAM (SRR31089072) [100%] 3 of 3 ✔
[34/0ba974] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:SAMTOOLS_STATS (SRR31089073)   [100%] 3 of 3 ✔
[05/d71257] NFCORE_SAREK:SAREK:CRAM_SAMPLEQC:CRAM_QC_RECAL:MOSDEPTH (SRR31089072)         [100%] 3 of 3 ✔
[93/45af9a] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:GATK4_HAPLOTYPECALLER (SRR31089072) [100%] 63 of 63 ✔
[ad/271104] NFC…L:BAM_VARIANT_CALLING_HAPLOTYPECALLER:MERGE_HAPLOTYPECALLER (SRR31089072) [100%] 3 of 3 ✔
[48/82e444] NFC…NG_GERMLINE_ALL:VCF_VARIANT_FILTERING_GATK:CNNSCOREVARIANTS (SRR31089072) [100%] 3 of 3 ✔
[9d/c31259] NFC…RMLINE_ALL:VCF_VARIANT_FILTERING_GATK:FILTERVARIANTTRANCHES (SRR31089072) [100%] 3 of 3 ✔
[60/383a35] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:BCFTOOLS_STATS (SRR31089072)      [100%] 3 of 3 ✔
[bf/b96405] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_COUNT (SRR31089072) [100%] 3 of 3 ✔
[4a/b6231e] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_TSTV_QUAL (SRR31089072)  [100%] 3 of 3 ✔
[6f/f7d323] NFCORE_SAREK:SAREK:VCF_QC_BCFTOOLS_VCFTOOLS:VCFTOOLS_SUMMARY (SRR31089072)    [100%] 3 of 3 ✔
[16/63902b] NFCORE_SAREK:SAREK:VCF_ANNOTATE_ALL:ENSEMBLVEP_VEP (SRR31089072)              [100%] 3 of 3 ✔
[0e/4fbb2c] NFCORE_SAREK:SAREK:MULTIQC (sarek)                                            [100%] 1 of 1 ✔
Plus 16 more processes waiting for tasks…
-[nf-core/sarek] Pipeline completed successfully-

Outputs:

  /root/resultsTKvep

  multiqc: multiqc/index.json

Completed at: 27-Jul-2026 11:27:55
Duration    : 1d 1h 55m 18s
CPU hours   : 272.4
Succeeded   : 303

