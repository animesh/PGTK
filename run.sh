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
# multiple lanes per sample, so need to rename laneNA to laneNA for Sarek to work
sed  's/NA/laneNA/' PGTK/samples.csv > PGTK/sample_laneNA.csv

nextflow run https://github.com/nf-core/sarek --revision 3.9.0 --input PGTK/sample_laneNA.csv --outdir resultsTKvep -profile singularity -c <(echo "process { withName: 'NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM' { cpus = 12; ext.args = { '-t 12' }; ext.args2 = { '--threads 12' } } }")  --tools haplotypecaller,vep --max_memory '64.GB' --max_cpus 12 -resume
