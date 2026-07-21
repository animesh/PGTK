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
nextflow run https://github.com/nf-core/sarek --revision 3.9.0 --input ./PGTK/samples.csv --outdir resultsTKvep -profile singularity -c <(echo "process { withName: 'NFCORE_SAREK:SAREK:FASTQ_PREPROCESS_GATK:FASTQ_ALIGN:BWAMEM1_MEM' { cpus = 12; ext.args = { '-t 12' }; ext.args2 = { '--threads 12' } } }")  --tools haplotypecaller,vep --max_memory '64.GB' --max_cpus 12 -resume
