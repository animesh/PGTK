#!/bin/bash
data=$1
#curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA1176350&usehistory=y" | grep -oP '(?<=<QueryKey>)\d+|(?<=<WebEnv>)[^<]+' | { read qk; read we; curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=$qk&WebEnv=$we&rettype=runinfo&retmode=text"; } | tail -n +2 | cut -d',' -f1
#SRR31089075
#SRR31089074
#SRR31089073
#SRR31089072
#SRR31089071
#SRR31089070
#SRR31089076
#curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=PRJNA1176350&usehistory=y" \
#| grep -oP '(?<=<QueryKey>)\d+|(?<=<WebEnv>)[^<]+' \
#| { read qk; read we; curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&query_key=$qk&WebEnv=$we&rettype=runinfo&retmode=text"; } \
#| tail -n +2 | cut -d',' -f1 \
#| xargs -n1 -I{} sh -c 'curl -s "https://www.ebi.ac.uk/ena/portal/api/filereport?accession={}&result=read_run&fields=fastq_ftp" | tail -n +2 | #awk -F"\t" "{print \$NF}" | tr ";," "\n"' \
#| awk '{gsub(/^ftp:\/\//,"https://"); if($0!~/^https?:\/\//) $0="https://"$0; print}' \
#| sort -u | tee data/fastq/ena_fastq_urls.txt \
#| xargs -n1 -P4 wget -c -P data/fastq
#data="/root/Download/rnafusion/data/fastq/"
#for i in data/fastq/*; do echo $i ; zcat $i | wc ; done
#data/fastq/SRR31089070_1.fastq.gz
#280145636 350182045 25946135266
#data/fastq/SRR31089070_2.fastq.gz
#280145636 350182045 25946135266
#data/fastq/SRR31089071_1.fastq.gz
#300948128 376185160 27873514716
#data/fastq/SRR31089071_2.fastq.gz
#300948128 376185160 27873514716
#data/fastq/SRR31089072_1.fastq.gz
#316158656 395198320 29282167108
#data/fastq/SRR31089072_2.fastq.gz
#316158656 395198320 29282167108
#data/fastq/SRR31089073_1.fastq.gz
#348329984 435412480 32262775449
#data/fastq/SRR31089073_2.fastq.gz
#348329984 435412480 32262775449
#data/fastq/SRR31089074_1.fastq.gz
#256067808 320084760 23714971126
#data/fastq/SRR31089074_2.fastq.gz
#256067808 320084760 23714971126
#data/fastq/SRR31089075_1.fastq.gz
#352588332 440735415 32656136952
#data/fastq/SRR31089075_2.fastq.gz
#352588332 440735415 32656136952
#data/fastq/SRR31089076_1.fastq.gz
#2270635568 2838294460 209086300796
#data/fastq/SRR31089076_2.fastq.gz
#2270635568 2838294460 209086300796
echo "sample,fastq_1,fastq_2,strandedness" > $data/samples.csv
rm S?
ls -1 $data/*1.f*q.gz  > S1
ls -1 $data/*2.f*q.gz  > S2
ls -1 $data/*1.f*q.gz | xargs -n 1 basename | awk -F '_' '{print $1}' > S0
printf 'unknown\n%.0s' $data/*1.f*q.gz > S3
paste -d ','  S? >> $data/samples.csv
cat $data/samples.csv
echo "nextflow run nf-core/rnafusion --max_memory '56.GB' --max_cpus 14  -profile docker --input $data/samples.csv --outdir RF  -resume"
#https://nf-co.re/denovotranscript/dev/docs/output/
#bash generateSampleSheet.sh /root/Download/rnafusion/data/fastq/
