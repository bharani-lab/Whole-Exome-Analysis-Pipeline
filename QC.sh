#######################
#$1="sample name: example "sample_S1"
#$2="work directory"
######################
#!/usr/bin/env bash
set -o errexit
set -o nounset
mkdir -p output
mkdir -p output/${1}
#fastqc
mkdir -p output/${1}/QC_check
docker run --rm -v "$2":/data pegi3s/fastqc -t 30 /data/data/${1}_R1.fastq.gz -o /data/output/${1}/QC_check
docker run --rm -v "$2":/data pegi3s/fastqc -t 30 /data/data/${1}_R2.fastq.gz -o /data/output/${1}/QC_check
mkdir -p output/$1/trimedreads_fastp
#sudo chmod 777 output/$1/trimedreads_fastp
docker run --rm -v "$2":/data biocontainers/fastp:v0.20.1_cv1 fastp -i data/${1}_R1.fastq.gz -I data/${1}_R2.fastq.gz -o output/$1/trimedreads_fastp/${1}_1.fq.gz -O output/$1/trimedreads_fastp/${1}_2.fq.gz -j output/$1/trimedreads_fastp/${1}.json -h output/$1/trimedreads_fastp/${1}.html -w 16
#Fastqc_fastp
mkdir -p output/${1}/Fastp_fastqc
docker run --rm -v"$2":/data pegi3s/fastqc -t 30 /data/output/$1/trimedreads_fastp/${1}_1.fq.gz -o /data/output/${1}/Fastp_fastqc
docker run --rm -v"$2":/data pegi3s/fastqc -t 30 /data/output/$1/trimedreads_fastp/${1}_2.fq.gz -o /data/output/${1}/Fastp_fastqc
