# Wole-Exome-Analysis-Pipeline_Manoj

eyeVarP is a modular pipeline

Install Docker for using eyeVarP pipeline 
https://docs.docker.com/engine/install/ubuntu/

Dowload the required Reference file 
https://readycloud.netgear.com/client/browselink.html#t=0pdy234lzn99ykuw8fhsew02f2y/eyeVarP/

also download the required docker conatiners 

docker pull pegi3s/fastqc

docker pull biocontainers/fastp:v0.20.1_cv1

docker pull biocontainers/bwa:v0.7.17_cv1

docker pull biocontainers/samtools:v1.9-4-deb_cv1 

docker pull broadinstitute/picard

docker pull broadinstitute/gatk

docker pull bioslimcontainers/tabix:1.7

docker pull mano2991/annovar:latest

docker pull mano2991/varp

docker pull mano2991/eyevarp


Runing eyeVarP pipline 

sh automated.sh

for moder pipline 

sh QC.sh filename(fastq) fastq location (sh qc.sh SRR123 ~/Desktop/data location)
