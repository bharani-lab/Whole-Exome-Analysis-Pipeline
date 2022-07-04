#######################
#$1="sample name: example "sample_S1"
#$2="work directory"
#$3="Reference location"
######################
#!/usr/bin/env bash
set -o errexit
set -o nounset
#set -o pipefail
#BWA_alignment
mkdir -p output/"$1"/"$1"_BWA_alignment
docker run --rm -v "$2":/data -v "$3":/Reference biocontainers/bwa:v0.7.17_cv1 bwa mem -t 30 ../Reference/hg38.fasta  /data/output/"$1"/trimedreads_fastp/"$1"_1.fq.gz /data/output/"$1"/trimedreads_fastp/"$1"_2.fq.gz > output/"$1"/"$1"_BWA_alignment/"$1"_alignment.sam
#SAM to BAM conversion
mkdir -p output/"$1"/"$1"_BWA_alignment_conversion
docker run --rm -v "$2":/data biocontainers/samtools:v1.9-4-deb_cv1 samtools view -@ 30 -bS output/"$1"/"$1"_BWA_alignment/"$1"_alignment.sam -o output/"$1"/"$1"_BWA_alignment_conversion/"$1"_alignment.bam
#SAM SORT
mkdir -p output/"$1"/"$1"_BWA_alignment_sort
docker run --rm -v "$2":/data biocontainers/samtools:v1.9-4-deb_cv1 samtools sort -@ 30 output/"$1"/"$1"_BWA_alignment_conversion/"$1"_alignment.bam -o output/"$1"/"$1"_BWA_alignment_sort/"$1"_alignment_Sort.bam 
#Remove PCR Duplicates
mkdir -p output/"$1"/"$1"_PCR_Duplicates_removal
docker run --rm  -v "$2":/data broadinstitute/picard java -Xmx100g -jar ../../usr/picard/picard.jar MarkDuplicates --VALIDATION_STRINGENCY SILENT -I ../../data/output/"$1"/"$1"_BWA_alignment_sort/"$1"_alignment_Sort.bam -O ../../data/output/"$1"/"$1"_PCR_Duplicates_removal/"$1"_alignment_PCR.bam --REMOVE_DUPLICATES true -M ../../data/output/"$1"/"$1"_PCR_Duplicates_removal/"$1"_alignment_PCR.Metrics
#identification_ADD
mkdir -p output/"$1"/"$1"_Add_Read_Group
docker run --rm -v "$2":/data broadinstitute/picard java -Xmx100g -jar ../../usr/picard/picard.jar AddOrReplaceReadGroups --VALIDATION_STRINGENCY SILENT -I ../../data/output/"$1"/"$1"_PCR_Duplicates_removal/"$1"_alignment_PCR.bam -O ../../data/output/"$1"/"$1"_Add_Read_Group/"$1"_alignment_RG.bam -SO coordinate --RGID SRR"$1" --RGLB SRR"$1" --RGPL illumina --RGPU SRR"$1" --RGSM SRR"$1" --CREATE_INDEX true
