#######################
#$1="sample name: example "sample_S1"
#$2="work directory"
#$3="Reference location"
#$4="bed files"
######################
#!/usr/bin/env bash
set -o errexit
set -o nounset
#set -o pipefail
#variant_Calling 
mkdir -p output/"$1"/"$1"_Variant_calling
docker run --rm -v "$2":/data -v "$3":/Reference broadinstitute/gatk gatk HaplotypeCaller --java-options "-Xmx100g" -R /Reference/hg38.fasta -I /data/output/"$1"/"$1"_Add_Read_Group/"$1"_alignment_RG.bam -O /data/output/"$1"/"$1"_Variant_calling/"$1"_GATK.vcf.gz
#variant sepration
docker run --rm -v "$2":/data -v "$3":/Reference bioslimcontainers/tabix:1.7 tabix -h -R Reference/"$4" /data/output/"$1"/"$1"_Variant_calling/"$1"_GATK.vcf.gz > output/"$1"/"$1"_Variant_calling/"$1"_GATK_Covered.vcf
