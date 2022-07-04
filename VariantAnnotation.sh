#######################
#$1="sample name: example "sample_S1"
#$2="work directory"
#$3="Reference location"
######################
#!/usr/bin/env bash
set -o errexit
set -o nounset
##set -o pipefail
##Variant_annotation
mkdir -p output/"$1"/"$1"_Variant_annotation
##annovar annotation
docker run --rm -it -v "$2":/data -v "$3":/Reference mano2991/annovar:latest perl annovar/table_annovar.pl data/output/"$1"/"$1"_Variant_calling/"$1"_GATK_Covered.vcf Reference/ -buildver hg38 -out data/output/"$1"/"$1"_Variant_annotation/"$1"_annovar -otherinfo -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,cosmic70,esp6500siv2_all,exac03,gnomad30_genome,nci60,clinvar_20210123,avsnp150,ljb26_all,dbnsfp41a,dbscsnv11,intervar_20180118,mcap,revel,regsnpintron -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
