#######################
#$1="sample name: example "sample_S1"
#$2="work directory"
######################
#!/usr/bin/env bash
set -o errexit
set -o nounset
#VFC_adding value for nulls 
mkdir -p output/"$1"/"$1"_varp
sed 's/=.;/=-999;/g' output/"$1"/"$1"_Variant_annotation/"$1"_annovar.hg38_multianno.vcf > output/"$1"/"$1"_Variant_annotation/"$1"_VarP.vcf
#VarP input
printf 'data/output/'"$1"/"$1"'_Variant_annotation/'"$1"'_VarP.vcf\tSRR'"$1"'' > output/"$1"/"$1"_Variant_annotation/"$1".txt
#VarP
docker run --rm -v "$2":/VarP/data mano2991/varp python3 VarP.py priority data/output/"$1"/"$1"_varp/"$1"_varp data/output/"$1"/"$1"_Variant_annotation/"$1".txt default_0.001_variants_parameters_PPF.txt
#combineVarP with heuristicmethod
docker run --rm -v "$2":/code/data mano2991/eyevarp Rscript Filtering.r data/output/"$1"/"$1"_Variant_annotation/"$1"_annovar.hg38_multianno.txt data/output/"$1"/"$1"_varp/"$1"_varp.txt data/output/"$1"/"$1"_varp/"$1"
grep -Ew 'Gene.refGene|splicing|stopgain|nonsynonymous SNV|frameshift deletion|frameshift insertion|stoploss|startloss' output/"$1"/"$1"_varp/"$1".csv > output/"$1"/"$1"_varp/"$1"_VarP_extracted.csv 
