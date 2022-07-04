#######################
#$1="sample name: example "sample_S1"
#$2="work directory"
#$3="Reference for exomiser"
#$4="eyeVarP output location"
#$5="Reference for Gene Ranking"
######################
#!/usr/bin/env bash
set -o errexit
set -o nounset
#VFC_adding value for nulls 
mkdir -p output/"$1"/"$1"_eyeVarP
cp Refence/test-analysis-exome.yml output/"$1"/"$1"_eyeVarP
cp output/"$1"/"$1"_Variant_calling/"$1"_GATK_Covered.vcf output/"$1"/"$1"_eyeVarP/
#exomiser
docker run --rm -v "$2":/"$1"_eyeVarP -v "$3":/exomiser-cli-12.1.0/data mano2991/exomiser java -jar exomiser-cli-12.1.0.jar --analysis ../"$1"_eyeVarP/test-analysis-exome.yml
#eyeVarP Model
docker run --rm -v "$4":/code/data -v "$5":/code/Reference mano2991/eyevarp Rscript eyeVarP.r data/output/"$1"/"$1"_eyeVarP/"$1"_GATK_Covered.variants.tsv data/output/"$1"/"$1"_varp/"$1"_VarP_extracted.csv Reference/Gene_ranking.csv data/output/TS49/TS49_eyeVarP/TS49_eyeVarP_final

