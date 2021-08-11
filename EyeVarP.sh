#Getting_input_from_file
for r1 in `cat sample.txt`;
do
#sam_File
bwa mem -M -t 20 /Documents/Reference/hg38.fa Data/"$r1"_R1_001.fastq.gz Data/"$r1"_R2_001.fastq.gz > output/"$r1"/"$r1".sam
#Bam_File										  
samtools view -bS output/"$r1"/"$r1".sam -o output/"$r1"/"$r1".bam
#Sort_Bam
java -Xmx100g -jar /Documents/Tools/picard-tools-1.141/picard.jar SortSam VALIDATION_STRINGENCY=SILENT I=output/"$r1"/"$r1".bam O=output/"$r1"/"$r1"_Sort.bam SORT_ORDER=coordinate
#ID_Addition
java -Xmx100g -jar /Documents/Tools/picard-tools-1.141/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT I=output/"$r1"/"$r1"_Sort.bam O=output/"$r1"/"$r1"_RG.bam SO=coordinate RGID=SRR"$r1" RGLB=SRR"$r1" RGPL=illumina RGPU=SRR"$r1" RGSM=SRR"$r1" CREATE_INDEX=true
#Variant calling Gatk
gatk HaplotypeCaller --java-options "-Xmx100g"  -R /Documents/Reference/hg38.fa -I output/"$r1"/"$r1"_RG.bam -O output/"$r1"/"$r1".vcf.gz
#variant Sepration
tabix -h -R /Documents/Reference/Covered_region.bed output/"$r1"/"$r1".vcf.gz > output/"$r1"/Covered_"$r1".vcf
#Annovar
perl /Documents/Tools/annovar/table_annovar.pl output/"$r1"/"$r1".vcf.gz /home/bioinformatics/Documents/COE/Miseq_LCA/COE_RUN3_32LCA/Tools/annovar/humandb/ -buildver hg38 -out output/"$r1"/"$r1"_GATK -otherinfo -remove -protocol refGene,esp6500siv2_all,1000g2015aug_all,cosmic70,esp6500siv2_all,exac03,gnomad30_genome,nci60,clinvar_20210123,avsnp150,ljb26_all,dbnsfp41a,dbscsnv11,intervar_20180118,mcap,revel,regsnpintron -operation g,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput
#Depth
java -jar /Documents/Tools/GenomeAnalysisTK.jar -T DepthOfCoverage -R /Documents/Reference/hg38.fa -I output/"$r1"/"$r1"_RG.bam -L /Documents/Reference/Region.list -omitBaseOutput -omitLocusTable -o output/"$r1"/"$r1"_Depth
#VFC_adding value for nulls 
sed 's/=.;/=-999;/g' output/"$r1"/"$r1"_GATK.hg38_multianno.vcf > output/"$r1"/"$r1"_eyeVarP.vcf
#eyeVarP
printf 'output'/"$r1"/"$r1"'_eyeVarP.vcf\tSRR'"$r1"'' > 'output'/"$r1"/"$r1".txt 
#eyeVarP
python3 eyeVarp/eyeVarP.py priority output/"$r1"/"$r1"_eyeVarP_filter output/"$r1"/"$r1".txt eyeVarP/test_data/default_variants_parameters.txt
#Filtering
Rscript Filtering_variants.r output/"$r1"/"$r1"_GATK.hg38_multianno.txt output/"$r1"/"$r1"_eyeVarP_filter.txt Reference/hc.csv output/"$r1"/"$r1"_eyeVarP_Filter 2>"$r1"_error_filter.txt
Rscript Fitering_gene.r output/"$r1"/"$r1"_eyeVarP_Filter Reference/hc.csv output/"$r1"/"$r1"_eyeVarP_Filter_gene 2>"$r1"_error_gene.txt
#Filtering_Indels and SNVs
grep -Esw 'AAChange.refGene|splicing|frameshift|nonsynonymous|stopgain|stoploss' output/"$r1"/"$r1"_eyeVarP_Filter.csv > output/annovar_output/"$r1"_filtered_variants.csv
#Exomiser
java -Xms20g -Xmx50g -jar exomiser-cli-12.1.0.jar --analysis output/"$r1"/"$r1"_Variant_Sort_analysis-exome.yml
done
