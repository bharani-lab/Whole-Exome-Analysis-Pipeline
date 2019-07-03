read -p " enter the Read R1 :" r1
read -p " enter the Read R2 :" r2
#sam_File
bwa mem -M -t 10 Reference/human_g1k_v37.fasta Data/"$r1"_R1.fastq.gz Data/"$r2"_R2.fastq.gz > output/"$r1"/"$r1".sam
#Bam_File
samtools view -@ 10 -bS output/"$r1"/"$r1".sam -o output/"$r1"/"$r1".bam
#Sort_Bam
java -Xmx10g -jar Tools/picard-tools-1.141/picard.jar SortSam VALIDATION_STRINGENCY=SILENT I=output/"$r1"/"$r1".bam O=output/"$r1"/"$r1"_Sort.bam SORT_ORDER=coordinate
#Pcr_Duplicates
java -Xmx10g -jar Tools/picard-tools-1.141/picard.jar MarkDuplicates VALIDATION_STRINGENCY=SILENT I=output/"$r1"/"$r1"_Sort.bam O=output/"$r1"/"$r1"_PCR.bam REMOVE_DUPLICATES=true M=output/"$r1"/"$r1"_pcr.metrics
#ID_Addition
java -Xmx10g -jar Tools/picard-tools-1.141/picard.jar AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT I=output/"$r1"/"$r1"_PCR.bam O=output/"$r1"/"$r1"_RG.bam SO=coordinate RGID=SRR"$r1" RGLB=SRR"$r1" RGPL=illumina RGPU=SRR"$r1" RGSM=SRR"$r1" CREATE_INDEX=true
#Variant_Calling
./Tools/gatk-4.0.6.0/gatk --java-options "-Xmx10g" HaplotypeCaller -R Reference/human_g1k_v37.fasta -I output/"$r1"/"$r1"_RG.bam -O output/"$r1"/BWA_GATK_"$r1".vcf.gz
#variant_Sepration
tabix -h -R Reference/Covered_region.bed output/"$r1"/BWA_GATK_"$r1".vcf.gz > output/"$r1"/Covered_"$r1".vcf
#Annovar
#preparing input for annovar
perl Tools/annovar/convert2annovar.pl -format vcf4 indels_variants.vcf -outfile annovar_indels_variants.avinput
#running annovar
perl Tools/annovar/table_annovar.pl annovar_indels_variants.avinput /home/bioinformatics/Documents/Mano/Tools/annovar/humandb/ -buildver hg19 -out myanno -remove -protocol refGene,cytoBand,genomicSuperDups,esp6500siv2_all,1000g2014oct_all,1000g2014oct_afr,1000g2014oct_eas,1000g2014oct_eur,snp138,ljb26_all -operation g,r,r,f,f,f,f,f,f,f -nastring . -csvout
