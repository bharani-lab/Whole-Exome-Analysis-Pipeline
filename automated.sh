for r1 in `cat data.txt`;
do
sh QC.sh "$r1" /home/bioinformatics/Music/eyeVarP
sh Alignment.sh "$r1" /home/bioinformatics/Music/eyeVarP /home/bioinformatics/Music/eyeVarP/Refence
sh VariantCalling.sh "$r1" /home/bioinformatics/Music/eyeVarP /home/bioinformatics/Music/eyeVarP/Refence Covered_region.bed
sh VariantAnnotation.sh "$r1" /home/bioinformatics/Music/eyeVarP /home/bioinformatics/Music/eyeVarP/Refence/humandb
sh VarP.sh "$r1" /home/bioinformatics/Music/eyeVarP
sh eyeVarP.sh "$r1" /home/bioinformatics/Music/eyeVarP/output/"$r1"/"$r1"_eyeVarP /home/bioinformatics/Music/eyeVarP/Refence/exomiser_data /home/bioinformatics/Music/eyeVarP /home/bioinformatics/Music/eyeVarP/Refence
done
