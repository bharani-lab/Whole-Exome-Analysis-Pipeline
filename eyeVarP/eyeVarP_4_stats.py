###########################################################################################################
# #
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
import eyeVarP_conf #
from shutil import copyfile #
#
#
###########################################################################################################
# Define global variables
##########################################################################################################
suffix=str(int(time.time())) # get a unique timestamp for suffix 
#
supplied_args=0 #
#
Sample_loc=-1 #
#
info_msg1_1="eyeVarP stats : Invalid number of inputs, must have two :" 
info_msg1_2="eyeVarP stats : 1) output destination directory + prefix" #
info_msg1_3="eyeVarP stats : 2) Input file - output from eyeVarP prioritisation process" 
info_msg1_4="eyeVarP stats : 3) percentage value - output gene breakdown for variant above this percentage (default is 75% ranking value)" 
#
Ast_ln=eyeVarP_conf.nl+"*********************************************************************************************"+eyeVarP_conf.nl #
############################################################################################################
###########################################################################################################
#
###########################################################################################################
def initial_setup():
#
	global supplied_args 
	global full_variant_file #
	global percentile_variant_file #
	global sample_percentile_variant_file #
	#
	print ("initial_setup():") #
	print (suffix) #
#	print (sys.argv) #
	supplied_args=len(sys.argv) #
#	print (supplied_args) #
#
	if (supplied_args > 5 ):  # arg [0] is the python program
		print (info_msg1_1+eyeVarP_conf.nl+info_msg1_2+eyeVarP_conf.nl+info_msg1_3+eyeVarP_conf.nl+info_msg1_4) #
		return 1 #
	else :
		eyeVarP_conf.output_dir=sys.argv[2] #
		eyeVarP_conf.input_file=sys.argv[3] #
		eyeVarP_conf.parameter_file=75 # set default
#		print (supplied_args)
		if (supplied_args > 4 ):  # there is a percentile supplied
			if (eyeVarP_conf.is_number(sys.argv[4])): # numeric
				eyeVarP_conf.parameter_file=sys.argv[4] #
				if ((int(eyeVarP_conf.parameter_file) < 1) or (int(eyeVarP_conf.parameter_file) > 99)): # not a valid percentile
					eyeVarP_conf.parameter_file=75 # set default to all
#		print (eyeVarP_conf.parameter_file)
	#
		eyeVarP_conf.final_output_file=eyeVarP_conf.output_dir+"variant_statistic_file_"+suffix+".txt" #
		eyeVarP_conf.working_file1=eyeVarP_conf.output_dir+"working_file1_"+suffix+"_tmp.txt" #
		eyeVarP_conf.working_file2=eyeVarP_conf.output_dir+"working_file2_"+suffix+"_tmp.txt" #
		eyeVarP_conf.full_file1=eyeVarP_conf.output_dir+"full_file1_"+suffix+"_tmp.txt" #
		eyeVarP_conf.full_file2=eyeVarP_conf.output_dir+"full_file2_"+suffix+"_tmp.txt" #
		eyeVarP_conf.sort_file1=eyeVarP_conf.output_dir+"sort_file1_"+suffix+"_tmp.txt" #
		eyeVarP_conf.sort_file2=eyeVarP_conf.output_dir+"sort_file2_"+suffix+"_tmp.txt" #
		full_variant_file=eyeVarP_conf.output_dir+"full_variant_file_"+suffix+"_tmp.txt" #
		percentile_variant_file=eyeVarP_conf.output_dir+"percentile_variant_file_"+suffix+"_tmp.txt" #
		sample_percentile_variant_file=eyeVarP_conf.output_dir+"Spercentile_variant_file_"+suffix+"_tmp.txt" #
		print ("output : ",eyeVarP_conf.final_output_file) #
	
	return 0 #
#
###########################################################################################################
#
###########################################################################################################
	
def file_lines_input(working_fn): #
#
#	global scorefile #
	i1=[0,0] #
#	print "Place input into array : " #
##
	lr=0 #
	for line in open(working_fn): 
		lr = lr+1 #
		if ( lr == 1 ):  # first line in file
			this_line=re.split('\t|\n|\r|=|,',line.rstrip()) #
			lc = len(this_line) 

	i1[0]=lr #
	i1[1]=lc
#	print (lc, i1[1])#
#
	COMMAND="awk '($2 !=\"0\") {print $0}' "+working_fn+" > "+full_variant_file # create a scoring variants file 
	subprocess.call(COMMAND, shell=True) # do it in shell
#	print("score file finished")
#	print (eyeVarP_conf.Sample_ids) #
#	print (eyeVarP_conf.txt_start, eyeVarP_conf.txt_end)
	return i1 #
###########################################################################################################
#
###########################################################################################################
	
def unique_genes(working_fn): #
#
	COMMAND="tail -n+2 "+working_fn+" | cut -f3 > "+eyeVarP_conf.sort_file1 # take only the gene names 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="sort -u "+eyeVarP_conf.sort_file1+" > "+eyeVarP_conf.sort_file2 # sort to find unique gene names 
	subprocess.call(COMMAND, shell=True) # do it in shell
	COMMAND="grep -v \",\" "+eyeVarP_conf.sort_file2+" > "+eyeVarP_conf.sort_file1 # remove all the multiple gene names entries - leave unique genes 
	subprocess.call(COMMAND, shell=True) # do it in shell
	return 0 #
###########################################################################################################
#
###########################################################################################################
def Total_variants_stats(Output_file): #
#
	global gene_loc #
	global infile_shape #
#
	Outln=Ast_ln+eyeVarP_conf.nl+"TOTAL STATS FOR - "+eyeVarP_conf.input_file+eyeVarP_conf.nl #
#
	j=[0,0] #
#	print "Place input into array : " #
##
	lr=0 #
	for line in open(full_variant_file):  # count how many scoring variants
		lr = lr+1 #
		if ( lr == 1 ):  # first line in file - header
			header_line=re.split('\t|\n|\r|=|,',line.rstrip()) #
			lc = len(header_line) 
		if ( lr == 2 ):  # first variant line in file
			variant_line=re.split('\t|\n|\r|=|,',line.rstrip()) #
			top_score = variant_line[1] 
			print ("top score :",top_score)

	j[0]=lr #
	j[1]=lc
#	print (lc, i1[1])#
#
	num_var=infile_shape[0]-1 # number of variants in input file - assume a header
	eyeVarP_conf.txt_start = 3 # samples id starts here
	eyeVarP_conf.txt_end = infile_shape[1] # -1 for the blank element 
	eyeVarP_conf.Sample_ids = header_line[eyeVarP_conf.txt_start:eyeVarP_conf.txt_end] #
#	print (eyeVarP_conf.Sample_ids)
	num_samples = k = eyeVarP_conf.txt_end-eyeVarP_conf.txt_start #
	Outln=Outln+"Total number of samples : "+str(num_samples)+eyeVarP_conf.nl #
#
# find unique genes that have scoring variants
	unique_genes(full_variant_file) #
	lr=0 #
	for line in open(eyeVarP_conf.sort_file1): 
		lr = lr+1 #
	Outln=Outln+"Total number of genes : "+str(lr)+eyeVarP_conf.nl #
	Outln=Outln+"Total number of variants : "+str(num_var)+eyeVarP_conf.nl #
	Outln=Outln+"Total number of non-scoring variants : "+str(infile_shape[0]-j[0])+eyeVarP_conf.nl #
	Outln=Outln+"Total number of scoring variants : "+str(j[0]-1)+eyeVarP_conf.nl+eyeVarP_conf.nl #
#
# determine breakdown of percentage of ranking with the scorefile (variants that have scores)
#	print("breakdown start "+eyeVarP_conf.nl)
	percentage=np.array([0.99,0.95,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0.01]) #
	for i in range(len(percentage)): #
		working_score=int(int(top_score)*percentage[i])
#		print (str(working_score))
		pc_val=percentage[i]*100 #
#		COMMAND="awk '{ if ($1 >= "+str(percentage[i])+") print $0}' "+full_variant_file+" > "+eyeVarP_conf.sort_file1 # create a working input file 
		COMMAND="awk '{ if ($2 >= "+str(working_score)+") print $0}' "+full_variant_file+" > "+eyeVarP_conf.sort_file1 # create a working input file 
#		print (COMMAND) # create a working input file 
		subprocess.call(COMMAND, shell=True) # do it in shell
		lr=-1 # start with header line
		for line in open(eyeVarP_conf.sort_file1): 
			lr = lr+1 # count how many lines
#		pc_val=percentage[i]*100 #
		Outln=Outln+"Variants scoring "+str(int(pc_val))+" percentage - "+str(lr)+eyeVarP_conf.nl #
#
# creat a variants for all samples is at the percentage level - percentile_variant_file 
#	percentile=float(eyeVarP_conf.parameter_file)/100 # setup the percentile argument given as a floating point value
	percentile=int(int(top_score)*(float(eyeVarP_conf.parameter_file)/100)) # setup the percentile argument given as a floating point value
#	print(str(percentile)) #
#	COMMAND="awk '{ if ($1 >= "+str(percentile)+") print $0}' "+full_variant_file+" > "+percentile_variant_file # create a working input file for input percentile for all samples 
	COMMAND="awk '{ if ($2 >= "+str(percentile)+") print $0}' "+full_variant_file+" > "+percentile_variant_file # create a working input file for input percentile for all samples 
#	print (COMMAND) # create a working input file 
	subprocess.call(COMMAND, shell=True) # do it in shell
#	
	print(Outln)
	Output_file.write(Outln) # write the line to final output file 
	# Print out the top 20 variants 
	Outln=eyeVarP_conf.nl+"Top 20 variants in the VPOL :"+eyeVarP_conf.nl #
	with open(eyeVarP_conf.input_file,'r',encoding="utf-8") as input_file : # 
		for x in range(21):
			Outln=Outln+input_file.readline() #
	#print(Outln)
	Output_file.write(Outln) # write the line to final output file 
	
	#
###########################################################################################################
#
#
###########################################################################################################
def Samples_stats(Output_file,sample,pos): #
#
# input file for filtering is the output from the eyeVarP process. This means the variant gene name is in a column named GENE_NAME.
#
	global gene_loc #
#	
	Outln=Ast_ln+eyeVarP_conf.nl+"STATS for Sample "+sample+" : "+eyeVarP_conf.nl+eyeVarP_conf.nl #
	#
#
#	work gene and variant counts for this sample in the full scoring variant file
#
	COMMAND="awk '($"+str(pos)+" !=\"0\") {print $0}' "+full_variant_file+" > "+eyeVarP_conf.full_file2 # create a working input file 
#	print (COMMAND+"/"+str(pos)+"/"+sample)
	subprocess.call(COMMAND, shell=True) # do it in shell
	#
	lr=-1 # start with header line
	for line in open(eyeVarP_conf.full_file2): 
		lr = lr+1 #
	Outln=Outln+"Total number of scoring variants : "+str(lr)+eyeVarP_conf.nl #
#
# eyeVarP_conf.full_file2 - contains variants for all samples is at the percentile level
	unique_genes(eyeVarP_conf.full_file2) #
	lr=0 #
	for line in open(eyeVarP_conf.sort_file1): 
		lr = lr+1 #
	Outln=Outln+"Total number of genes with scoring variants : "+str(lr)+eyeVarP_conf.nl #
#
#  now list variants per gene for ones that meet the percentile criteria
#
	Outln=Outln+eyeVarP_conf.nl+"Variants within Gene breakdown based on supplied " #
	Outln=Outln+str(eyeVarP_conf.parameter_file)+" percentage parameter :"+eyeVarP_conf.nl # a supplied value
#
#	eyeVarP_conf.full_file1 - contains variants at the required percentiles for all samples #
#	process to create a subset of variants that occurs in the sample for the percentile entered eyeVarP_conf.full_file2 #
#
	COMMAND="awk '($"+str(pos)+" !=\"0\") {print $0}' "+percentile_variant_file+" > "+sample_percentile_variant_file # create a working input file 
	subprocess.call(COMMAND, shell=True) # do it in shell
	#
# Gene level count
	unique_genes(sample_percentile_variant_file) #
#
	for line in open(eyeVarP_conf.sort_file1): 
		lr = lr+1 #
		COMMAND="awk '($3 ~ \""+line.rstrip()+"\") {print $0}' "+sample_percentile_variant_file+" > "+eyeVarP_conf.full_file2 # create a working input file 
#		print (COMMAND+"/"+str(pos)+"/"+sample)
		subprocess.call(COMMAND, shell=True) # do it in shell
		lr = 0 #
		for vline in open(eyeVarP_conf.full_file2): 
			lr = lr+1 #
		Outln=Outln+eyeVarP_conf.nl+"Number of variants in "+line.rstrip()+" : "+str(lr) #
#
#	print(Outln) #
	Output_file.write(Outln) # write the line to final output file 
	# Print out the top 20 variants 
#	Outln=eyeVarP_conf.nl+"Top 20 variants in the VPOL for sample:"+eyeVarP_conf.nl #
#	with open(eyeVarP_conf.fullinput_file,'r',encoding="utf-8") as input_file : # 
#		for x in range(21):
#			Outln=Outln+input_file.readline() #
	#print(Outln)
#	Output_file.write(Outln) # write the line to final output file 
	#
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	global full_variant_file # variants with score
	global percentile_variant_file # variants that is above the percentile parameter
	global sample_percentile_variant_file # variants that is above the percentile parameter
	global gene_loc #
	global infile_shape # contains the array dimension of the input file, so you can get number of variants and number of samples
#
	eyeVarP_conf.init() #
	#
	print ("Variant Statistics - Main") #
	#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
	#
	COMMAND="cut --complement -f3-10 "+eyeVarP_conf.input_file+" > "+eyeVarP_conf.working_file1 # create a working input file - by reducing to only the cols we need 
	subprocess.call(COMMAND, shell=True) # do it in shell
#
# Now filter the input file by gene list 
	infile_shape = file_lines_input(eyeVarP_conf.working_file1) # set up the numpy array
#	print (infile_shape)
	gene_loc = 2 # standard col location for gene name
	with open(eyeVarP_conf.final_output_file,'w',encoding="utf-8") as Output_file : # 
#		print("do total variant "+eyeVarP_conf.nl)
		Total_variants_stats(Output_file) #
		i = eyeVarP_conf.txt_start+1 # where the samples start
		j = 0 #
		while (i <= eyeVarP_conf.txt_end): # for number of samples
#			print (i)
			Samples_stats(Output_file,eyeVarP_conf.Sample_ids[j],i) #
			i = i+1
			j = j+1
	#
#
#
############################################################################################################
