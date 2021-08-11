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
tab='\t' # 
nl='\n' #
#
Sample_loc=-1 #
#
info_msg1_1="eyeVarP samplef : Invalid number/type of inputs, :" 
info_msg1_2="eyeVarP samplef : 1) output destination directory + prefix" #
info_msg1_3="eyeVarP samplef : 2) Input file - output from eyeVarP prioritisation process" 
info_msg1_4="eyeVarP samplef : 3) Selection Criteria list (ped format)" #
info_msg1_5="eyeVarP samplef : 4)  optional - Sample to apply inheritance model filter (same as in ped file)" #
info_msg1_6="eyeVarP samplef : 5)  optional - Inheritance model - DN - DeNovo, AD - Autosomal Dominant, AR - Autosomal Recessive, CH - compound Hete" #
#
format_msg=info_msg1_1+nl+info_msg1_2+nl+info_msg1_3+nl+info_msg1_4+nl+info_msg1_5+nl+info_msg1_6 #
#
############################################################################################################
###########################################################################################################
#
###########################################################################################################
def initial_setup():
#
	global supplied_args 
	#
	print ("initial_setup():") #
	print (suffix) #
#	print (sys.argv) #
	supplied_args=len(sys.argv) #
	print (supplied_args) #
#
	if (supplied_args != 5 and supplied_args != 7 ):  # arg [0] is the python program
		print (format_msg) #
		return 1 #
	else :
		eyeVarP_conf.output_dir=sys.argv[2] #
		eyeVarP_conf.input_file=sys.argv[3] #
		eyeVarP_conf.selection_list=sys.argv[4] #
#		eyeVarP_conf.final_output_file=eyeVarP_conf.output_dir+"variant_filtered_output_file_"+suffix+".txt" #
		eyeVarP_conf.inh_model="NONE" #
		if (supplied_args == 7 ):
			eyeVarP_conf.inh_sample=sys.argv[5] #
			eyeVarP_conf.inh_model=sys.argv[6] #
			eyeVarP_conf.output_dir=sys.argv[2]+eyeVarP_conf.inh_sample+"_"+eyeVarP_conf.inh_model+"_" #
#			eyeVarP_conf.final_output_file=eyeVarP_conf.output_dir+"variant_filtered_output_file_"+eyeVarP_conf.inh_sample+"_"+eyeVarP_conf.inh_model+"_"+suffix+".txt" #
		eyeVarP_conf.final_output_file=eyeVarP_conf.output_dir+"variant_filtered_output_file_"+suffix+".txt" #
		print ("output : ",eyeVarP_conf.final_output_file) #
	
	return 0 #
#
#
###########################################################################################################
def setup_samples(): #
#
#	print "setup_samples(): " #
#
	eyeVarP_conf.Sample_ids=[] #
	with open(eyeVarP_conf.selection_list,'r',encoding="utf-8") as select_input: # 
		for line1 in select_input: # work each line of ped file 
			this_line=re.split('\t|\n|\r',line1,7) # split into sample id and status (1/0)
			this_line[6]=0 # 
			print (this_line) #
			this_line[5]=str(int(this_line[5])-1) # set on else leave as 2
#			print (this_line[5]) #
			eyeVarP_conf.Sample_ids.append(this_line) 
		#
#	print (eyeVarP_conf.Sample_ids) #
#
###########################################################################################################
def setup_working_ped(): # create a ped file suitable for the inheritance model we are looking for
#
	matched=False
	p1_matched=False
	p2_matched=False
	p1_id="" #
	p2_id="" #
#	print ("setup_working_ped(): ") #
#
	eyeVarP_conf.working_file1=eyeVarP_conf.output_dir+"working_file1_"+suffix+"_tmp.txt" # 
	eyeVarP_conf.working_file2=eyeVarP_conf.output_dir+"working_file2_"+suffix+"_tmp.txt" #
	with open(eyeVarP_conf.selection_list,'r',encoding="utf-8") as select_input, \
		open(eyeVarP_conf.working_file1,'w',encoding="utf-8") as temp_ped_file1, \
		open(eyeVarP_conf.working_file2,'w',encoding="utf-8") as temp_ped_file2 : # 
		for line1 in select_input: # work each line of ped file 
			this_line=re.split('\t|\n|\r',line1,7) # split into sample id and status (1/0)
#			print (this_line[1]) #
			if (this_line[1] == eyeVarP_conf.inh_sample ):  # is this the sample we want to apply model to?
				matched=True
				p1_id=this_line[2] #
				p2_id=this_line[3] #
				if (eyeVarP_conf.inh_model == "AR" ):  # for recessive?
					this_line[5]="3" # set as affected - homozygous reguired
				temp_ped_file1.write("\t".join(this_line)+nl) # write the line to final output file 
				temp_ped_file2.write("\t".join(this_line)+nl) # write the line to final output file 
			if (this_line[1] == p1_id ):  # is this a parental id?
				p1_matched=True
				if (eyeVarP_conf.inh_model == "CH" ):  # for compound hete?
					this_line[5]="2" # set as affected
					temp_ped_file1.write("\t".join(this_line)+nl) # write the line to final output file 
					this_line[5]="1" # set as unaffected
					temp_ped_file2.write("\t".join(this_line)+nl) # write the line to final output file 
				elif (eyeVarP_conf.inh_model == "AR" ):  # for recessive?
					this_line[5]="2" # set as affected
					temp_ped_file1.write("\t".join(this_line)+nl) # write the line to final output file 
				else : # any other model
					temp_ped_file1.write(line1) # write the line to final output file 
			if (this_line[1] == p2_id ):  # is this a parental id?
				p2_matched=True
				if (eyeVarP_conf.inh_model == "CH" ):  # for compound hete?
					this_line[5]="2" # set as affected
					temp_ped_file2.write("\t".join(this_line)+nl) # write the line to final output file 
					this_line[5]="1" # set as unaffected
					temp_ped_file1.write("\t".join(this_line)+nl) # write the line to final output file 
				elif (eyeVarP_conf.inh_model == "AR" ):  # for recessive?
					this_line[5]="2" # set as affected
					temp_ped_file1.write("\t".join(this_line)+nl) # write the line to final output file 
				else : # any other model
					temp_ped_file1.write(line1) # write the line to final output file 
		#
	if (matched and p1_matched and p2_matched ):  # all samples of the trio found
		return matched #
	else :
		return False # 
#
###########################################################################################################
def extract_the_variants(): #
#
	setup_samples() #
	filter_the_variants() #
#
###########################################################################################################
def CH_inh_model(): #
#
	eyeVarP_conf.temp_output_file=eyeVarP_conf.final_output_file # save final output file name
	# do one parent
	eyeVarP_conf.working_file3=eyeVarP_conf.output_dir+"working_file3_"+suffix+"_tmp.txt" # 
	eyeVarP_conf.sort_file1=eyeVarP_conf.output_dir+"sort_file1_"+suffix+"_tmp.txt" # 
	eyeVarP_conf.final_output_file=eyeVarP_conf.working_file3 # 
	eyeVarP_conf.selection_list=eyeVarP_conf.working_file1 #
	extract_the_variants() # for CH only pick heteozygous variants
	COMMAND="cut -f 11 "+eyeVarP_conf.working_file3+" | sort -u > "+eyeVarP_conf.sort_file1 #  
	subprocess.call(COMMAND, shell=True) # do it in shell
	# now the other
	eyeVarP_conf.working_file4=eyeVarP_conf.output_dir+"working_file4_"+suffix+"_tmp.txt" # 
	eyeVarP_conf.sort_file2=eyeVarP_conf.output_dir+"sort_file2_"+suffix+"_tmp.txt" # 
	eyeVarP_conf.final_output_file=eyeVarP_conf.working_file4 # 
	eyeVarP_conf.selection_list=eyeVarP_conf.working_file2 #
	extract_the_variants() # for CH only pick heteozygous variants
	COMMAND="cut -f 11 "+eyeVarP_conf.working_file4+" | sort -u > "+eyeVarP_conf.sort_file2 #  
	subprocess.call(COMMAND, shell=True) # do it in shell
	# merge the output together
	COMMAND="tail -n+2 "+eyeVarP_conf.working_file4+" >> "+eyeVarP_conf.working_file3 #  
	subprocess.call(COMMAND, shell=True) # do it in shell
	# merge also the gene name files
	eyeVarP_conf.sort_file3=eyeVarP_conf.output_dir+"sort_file3_"+suffix+"_tmp.txt" # 
	COMMAND="cat "+eyeVarP_conf.sort_file1+" "+eyeVarP_conf.sort_file2+" | sort | uniq -d > "+eyeVarP_conf.sort_file3 #  
	subprocess.call(COMMAND, shell=True) # do it in shell
	#need to sort it
	eyeVarP_conf.final_output_file=eyeVarP_conf.temp_output_file # restore final output file name
	eyeVarP_conf.temp_output_file=eyeVarP_conf.output_dir+"tmp_out_file_"+suffix+"_tmp.txt" # 
	COMMAND="sort -k 1,2r -k 11,11 "+eyeVarP_conf.working_file3+" > "+eyeVarP_conf.temp_output_file #  
	subprocess.call(COMMAND, shell=True) # do it in shell
#
###########################################################################################################
def filter_the_variants(): #
#
# input file for filtering is the output from the eyeVarP process. 
	global Sample_loc #
#
#	print "filter_the_variants(): " #
	#
	with open(eyeVarP_conf.input_file,'r',encoding="utf-8") as variants_file, open(eyeVarP_conf.final_output_file,'w',encoding="utf-8") as filtered_file : # 
		for line1 in variants_file: # work each line of new sample vcf file 
			write_it=False # initialise score 
			line_parts=re.split('\t|\n|\r',line1) # split the variant up
#			print "line part 0 : ",line_parts[0] #
			if ("#CHROM" != line_parts[2]): #
#				print src_line1 #
				write_it=filter_variants_by_Seg(line_parts) # check get sample values
				#
			else : # save the header line	
				write_it=True # initialise score 
				for i, content in enumerate(line_parts): # return the value and index number of the sample id item in the line array 
#					print "content-",content,"/",i				#
					for j in range(len(eyeVarP_conf.Sample_ids)): #
						if (content == eyeVarP_conf.Sample_ids[j][1]) : # look for sample id in header 
							eyeVarP_conf.Sample_ids[j][6]=i #save sample location
#						print "INFO_loc: ",INFO_loc #
					#	
				print (eyeVarP_conf.Sample_ids) #
			if (write_it): #
				filtered_file.write(line1) # write the line to final output file 
#
###########################################################################################################
#
###########################################################################################################
def filter_variants_by_Seg(INFO_details): #
#
#	print "filter_variants_by_Seg(INFO_details): " #
#
	val=True #
	for j in range(len(eyeVarP_conf.Sample_ids)): #
		if ( eyeVarP_conf.Sample_ids[j][6] != 0 ) : # check if this samples is in the input file 
			working_value=INFO_details[eyeVarP_conf.Sample_ids[j][6]] # copy the sample's genotype value for the variant
# match it with QC failed values
#			if ( working_value=="2" and eyeVarP_conf.inh_model != "AR" ) : # if homozygous and this is not AR, then set the sample's variant value as hete to match the sample's ped value
			if ( working_value=="2" and eyeVarP_conf.inh_model != "AR" and eyeVarP_conf.inh_model != "CH" ) : # if homozygous and this is not AR/CH, then set the sample's variant value as hete to match the sample's ped value
				working_value="1" #
#			if ( eyeVarP_conf.Sample_ids[j][5] != INFO_details[eyeVarP_conf.Sample_ids[j][6]] ) : # yes - check the sample's value 
#			if ( eyeVarP_conf.Sample_ids[j][5] != working_value ) : # yes - check the sample's value 
#			QC_working_value=working_value # setup QC
#			if ( QC_working_value=="9" ) : # if homozygous and this is not AR, then set it as hete 
#				working_value="2" #
#			elif ( QC_working_value=="8" ) : # if homozygous and this is not AR, then set it as hete 
#				QC_working_value="1" #
#			print (eyeVarP_conf.Sample_ids[j][5],"/",INFO_details[eyeVarP_conf.Sample_ids[j][6]],"/",working_value,"/",QC_working_value)
#			if (( eyeVarP_conf.Sample_ids[j][5] != working_value ) and ( eyeVarP_conf.Sample_ids[j][5] != QC_working_value )) : # yes - check the sample's ped value against the sample's VPLOL value  - test
			if ( eyeVarP_conf.Sample_ids[j][5] != working_value ) : # yes - check the sample's ped value against the sample's VPLOL value  - test
				val=False # not what is needed 
				break # then get out and move to next variant
#
	return val #
#	
##
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	eyeVarP_conf.init() #
	#
	print ("Variant Sample and Inheritance Model Filter - Main") #
	#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
	#
# Now filter the input file by gene list 
	if (eyeVarP_conf.inh_model=="NONE"): #
		extract_the_variants() #
	elif (eyeVarP_conf.inh_model in eyeVarP_conf.Inheritance_model) : # check inheritance model
#		print (" inheritance model found") 
		if not setup_working_ped() :  #
			print ("Sample ID or Parent IDs not in ped file") #
		else : #
			if (eyeVarP_conf.inh_model in ["CH"]) : #
				CH_inh_model() #
			else : # for DN,AR,AD 
				eyeVarP_conf.selection_list=eyeVarP_conf.working_file1 # setup the temp ped file
				extract_the_variants() #
	elif (eyeVarP_conf.inh_model=="SELECT") : # subset selected samples from VPOL input
#		print (" select samples found") 
		if not setup_working_ped() :  #
			print ("Sample ID or Parent IDs not in ped file") #
	else :
		print (format_msg) #
#
#	print (eyeVarP_conf.Sample_ids) #
	#
#	clean_up() #
#
############################################################################################################
