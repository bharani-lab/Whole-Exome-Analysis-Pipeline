###########################################################################################################
# #
###########################################################################################################
#
import logging, sys, re, glob, os, subprocess, time #
import numpy as np #
import eyeVarP_conf #
import gzip #
from shutil import copyfile #
#from __main__ import *
tab='\t' # 
nl='\n' #
#logging.basicConfig(stream=sys.stderr, level=logging.DEBUG) # to show debug message remove '#' at the beginning of line
#
###########################################################################################################
#
###########################################################################################################
def parameters_p1(input_file):
#
#	global pop_array 
#	global pred_array 
#
#		print ("parameters_p1") 
		predictors=False #
		while True:
			line=input_file.readline() # go thru the INFO annotation section of the VCF header
#		print line 
			this_line=re.split('\t|\n|\r|=|,',line) #
#			print (this_line) 
			if (this_line[0] =="#CHROM") : 
				break #
			elif (this_line[2] =="ANNOVAR_DATE") : #
				predictors=True #
			elif (this_line[2] =="ALLELE_END") : #
				predictors=False #
#
			if (predictors and this_line[2] !="ANNOVAR_DATE") : # true
#				print this_line[2] #
#				print eyeVarP_conf.Population # when filtering for QC value 
				if any( s in this_line[2] for s in eyeVarP_conf.Population ) : # when filtering for QC value 
					eyeVarP_conf.pop_array.append([eyeVarP_conf.PF,this_line[2],eyeVarP_conf.MAFval,"LE"]) 
				else :
					eyeVarP_conf.pred_array.append([eyeVarP_conf.PD,this_line[2],"N","","0","","1","","2"]) 
#
###########################################################################################################
#
###########################################################################################################
def parameters(input_file):
#
#	global pop_array 
#	global pred_array 
#
	print ("eyeVarP_VCF.parameters: ") #
#
	#
	if input_file.endswith('.gz'):
		with gzip.open(input_file,'rt') as first_fn : #
			parameters_p1(first_fn) #
	else:
		with open(input_file,'r',encoding="utf-8") as first_fn : #
			parameters_p1(first_fn) #
###########################################################################################################
#
###########################################################################################################
def setup_default_pred_values_1(first_fn): #
#
#	global pred_array 
#
	for line in first_fn: #
		this_line=re.split('\t|\n|\r|=|;',line) #
#			print line 
		try: #
			for i, content in enumerate(this_line): # return the value and index number of each item in the line array 
#				print "content-",content,"/",i				#
				for j in range(len(eyeVarP_conf.pred_array)): #
					if (content == eyeVarP_conf.pred_array[j][1]) : # when filtering for QC value 
#						print eyeVarP_conf.pred_array[j][1], content, this_line[i+1] #
						#if ((this_line[i+1] != ".") and (this_line[i+1] != "-999")): 
						if ((this_line[i+1] != ".") and (this_line[i+1] != "-") and (this_line[i+1] != "-999")): # numeric		
#							if (eyeVarP_conf.is_number(this_line[i+1])): # numeric
							if (eyeVarP_conf.is_number(this_line[i+1]) and (eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_type] != "A" )): # numeric
								this_numeric=True
								if ((eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low] == "") or (float(eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low]) > float(this_line[i+1]))): 
									eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low] = this_line[i+1] 
								if ((eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high] == "") or (float(eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high]) < float(this_line[i+1]))): 
									eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high] = this_line[i+1] 
							else : # if not a value that can be expressed as a floating point then is alpha or field is alpha
#								print "except", this_line[i] 
								this_numeric=False
								if ((eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low] == "") or (eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low] > this_line[i+1])): 
									eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low] = this_line[i+1] 
								if ((eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high] == "") or (eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high] < this_line[i+1])): 
									eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high] = this_line[i+1] 
								eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_type] = "A" 
								aa=len(eyeVarP_conf.pred_array[j]) 
#								print aa #
								if (aa < eyeVarP_conf.Maxval): #still can add some more alpha options
									k=eyeVarP_conf.startlen # point to 1st option slot 
#									print k, aa #
									while (k <= aa) :
										if (k == aa): # completed search and not found	
											eyeVarP_conf.pred_array[j].append(this_line[i+1]) # then add it
#											print "add- ",this_line[i+1],"/",this_line[i] 
										else : # keep searching
											if (eyeVarP_conf.pred_array[j][k] == this_line[i+1] ) : # current value already added 
												break #
										k+=1 # move to eyeVarP_conf.pred_array slot
#									print eyeVarP_conf.pred_array[j] #
		except: # debug messages
			print("Error occurred at line :", line) #
			print("pred value in line :", this_line[i+1],":", this_numeric) #
#			print("pred_array when error occurred :", eyeVarP_conf.pred_array) #
#			print("pred_array index when error occurred :", j) #
			print("pred_array when error occurred :", eyeVarP_conf.pred_array[j]) #
			sys.exit(1) #
#
###########################################################################################################
#
###########################################################################################################
def setup_default_pred_values(file1): #
#
#	global pred_array 
#
	print ("eyeVarP_VCF.setup_default_pred_values: ") #
##
	if file1.endswith('.gz'):
		print ("using input file: ", file1input_file) #
		with gzip.open(file1input_file,'rt') as first_fn : #
			setup_default_pred_values_1(first_fn) #
	else:
		print ("using input file: ", file1) #
		with open(file1,'r',encoding="utf-8") as first_fn : #
			setup_default_pred_values_1(first_fn) #
#
#	print pop_array #
#	print eyeVarP_conf.pred_array #
	#
# determine a middle value for the numeric predictors 
	for j in range(len(eyeVarP_conf.pred_array)): #
#		print eyeVarP_conf.pred_array[j] 
		if (eyeVarP_conf.pred_array[j][1] == "N") : # this predictor is numeric 
			if (eyeVarP_conf.is_number(eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high]) and eyeVarP_conf.is_number(eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low])): #
				mid_value=float((float(eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_high])+float(eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_low]))/2) #
#			print "mid-",mid_value #
				eyeVarP_conf.pred_array[j][eyeVarP_conf.PD_mid] = mid_value 
#
#	print eyeVarP_conf.pred_array 

#
###########################################################################################################
#
###########################################################################################################
def read_variant_source_file(): #
#
##
#	print "read_variant_source_file_for_VCF(): " #
#	print info_msg2 # parameter file supplied then use it
	print (eyeVarP_conf.parameter_file) #
	#
	with open(eyeVarP_conf.input_file,'r',encoding="utf-8") as input: # 
		for line in input: # work each input vcf file 
			this_line=re.split('\t|\n|\r',line) # split into file location and sample id
#			print this_line #
			if ( setup_for_this_src_file(this_line) != 0 ): #
				print ("Issue with input file 1: ", line) #
				return 1 # error

			else : #
				print ("processing input file : ", this_line) # 
				if ((eyeVarP_conf.sample_loc >=0) and (eyeVarP_conf.INFO_loc >=0) and (eyeVarP_conf.FORMAT_loc >=0) and (eyeVarP_conf.sample_coverage_loc[0] >=0)): #
#					print "OK" #
					eyeVarP_conf.header_ln = eyeVarP_conf.header_ln+tab+eyeVarP_conf.sample_ID # add sample ID to header line
#					print "header- ",header_ln #
#					outline1 = eyeVarP_conf.header_ln+nl #
					eyeVarP_conf.blank_variant_ln = eyeVarP_conf.blank_variant_ln+tab+"0" #
#					print "blank_variant_ln- ",blank_variant_ln #
#					print sample_coverage, "/", sample_coverage_loc #
					work_this_src_file(this_line)  #
					if (os.path.isfile(eyeVarP_conf.full_file2) is False ) or (os.stat(eyeVarP_conf.full_file2).st_size == 0): # are there any variants in the final file? 
						eyeVarP_conf.update_existing_variants(eyeVarP_conf.working_file1,"1") # No - then create it based on sample file 
					elif (os.stat(eyeVarP_conf.working_file1).st_size == 0): # are there any variants for this sample? 
						copyfile(eyeVarP_conf.full_file2,eyeVarP_conf.working_file1) # copy the new full file
						eyeVarP_conf.update_existing_variants(eyeVarP_conf.working_file1,"0") # then for existing variants, just add 0 to the end for this sample 
					else : # yes, then incorporate them 
						eyeVarP_conf.incorporate_this_src_into_full_file() #
					copyfile(eyeVarP_conf.full_file2,eyeVarP_conf.full_file1) # copy the new full file
					
					
				else : #
					print ("Issue with input file 2: ", line) #
					return 1 # error

###########################################################################################################
#
###########################################################################################################
def setup_for_this_src_file_1(source_vcf,file_line): #
#
	for src_line in source_vcf: # work each line of source vcf file 
#			print (src_line) #
		src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
		if ("##" not in src_line1[0]): #
#				print src_line1 #
			if ("#CHROM" in src_line1[0]): # find sample location
				for i, content in enumerate(src_line1): # return the value and index number of each item in the line array 
#								print "content-",content,"/",i				#
					if (content == file_line[1]) : # when filtering for sample ID 
						eyeVarP_conf.sample_loc=i #save sample location
						eyeVarP_conf.sample_ID=file_line[1] # save sample ID
#							print "SAMPLE: ",sample_loc #
					elif (content == "INFO") : # 
						eyeVarP_conf.INFO_loc=i #save sample location
#							print "INFO_loc: ",INFO_loc #
					elif (content == "FORMAT") : # 
						eyeVarP_conf.FORMAT_loc=i #save sample location
#							print "FORMAT_loc: ",FORMAT_loc #
					#	
				if ("#CHROM" not in eyeVarP_conf.header_ln): # is this first time
#						print "setup header" #
					eyeVarP_conf.header_ln='\t'.join(src_line1[:eyeVarP_conf.INFO_loc+1])+tab+eyeVarP_conf.GENE_NAME # save heading up to INFO (using INFO_loc +1 to denote the stop field)
					eyeVarP_conf.blank_variant_ln=eyeVarP_conf.header_ln # save heading up to INFO (using INFO_loc +1 to denote the stop field)
#						eyeVarP_conf.header_ln='\t'.join(src_line1[:eyeVarP_conf.INFO_loc+1]) # save heading up to INFO (using INFO_loc +1 to denote the stop field)
#						eyeVarP_conf.blank_variant_ln='\t'.join(src_line1[:eyeVarP_conf.INFO_loc+1]) #
			else : # variants lines 
				FORMAT1=re.split(':',src_line1[eyeVarP_conf.FORMAT_loc]) # split into file location and sample id
#				print ("format : ",FORMAT1) #
				for i, content in enumerate(FORMAT1): # return the value and index number of each item in the line array 
#					print ("content-",content,"/",i)				#
					for j in range(len(eyeVarP_conf.sample_coverage)): #
#						print (eyeVarP_conf.sample_coverage[j],"/",content) #
						if (content == eyeVarP_conf.sample_coverage[j]) : 	# look for FORMAT field 
							eyeVarP_conf.sample_coverage_loc[j]=i 			# save sample location	
							break #
#					print ("coverage : ",eyeVarP_conf.sample_coverage_loc,i) #
				source_vcf.close() # finish with source vcf file 
				return 0 # have setup all location values - ok to go back
#							print INFO1 #
#							print FORMAT1 #
#							print SAMPLE1 #

#ed	source_vcf.close() # finish with source vcf file 
	return 1 # not a clean exit
#
###########################################################################################################
#
###########################################################################################################
def setup_for_this_src_file(file_line): #
#
#	global GT_loc #
	
#
#	print "setup_for_this_src_file): ", file_line #
#	print info_msg2 # parameter file supplied then use it
	eyeVarP_conf.sample_loc=-1 #
	eyeVarP_conf.sample_ID=""
	eyeVarP_conf.INFO_loc=-1 #
	eyeVarP_conf.FORMAT_loc=-1 #
#	eyeVarP_conf.sample_coverage_loc = [-1,-1,-1,-1,-1] # location of VCF format codes for sample 
#	eyeVarP_conf.sample_coverage_loc = [-1,-1,-1,-1,-1,-1] # location of VCF format codes for sample 
	#
#	with open(file_line[0],'r', encoding="utf-8") as source_vcf : #
##
	return_val=0 #
	if file_line[0].endswith('.gz'):
		with gzip.open(file_line[0],'rt') as source_vcf : #
			if ( setup_for_this_src_file_1(source_vcf,file_line) != 0 ): #
				return_val=1 #
	else:
		with open(file_line[0],'r',encoding="utf-8") as source_vcf : #
			if ( setup_for_this_src_file_1(source_vcf,file_line) != 0 ): #
				return_val=1 #
#
	return return_val #
#
###########################################################################################################
#
###########################################################################################################
def work_this_src_file_1(source_vcf, wrkf1): #
#
##
#	print "work_this_src_file(file_line): " #
#	print working_file1 #
		for src_line in source_vcf: # work each line of source vcf file 
#			print (src_line) #
			src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
#			if ("#" not in src_line1[0]): # skip the header lines
			if (("#" not in src_line1[0]) and ("_" not in src_line1[0])): # skip the header lines and any alternate chromosome alignments variants
#				print (src_line1) #
#				print (src_line1[eyeVarP_conf.FORMAT_loc]) #
#				print (src_line1[eyeVarP_conf.sample_loc]) #
# variants lines 
				Sample_coverage=0 # initial sample coverage total 
				Sample_GQ=0 # initial sample genotype quality score 
				SAMPLE1_F=re.split(':',src_line1[eyeVarP_conf.FORMAT_loc]) # split the sample's FORMAT fields 
				eyeVarP_conf.sample_coverage_loc = [-1,-1,-1,-1,-1,-1]
				for j in range(len(SAMPLE1_F)) : #
#						print (SAMPLE1_F[j]) #
						for k in range(len(eyeVarP_conf.sample_coverage)) : #
							if (SAMPLE1_F[j] == eyeVarP_conf.sample_coverage[k]): # VCF format codes to look out for 
#								print (SAMPLE1_F[j], "/", eyeVarP_conf.sample_coverage[k], "/", j) #
								eyeVarP_conf.sample_coverage_loc[k]=j 
								break #
				SAMPLE1=re.split(':',src_line1[eyeVarP_conf.sample_loc]) # split the sample's FORMAT fields 
#				print sample_coverage_loc[1],"/", sample_coverage_loc[2] #
#				print SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]],"/", SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]] #
# have a coverage check on the sample 
#				print ("MaxCOverage : ",eyeVarP_conf.Maxcoverage) #
#				print ("Hete_balance : ",eyeVarP_conf.Hete_Balance) #
#				print ("Genotype Quality : ",eyeVarP_conf.Genotype_Quality) #
#				print "coverage_loc : ",eyeVarP_conf.sample_coverage_loc #
#				print (SAMPLE1_F) 
#				print (SAMPLE1) 
#				print ("GT: ",eyeVarP_conf.GT_val,"/ NR: ",eyeVarP_conf.NR_val,"/ NV: ",eyeVarP_conf.NV_val,"/ DP: ",eyeVarP_conf.DP_val,"/ AD: ",eyeVarP_conf.AD_val,"/ GQ: ",eyeVarP_conf.GQ_val) 
#				print ("GT: ",eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val],":",SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]]) 
#				print ("NR: ",eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val],"/",SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]] )
#				print ("NV: ",eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val],"/",SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val]] )
#				print ("DP: ",eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val],"/",SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]]) 
#				print ("AD: ",eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.AD_val],"/",SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.AD_val]]) 
#				print ("GQ: ",eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GQ_val],"/",SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GQ_val]]) 
# check to see if we want this variant line for this sample
				if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]] != "./.") and (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]] != "0/.") and (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]] != "./0") and (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]] != "0/0") : # a valid alternate genotype 
#					print ("DP: ",eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]) 
#    				check coverage depth
					if (eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val] != -1 ) : #this sample have a NR_value
						if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]] == ".") : # no NR_val 
							SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]] = "0"  # set it as zero 
					if (eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val] != -1 ) : #this sample have a NV_value
						if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val]] == ".") : # no NV_val 
							SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val]] = "0"  # set it as zero 
					if (eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val] != -1 ) : #this sample have a coverage depth
						if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]] == ".") : # no DP_val 
							SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]] = "0"  # set it as zero 
						Sample_coverage=int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]]) # save DP value
						Alt_reads=int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]])/2 # save DP value
#						print ("DP") #
					if (eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.AD_val] != -1 ) : #this sample have an allele depth
						if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.AD_val]] != ".") : # has AD_val 
							AD_values=re.split(',',SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.AD_val]]) # get the alleles depth
							if (len(AD_values) > 1 ) : #  
								Alt_reads=int(AD_values[1]) # save alternate read count
							else : # this AD only has alternate count
								Alt_reads=int(AD_values[0]) # save alternate read count
#						print ("AD") #
					if (eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GQ_val] != -1 ) : #this sample have a genotype quality score
						if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GQ_val]] == ".") : # no GQ_val 
							SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GQ_val]] = "0"  # set it as zero 
						Sample_GQ=int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GQ_val]]) # save GQ value
#						print ("GQ") #
					if (eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val] != -1 ) and (eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val] != -1 ) : #this sample have a coverage depth from NR and NV
						if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]] == ".") : # no NR_val 
							SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]] = "0"  # set it as zero 
						if (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val]] == ".") : # no NV_val 
							SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val]] = "0"  # set it as zero 
#						Sample_coverage=int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]])+int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val]]) # save DP value
						Sample_coverage=int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]]) # save DP value
						Alt_reads=int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NV_val]]) # save alt reads value
#						print ("NR+NV") #
#					print ("TOT: ",str(Sample_coverage)) # 
#					print ("ALT_READ: ",str(Alt_reads)) # 
#					print ("pass") #
#					if (((eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val] == -1) or 
#						((eyeVarP_conf.is_number(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]])) and (int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.NR_val]]) >= int(eyeVarP_conf.Maxcoverage)))) and  #  NR
#						((eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val] == -1) or 
#						((eyeVarP_conf.is_number(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]])) and (int(SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.DP_val]]) >= int(eyeVarP_conf.Maxcoverage))))) : # DP
					eyeVarP_conf.QC_PASS=False #
#	QC checks
#					if (Sample_coverage >= int(eyeVarP_conf.Maxcoverage)) : 
					if ((Sample_coverage >= int(eyeVarP_conf.Maxcoverage)) and (Sample_GQ >= int(eyeVarP_conf.Genotype_Quality))) : 
						if ( Sample_coverage == 0 ) : # can't have divide by 0
							Sample_coverage=1 # set a dummay value
						if (int((Alt_reads/Sample_coverage)*100) >= int(eyeVarP_conf.Hete_Balance)) : # Pas QC for coverage and balance
							eyeVarP_conf.QC_PASS=True # Yes
#							print ("QC_PASS",eyeVarP_conf.QC_PASS) #
					GT_values=re.split('/',SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]]) # get the genotype fields
#				print GT_values #
					for j in range(len(GT_values)) : #
#					print (GT_values[j]) #
						if ( GT_values[j] not in eyeVarP_conf.Non_alt_GT_types ) : # when filtering for QC value 
#							print ("keep this variant1") #
							check_this_variant(src_line, wrkf1) #
							break # get out of for loop (GT_values)
#				else : #  NR
#					print ("not adding") #
#
#	source_vcf.close() # finish with source vcf file 
###########################################################################################################
#
###########################################################################################################
#
def work_this_src_file(file_line): #
#
##
#	print "work_this_src_file(file_line): " #
#	print working_file1 #
	wrkf1=open(eyeVarP_conf.working_file1,'w', encoding="utf-8") # 
#
	if file_line[0].endswith('.gz'):
		with gzip.open(file_line[0],'rt') as source_vcf : #
			work_this_src_file_1(source_vcf,wrkf1) #
	else:
		with open(file_line[0],'r',encoding="utf-8") as source_vcf : #
			work_this_src_file_1(source_vcf,wrkf1) #
#
	wrkf1.close() # finish with the output file 
#
#	print "sort unique" #
	COMMAND="sort -V -u -k1,5 "+eyeVarP_conf.working_file1+" > "+eyeVarP_conf.sort_file1 #  
	subprocess.call(COMMAND, shell=True) #
	copyfile(eyeVarP_conf.sort_file1,eyeVarP_conf.working_file1) # copy back

#	source_vcf.close() # finish with source vcf file 
#
###########################################################################################################
#
###########################################################################################################
def check_this_variant(src_line, wrkf1):  #
#
#	print "check_this_variant(src_line, wrkf1):  ",src_line #
#
	src_line1=re.split('\t|\n|\r',src_line) # split into file location and sample id
	SAMPLE1=re.split(':',src_line1[eyeVarP_conf.sample_loc]) # split the sample's FORMAT fields 
	GENOTYPE1=re.split('/',SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]]) # split the sample's GENOTYPE fields 
#	print (src_line1) #
#	print (SAMPLE1) #
#	print (GENOTYPE1) #
	if (population_frequency(src_line1[eyeVarP_conf.INFO_loc]) == 0 ): # check if variant with in filter
#		print ("saving this variant1") #
		eyeVarP_conf.gene_ref="NONE" #
		find_gene_ref(src_line1[eyeVarP_conf.INFO_loc])    # check if variant found != means yes 
#
		if (not eyeVarP_conf.QC_PASS) : # sample did not pass QC  
			gtype="." #
#		elif (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]] == "1/1") : # a homozygous alt genotype 
		elif (GENOTYPE1[0] == GENOTYPE1[1]) : # a homozygous alt genotype 
			gtype="2" #
		else : #
			gtype="1" 
#		print (SAMPLE1[eyeVarP_conf.sample_coverage_loc[eyeVarP_conf.GT_val]],"-",gtype)#
		outline='\t'.join(src_line1[:eyeVarP_conf.INFO_loc+1])+tab+eyeVarP_conf.gene_ref+tab+gtype+nl #
		wrkf1.write(outline) #
#	else : #
#		print "did not pass PF" #
#
###########################################################################################################
#
###########################################################################################################
def population_frequency(info_ln): #
#
#
#	print "population_frequency(info_ln): " #
	#
	val=0 #
	INFO1=re.split(';|=',info_ln) # split into file location and sample id
#	print (INFO1) #
#	print (eyeVarP_conf.PF_array) #
	for j in range(len(eyeVarP_conf.PF_array)): #
#		print (eyeVarP_conf.PF_array[j][1]) # move to pred_array slot
		for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#			print ("content-",content,"/",i)				#
#			print (eyeVarP_conf.PF_array[j][1]) # move to pred_array slot
			if (content == eyeVarP_conf.PF_array[j][1]) : # the population freq annotation we want? 
				if ( INFO1[i+1] != ".") :	 # if numeric check, else ok as most likely "." to state no annotation 
					temp_val=INFO1[i+1]
				else : # not number, so set a default as no annotation in population_frequency check would mean novel	
					temp_val=-999
#				print ("temp_val :",temp_val) #
#				if ( eyeVarP_conf.is_number(INFO1[i+1]) ) : # numeric
				#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
#					if ( float(INFO1[i+1]) > float(eyeVarP_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
				if ( eyeVarP_conf.is_number(temp_val) ) : # numeric
				#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
#					print ("temp_val yes :",temp_val) #
#					print (eyeVarP_conf.PF_array[j][3]) #
					if ( eyeVarP_conf.PF_array[j][3] == "LE" ) :	 # lower or each to PF limit 
						if ( float(temp_val) > float(eyeVarP_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
							val=1 # do not want this variant 
							break #
					else : # GT  look for greater or equal to PF limit	
						if ( float(temp_val) < float(eyeVarP_conf.PF_array[j][2]) ) :	 # when number is < 0.0001 it is expressed as e-0x 
#					print INFO1[i],INFO1[i+1],"/",PF_array[j][2] #
							val=1 # do not want this variant 
							break #
	#
#	print ("val :",val)
	return val #
#	
###########################################################################################################
#
###########################################################################################################
def find_gene_ref(info_ln): #
#
#
#	print "find_gene_ref(info_ln): " #
	#
	val=0 #
	INFO1=re.split(';|=',info_ln) # split into file location and sample id
#	print INFO1 #
#	print PF_array #
#	print PF_array[j][1] # move to pred_array slot
	for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#		print "content-",content,"/",i				#
		if (content == eyeVarP_conf.GN_value) : # the population freq annotation we want? 
			eyeVarP_conf.gene_ref=INFO1[i+1] # save gene name 
			break #
#
###########################################################################################################
#
###########################################################################################################
def score_the_variants(): #
##
#	print "score_the_variants(): " #
	#
	with open(eyeVarP_conf.full_file1,'r',encoding="utf-8") as variants_file, open(eyeVarP_conf.working_file1,'w',encoding="utf-8") as score_file : # 
		for line1 in variants_file: # work each line of new sample vcf file 
			priority_score=0 # initialise score 
			line_parts=re.split('\t|\n|\r',line1) # split the variant up
#			print "line part 0 : ",line_parts[0] #
			if ("#CHROM" != line_parts[0]): #
				logging.debug('variant line : %s', str(line1))
				#print ("variant line : ",line1) #
#				
				if (len(eyeVarP_conf.PD_array) > 0 ) : #  
					priority_score=prioritise_variants_by_predictors(line_parts[eyeVarP_conf.INFO_loc]) # get priority score for the PDs in variant
				#
				if (len(eyeVarP_conf.VT_array) > 0 ) : #   
#				print "check VT" #
					priority_score=priority_score+prioritise_variants_by_VT_types(line_parts[eyeVarP_conf.INFO_loc]) # get priority score for the VT in the variant
				#
				if ( float(priority_score) <= 0 ) : # check priority score
					priority_score=0 #
				#
				if ( int(priority_score) >= int(eyeVarP_conf.VariantScoreThreshold) ) : # check priority score,is it larger than threshold
					outline=str(priority_score)+tab+line1 # yes - then output it
					score_file.write(outline) #
#				print outline #
#				score_file.write(outline) #
			else : # save the header line	
				outline = "Priority_score"+tab+line1 #
				score_file.write(outline) #
#
###########################################################################################################
#
###########################################################################################################
def prioritise_variants_by_predictors(INFO_details): #
#
#	global PD_array 
#
#	print "prioritise_variants_by_predictors(INFO_details): " #
#
	type=2 #
	a_value=a_score=0 #
#
	val=wrkval=pdval=0 # this is the predictor priority score
	INFO1=re.split(';|=',INFO_details) # split into file location and sample id
#	print ("INFO1 - ",INFO1) #
	for j in range(len(eyeVarP_conf.PD_array)): #
		for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#			print ("content-",content,"/",i)				#
			if (content == eyeVarP_conf.PD_array[j][1]) : # the predictor annotation we want? 
				limit=len(eyeVarP_conf.PD_array[j]) # 
#				print ("PD_array: ", limit, eyeVarP_conf.PD_array[j]) #
#				print ("PD_array: ", limit, eyeVarP_conf.PD_array[j]) #
				logging.debug('PD INFO1 : %s / %s', INFO1[i], INFO1[i+1])
#				print ("PD INFO1: ", INFO1[i], INFO1[i+1]) #
				a_value=type+1 #
				a_score=a_value+1 #
				pdval=0 #
				predval=re.split('\||&',INFO1[i+1]) # split into individual values using delimiter | and & 
#				print ("predval : ",predval) #
#
				for k in range(len(predval)): #
					b_value=a_value #
					b_score=a_score #
					wrkval=0 #
#
					if (eyeVarP_conf.PD_array[j][type] == "A") : # character field 
						while (b_value < len(eyeVarP_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
#							print ("compare A-PD : vcf annotation - ", predval[k], "PPF val - ",eyeVarP_conf.PD_array[j][b_value]) #
							if ( predval[k] == eyeVarP_conf.PD_array[j][b_value] ) :		# is this the variant type we are looking for   
								wrkval=int(eyeVarP_conf.PD_array[j][b_score])	 		# yes - return the score for the variant 
								break 										# finish 
							b_value+=2	 # bump to next value option   
							b_score+=2	 #   
						if (wrkval > pdval ) : #check for higher predictor score
							pdval=wrkval	
#						while (a_value < len(eyeVarP_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
#							if ( INFO1[i+1] == eyeVarP_conf.PD_array[j][a_value] ) :		# is this the variant type we are looking for   
#								val=val+int(eyeVarP_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
#								break 										# finish 
#							a_value+=2	 # bump to next value option   
#							a_score+=2	 #   
					else : # numeric field 
#						if ( INFO1[i+1] != ".") :	 # if annotation then continue   
						if ( eyeVarP_conf.is_number(predval[k]) ) : # numeric
							while (b_value < len(eyeVarP_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
#								print ("compare N-PD : vcf annotation - ", predval[k], "PPF val - ",eyeVarP_conf.PD_array[j][b_value]) #
								if ( float(predval[k]) < float(eyeVarP_conf.PD_array[j][b_value])) :	 # is this the variant type we are looking for   
									wrkval=int(eyeVarP_conf.PD_array[j][b_score])	 		# yes - return the score for the variant 
									break 										# finish 
								if ( b_value+2 >= len(eyeVarP_conf.PD_array[j])-1) : # end of this PD limits  
#									print ("checking last PD:") #
									if ( float(predval[k]) > float(eyeVarP_conf.PD_array[j][b_value])) :	 # check the last value, which is a greater then check   
										wrkval=int(eyeVarP_conf.PD_array[j][b_score])	 		# yes - return the score for the variant 
								b_value+=2	 # bump to next value option   
								b_score+=2	 # 
							if (wrkval > pdval ) : #check for higher predictor score
								pdval=wrkval	
#						if ( eyeVarP_conf.is_number(INFO1[i+1]) ) : # numeric
#							print INFO1[i+1] #
#							while (a_value < len(eyeVarP_conf.PD_array[j])-1) : # loop thru the options for predictor scores 
#								print eyeVarP_conf.PD_array[j][a_value] #
	#							print "INFO1:", INFO1[i+1] #
#								if ( float(INFO1[i+1]) < float(eyeVarP_conf.PD_array[j][a_value])) :	 # is this the variant type we are looking for   
#									val=val+int(eyeVarP_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
#									break 										# finish 
	#							work1=a_value+2 # end of this PD limits  
#								print "work1:", work1 # end of this PD limits  
#								if ( a_value+2 >= len(eyeVarP_conf.PD_array[j])-1) : # end of this PD limits  
#									print ("checking last value") #
#									if ( float(INFO1[i+1]) > float(eyeVarP_conf.PD_array[j][a_value])) :	 # check the last value, which is a greater then check   
#										val=val+int(eyeVarP_conf.PD_array[j][a_score])	 		# yes - return the score for the variant 
#								a_value+=2	 # bump to next value option   
#								a_score+=2	 # 
#							
	#						print "a_score:", a_score #
#							print "len:", len(eyeVarP_conf.PD_array[j]) #
#				print ("variant's PD score : ",pdval) #
				val=val+pdval	 		# add to overall variant score 
				break # done for this predictor annotation 
#
	return val #
#	
##
#	print info_msg2 # parameter file supplied then use it
###########################################################################################################
#
###########################################################################################################
def prioritise_variants_by_VT_types(INFO_details): #
#
#	global VT_array 
	variant_INFO_array=[]
	variant_INFO_val_array=[]
#
#	print "prioritise_variants_by_VT_types(INFO_details): " #
#	print VT_array #
	val=wrkval=pdval=-999 # this is the predictor priority score
	VT_already_stored=False# 
	VT_found_for_variant=False #
	VT_k_found=False #
	VT_k_has_annotation=False #
	VT_has_annotation_value_in_variant=False #
	INFO1=re.split(';|=',INFO_details) # split into field name and values
#	print ("Variant INFO1 :",INFO1) #
	VT_array_size=len(eyeVarP_conf.VT_array)
#	print ("number of VT lines in PPF :",VT_array_size) #
#
#   loop thru the PPF VT lines
	for j in range(len(eyeVarP_conf.VT_array)): # move thru the VT lines in PPF
#		print ("VT array :", eyeVarP_conf.VT_array[j]) # move to pred_array slot
#
# 		loop thru the variant's INFO details, working one VT PPF line at a time
		for i, content in enumerate(INFO1): # return the value and index number of each item in the line array 
#			print ("content-",content,"/",i)				#
			if (content == eyeVarP_conf.VT_array[j][1]) : # the variant functional annotation we want? then save the annotation details for later use
				VT_already_stored=False# 
				if (len(variant_INFO_array) > 0) : # have I saved something already- yes
					for v in range(len(variant_INFO_array)): # move thru the saved VT INFO fields
						if (INFO1[i] == variant_INFO_array[v]) : # my current VT INFO is already saved? - yes, no need to save it again
							VT_already_stored=True# 
							break #
				if (not VT_already_stored) : # if not already saved - then save the INFO field name and values
						variant_INFO_array.append(INFO1[i])
						variant_INFO_val_array.append(INFO1[i+1])
	if (len(variant_INFO_array) < 1) : # NO VT INFO fields found
#		print ("no VT INFO field to match the VT PPF lines") #
		val=0	# then set a 0 score,
		return val # no further processing
######
	logging.debug('VT WORK array : %s / %s', str(variant_INFO_array), str(variant_INFO_val_array))
#	print ("VT WORK array :",variant_INFO_array, variant_INFO_val_array,"\n") #
#####
# now work these INFO fields
	for v in range(len(variant_INFO_array)): # move thru the saved VT INFO fields
#		print ("VT array :", eyeVarP_conf.VT_array[j]) # move to pred_array slot
#		print ("content-",variant_INFO_array[v])				#
#		print ("INFO_details : ",variant_INFO_val_array[v]) #
		predval=re.split('\||&',variant_INFO_val_array[v]) # split into individual values using delimiter | and & 
#
		for k in range(len(predval)): # for the variant info fields
			VT_k_has_annotation=False # initalise for this set of VT annotaion
			VT_k_found=False # initalise for this set of VT annotaion
#
			for j in range(len(eyeVarP_conf.VT_array)): # move thru the VT lines in PPF
				if (variant_INFO_array[v] == eyeVarP_conf.VT_array[j][1]) : # the variant functional annotation we want? 
#					print ("compare VT : vcf annotation - ",predval[k],"/ PPF VT- ",eyeVarP_conf.VT_array[j])  #
					if (predval[k] != ".") : # an annotated value
						VT_k_has_annotation=True #
						if ( predval[k] == eyeVarP_conf.VT_array[j][2]) :	 # is this the specific variant type we are looking for   
							VT_k_found=True # found a match
#							print ("vtfound") #
							wrkval=int(eyeVarP_conf.VT_array[j][3])	 # yes - return the score for the variant
							if (wrkval > val ) : #check for higher predictor score
#								print ("new variant's VT score : ",wrkval) #
								val=wrkval	
							break # stop loop as this one has found a match
					else : # there is no annotation for this INFO field
						break
				# end of for loop
				#
#			print ("setting the inds - found -", VT_k_found, "anno -", VT_k_has_annotation,"\n") # This annotation was found in PPF and scored
			if (VT_k_found and VT_k_has_annotation) : # This annotation was found in PPF and scored
				VT_found_for_variant=True
			elif (VT_k_has_annotation) : # only had a valid annotation but not in PPF
				VT_has_annotation_value_in_variant=True
#
# have finished checking all the VT lines in the PPF against variant
# now determine what is the prority score to pass back for the variant.
#
#	print ("variant's before correction VT score : ",val) # this is the score based on matches to the PPF VT lines
#
	logging.debug('final inds - VT_found: %s / ANNO exit in variant - %s', VT_found_for_variant, VT_has_annotation_value_in_variant)
#	print ("final inds - VT_found -", VT_found_for_variant, "ANNO exit in variant -", VT_has_annotation_value_in_variant,"\n") # This annotation was found in PPF and scored
	if (not VT_found_for_variant and VT_has_annotation_value_in_variant and (val == -999 )) : # variant's annotations were not in any of the VT PPF 
#		print ("no VT annotations matched in the PPF VT lines for the variant") #
		val=0	# then set a 0 score,
	elif (not VT_found_for_variant and not VT_has_annotation_value_in_variant ) : # no annotation found in the variant 
#		print ("no VT annotations exist in the variant") #
		val=0	# then set a 0 score,
	elif (VT_found_for_variant and (val < 0 )) : # there was a match but only for low score
#		print ("there was a match but only for low score") #
		if (VT_has_annotation_value_in_variant) : # if there were other annotated values in the variant
#			print ("but there were other VT annotation for variants not in PPF") #
			val=0	# then set a 0 score,
#			return val	# then OK
#		else : # no other annotated field, so return the low value	
#			print ("so return matched low value VT fields") #
#			return val	# then OK
#	elif (VT_found_for_variant and (val > 0 )) : # there was a match VT value from the PPF to the variant, and with a positive score 
#		print ("matched VT fields") #
#		return val	# then OK
#	print ("variant's final VT score : ",val) #
##
	return val #
##
#
