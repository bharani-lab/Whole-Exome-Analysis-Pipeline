###########################################################################################################
# #
###########################################################################################################
#
import sys, re, glob, os, subprocess, time #
import numpy as np #
import eyeVarP_conf, eyeVarP_6_1_convertVEP #
from shutil import copyfile #
#
tab='\t' # 
nl='\n' #
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
info_msg1_1="eyeVarP utility : Invalid number of inputs." 
#
Ast_ln=eyeVarP_conf.nl+"*********************************************************************************************"+eyeVarP_conf.nl #
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
#	print (supplied_args) #
#
	if (supplied_args < 5 ):  # arg [0] is the python program, must have at least 5 args, 
		print (info_msg1_1) #
		return 1 #
	else :
		eyeVarP_conf.ufunction=sys.argv[2] #
		eyeVarP_conf.input_file=sys.argv[3] # 
		eyeVarP_conf.final_output_file=sys.argv[4] # 
#		print (eyeVarP_conf.input_file)
	#
#		eyeVarP_conf.final_output_file=eyeVarP_conf.output_dir+"variant_merge_file_"+suffix+".txt" #
#		eyeVarP_conf.temp_output_file=eyeVarP_conf.output_dir+"variant_merge_file_"+suffix+"_tmp.txt" #
#		eyeVarP_conf.working_file1=eyeVarP_conf.output_dir+"working_file1_"+suffix+"_tmp.txt" #
#		eyeVarP_conf.working_file2=eyeVarP_conf.output_dir+"working_file2_"+suffix+"_tmp.txt" #
#		eyeVarP_conf.full_file1=eyeVarP_conf.output_dir+"full_file1_"+suffix+"_tmp.txt" #
#		eyeVarP_conf.full_file2=eyeVarP_conf.output_dir+"full_file2_"+suffix+"_tmp.txt" #
#		eyeVarP_conf.sort_file1=eyeVarP_conf.output_dir+"sort_file1_"+suffix+"_tmp.txt" #
#		eyeVarP_conf.sort_file2=eyeVarP_conf.output_dir+"sort_file2_"+suffix+"_tmp.txt" #
		print ("output : ",eyeVarP_conf.final_output_file) #
	
	return 0 #
#
#
###########################################################################################################
#
###########################################################################################################
def main(): #
##
	global scorefile #
	global gene_loc #
	global infile_shape #
#
	eyeVarP_conf.init() #
	#
	print ("eyeVarP utilities - Main") #
	#
	if (initial_setup() != 0): #
#		print "no good" #
		return #
#
	if (eyeVarP_conf.ufunction=="convertVEP"): #
#		print ("opt1") 
		print ("eyeVarP : utility option - ", eyeVarP_conf.ufunction,"- selected.")
		eyeVarP_6_1_convertVEP.main()  #
	else :
		print ("eyeVarP : utility option - ", eyeVarP_conf.ufunction,"- not found ")
#		
#
#
############################################################################################################
