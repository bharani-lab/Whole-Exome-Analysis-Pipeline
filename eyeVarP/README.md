# eyeVarP
eyeVarP - Variant Prioritisation Ordering Tool. 

eyeVarP is a Python tool written to allow prioritisation of variants in ANNOVAR annotated VCF files. eyeVarP provides various functions for the purpose of speeding up variant discovery.
* priority  - priority tool
* genef     - gene filter
* samplef   - samples and inheritance model filtering
* stats     - general statistics on the eyeVarP priority output file
* merge     - merge multiple eyeVarP priority output files into a single eyeVarP priority output file
* utility     - eyeVarP utilities 

### Requirements
* Python 3.6.+ and NumPy
* Linux environment, or access to linux via ssh

### Installation
* Navigate to desired install directory and clone this repository.

 `git clone https://github.com/VCCRI/eyeVarP.git eyeVarP`

* Ensure that requirements are met. 
* Test that eyeVarP is working using the test_data and README.MD/tutorial provided.
* all done!

### Usage - eyeVarP.py  < option >

* **help**  -   will return a help screen                                                                                 
                                                                           
 * **priority**  - priority tool
    this option performs the variant proritisation process on the input samples VCF files. It will score each variant found for the supplied samples based on the weighting affixed to the annotations supplied in the Prioritisation Parameter File, PPF.
    
    Inputs :
    * location for output file+prefix
    * file of input VCF files and sample IDs
    * Prioritisation Parameter File (PPF)
    
    Output :
    * eyeVarP Priority Output List (VPOL)

* **genef**     - gene filter
   this option performs variant filtering of the VPOL based on genes supplied as input.
    
    Inputs :
    * location for output file+prefix
    * eyeVarP Priority Output List (VPOL)
    * gene list
    
    Output :
    * gene filtered eyeVarP Priority Output List (VPOL)     

* **samplef**   - samples and inheritance model filtering
   this option performs variant filtering of the VPOL based on a suplied ped format file. It can be used for a simple case-control filtering or an inheritance model filtering for a family trio.
    
    Inputs :
    * location for output file+prefix
    * eyeVarP Priority Output List (VPOL)
    * sample selection file (ped format)
    * proband sample ID (for inheritance model filtering)
    * inheritance model (DN/AD/AR/CH for inheritance model filtering)
    
    Output :
    * sample filtered eyeVarP Priority Output List (VPOL) 

* **stats**     - general statistics on the eyeVarP priority output file
       this option returns a summary statistic file for the VPOL supplied. It provide a small report listing the number of variants (the total number of variants, the number of scoring variants, the number of non-scoring variants), the number of genes and the number of samples. There is a breakdown of the number of variants found in each score 10% percentile range. The top 20 variants are also listed. A breakdown for each sample is also provided, with a table containing the number of variants in genes found above the percentile value supplied.
    
    Inputs :
    * location for output file+prefix
    * eyeVarP Priority Output List (VPOL)
    * a percentile value [1-99] for quick summary statistic
    
    Output :
    * statistic eyeVarP Priority Output List (VPOL) 

* **merge**     - merge multiple eyeVarP priority output list files into a single eyeVarP priority output list file
       this option provide the ability to merge various VPOLs into one single VPOL. This function allows large cohorts to be split into small groups to speed up proritisation processing and then output to be re-consolidated back into one single large cohort VPOL for downstream analysis or filtering.     
       
    Inputs :
    * location for output file+prefix
    * file containing a list of eyeVarP Priority Output List (VPOL)
        
    Output :
    * merged eyeVarP Priority Output List (VPOL) 
    
 * **utilities**     - eyeVarP utilities.     
       
		    Utility 1 : convertVEP - this utility convert VEP annotated VCF into the standard VCF format.
				    Inputs :
					* full pathname of input VEP annotated VCF file
				    * full pathname of output VCF file, including directory path 
        
				    Output :
					 * standard VCF 
                                                                                                                                  
## see README.MD in the test_data directory for more details on each function.

## 
