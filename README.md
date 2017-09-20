# wgs_tool

Introduction

wfgs is a tool for estimation of within family genomic breeding values. It uses gmatrix (Guosheng Su and Per Madsen) program for setting up G-matrix and WOMBAT (Meyer, 2007) for fitting mixed models. 

Usage of wfgs tool 

wfgs is expected to run form a command line in Linux environment by typing the appropriate command and hitting return to start the program running. The general form of the command is:

  ./wfgs	parameter_file nfam
  
where parameter_file is the name of the ‘parameter file’ and the nfam is the number of families in the dataset. If there is only one family in the data, nfam can be omitted and if nfam is not specified the program assumes only single family information is provided in the data sets.

Input files 

The minimal requirement of wfgs are parameter file, data file and genotype file. 

Parameter file: this file contains all information that specifies everything wfgs needs to know 	about the input files and their outline and user running preferences (see parameter 	file 	section).

Data file: this file should at least contains ID of individuals, family ID, Sire ID, dam ID and at 	least one trait’s phenotype. The file should be without headers and name of the column 	should be provided in the parameter file.

	The 1st column should contain the individual ID
	The 2nd column should contain the family ID
	The 3nd column should contain the sire ID
	The 4th column should contain the dam ID
	The order of the remaining columns can be arbitrary 
  
Genotype file: the marker file should have one line for one individual 
	First column is ID number
	The 2nd and 3rd columns are alleles of locus 1
	The 4th and 5th columns are alleles of locus 2
	So on …
  
Parameter file

Wfgs uses a parameter file in a format of R language. The parameter file is read line by line and any line beginning with a ‘#’ is considered as comment line and therefore skipped. Wgs relies on specific keywords (codes) at the beginning of each line to distinguish between different types of information given. The ordering of the lines in the ‘parameter file’ is arbitrary. The format of a line in the ‘parameter file’ looks like:

KEYWORD = c (parameter1, parameter2 …)
The keywords (codes) are: 
PHENO = c(“name of the phenotype file”, “ID”, “FamID”, “Sire”, “Dam”, …)
Specifies the names of the phenotype file and describes its content. The first content tells the name of the data file and the next four columns should be constant as stated above. The remaining columns are arbitrary and no limit in the number of columns to include in the data file. 
GENO = “Genotype file name”
This code specifies the name of the genotype file. 
COM = “User define comment about the analysis”
Users can choose to give a comment to an analysis. This is optional.  
NTRAIT = value
This code specifies the number of traits for the analysis 
VAR = c (value1, value2, value3, … valueN)
This code specifies the covariance components for the residual and all the random effects. The order of the covariance components is first for the residual effect and then follows the order of the random effects as given in the RAN code. For q traits, the q(q + 1)/2 elements of the upper triangle of the covariance matrix, given row-wise (i.e. σ_1^2, σ_12, …, σ_1q, σ_2^2, σ_23, …,σ_q^2). This should be repeated for all random effects.
RAN = c (value1, value2, value3, … valueN)
This code specifies the location of random effects in the data file. 
REL = c (“GIN”/ “NRM”/ “IDE”)
This code describes the covariance structure of the random effects. Valid codes for the covariance structures are: 
 GIN: denotes that the random effect is distributed proportional to a relationship or correlation matrix. The user must supply the inverse of this matrix in the form outlined …
IDE: denotes that the different levels of the random effect are uncorrelated. 
NRM: denotes that the random effect is distributed proportionally to the numerator relationship matrix. If this code is given, a pedigree file must be provided.
FIX = c (value1, value2, value3, … valueN)
This code specifies the location of fixed effects in the data file. 
COV = c (value1, value2, value3, … valueN)
This code specifies the location of covariable in the data file. 
TRLOC = c (value1, value2, value3, … valueN)
This code specifies the location (columns) of phenotypes (traits) in the data file. 
MAF = value
This code specifies the minimum allele frequency. Markers with MAF lower than specified value will be removed from the analysis. If not specified, the default values of 0.05 is used.    
FREQMETHOD = 1 or 2 
This code specifies (fixes) the method to get marker allele frequency. 
1: calculate from the data
2: M-matrix divide by sqrt(2pq) for each locus
DIAG = value
This code specifies a value to be added to the diagonal of the G-matrix in order to male G-matrix positive definitive. The value usually varies form 0 – 0.02. If not specified, a default value of 0.01 is added to the diagonal.
GMATOUT = 1 or 2
This code specifies the user preference to whether to print out the G-matrix or not. If not specified, the default is not to print out. 
0: do not print out G-matrix
1: print out G-matrix
IGMATOUT = 1 or 2
This code specifies the user preference to whether to print out inverse of the G-matrix or not. If not specified, the default is to print out. 
0: do not print out G-matrix
1: print out G-matrix
SCALEMETHOD  = 1 or 2
Method to scale G-matrix
1: divide by sum of 2pq.
2: M-matrix divide by sqrt(2pq) for each locus

Output 

At the end of each run, wgs generates a folder with the name of the parameter file st two outputs, but produce more file depending on the running options. The output files are:

SumModel.out: This file gives a summary about the model of analysis specified and the 	corresponding features of the data found. Statistics given include means, standard 	deviations and ranges for traits and covariables, and numbers of levels found for the 	effects in the model of analysis.

FixSolutions.out: This file lists the generalized least-squares solutions for all fixed effects fitted.

RnSoln_random effect name.dat: Solutions for each random effect are written to a separate file. 	These files have names RnSoln_random effect name.dat, with random effect name representing the name of the random effect. Columns in these files are:

Gmat.log: Log information from ‘gmatrix’ run
wombat.log: Log information from WOMBAT run

 
References 

Meyer, K. 2007. WOMBAT – A tool for mixed model analyses in quantitative genetics by 	REML. J. Zhejiang Uni. SCIENCE B 8:815-821.
Guosheng Su and Per Madsen (http://webcache.googleusercontent.com/search?q=cache:7Y9Qo9fPo5oJ:dmu.agrsci.dk/Gmatrix/Doc/Previous/Gmatrix-User-Guide.pdf+&cd=2&hl=en&ct=clnk&gl=no)
