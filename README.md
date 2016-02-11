![]( https://cdn.rawgit.com/boutroslab/cld/master/logo.png )

#CRISPR Library Designer (CLD): a software for the multispecies design of sgRNA libraries

ABSTRACT

Genetic screens using CRISPR/Cas9 are a powerful method for the functional analysis of genomes. Here we provide a fully integrated bioinformatics workflow for the design of custom single guide (sg) RNA libraries for a broad spectrum of organisms, termed CRISPR library designer (CLD). CLD can predict a high fraction of functional sgRNAs. An analysis on parameters that determine on-target efficiency of sgRNAs indicates that their prediction by CLD gives valuable insights into their efficiency in experiment. CLD enables the design of custom scalable, high-coverage sgRNA libraries for many species.

Quick-Start:

On Mac: If you have'nt, install Xquartz from http://www.xquartz.org/

When logging in remotely: log into your remote server by ssh -X

Download the GUI application according to your Operating system and unzip it.

Double click the application or open it by ./CLD in the Terminal.

Download the database for your organism of interest.

Enter its name in the field below.

Enter a gene list and go to the Parameter Tab to the End and Start your analysis.

Command-Line-Start:

cld can be called either with “--version”, printing its version number and copyrights, 
“--help” printing a more elusive help documentation and with “--task”. 

EXAMPLE to execute from the path containing all needed files:

cld -—task=end_to_end —-output-dir=. --parameter-file=./params.txt --gene-list=./gene_list.txt		    

cld can run 2 distinct tasks, database creation and 
library design.

Database creation is called using the “--task=make_database” command 
	giving the organism name of interest, as it is denoted in ENSEMBLs ftp folder structure
	e.g. homo_sapiens, and the rsync url to the current ftp server of ENSEMBL, examples 
 	can be found when cld  --help is called. After calling this function CRISPR will 
 	automatically download the latest toplevel FASTA, GFF and GTF files for the organism 
 	of interest and compile a database containing bowtie indexes, mygff files and 
 	reformatted sequence files. If not enough computing power is available to the user, 
 	these databases also might be downloaded from http://www.dkfz.de/signaling/crispr-downloads/. 

Library design can either be done in two steps: “cld 
	 --task=target_ident” and then “cld  --task=library_assembly” if the user wants 
 	to separate the two steps for example in order to only identify target sites without 
 	compiling a clonable library. 
 	Else “cld  --task=end_to_end” which automatically will perform the steps mentioned before 
 	after another and present the end-result in a user defined output folder. 
 	For reasons of manageability for high throughput design, output files are kept 
 	as simple and standardised as possible. However a genome wide library targeting 
 	the human genome quickly spans several GB depending on how strict the parameters 
 	are chosen. Since the end_to_end task takes most time we benchmarked its time 
 	consumption to be approximately 1 h wall-time for an 8-core cpu node.

 	
For running cld from the command line the syntax as outlined in the MANUAL must be used.

