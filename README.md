#CRISPR Library Designer (CLD): a software for the multispecies design of sgRNA libraries

ABSTRACT

Genetic screens using CRISPR/Cas9 are a powerful method for the functional analysis of genomes. Here we provide a fully integrated bioinformatics workflow for the design of custom single guide (sg) RNA libraries for a broad spectrum of organisms, termed CRISPR library designer (CLD). CLD can predict a high fraction of functional sgRNAs. An analysis on parameters that determine on-target efficiency of sgRNAs indicates that their prediction by CLD gives valuable insights into their efficiency in experiment. CLD enables the design of custom scalable, high-coverage sgRNA libraries for many species.


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
 	
For running cld from the command line the following syntax must be used.

Usage: cld  --task=end_to_end [options=value] ...

Options:

	    --task=<task option>
		 make_database 				to provide an cld ready data base.
		    --organism=<string>		Specify an organism to build the database for.
								    it must be one of the organisms available in ENSEMBLs ftp repository.
								    And in the same format as its ENSEMBL ftp directoy name.
								    E.g.: drosophila_melanogaster or homo_sapiens

		    --rsync-link=<rsync://path/to/dir>	         Specify an ftp repository to build the database from.
								    it must be one of the organisms available in ENSEMBLs ftp repository.
								    And in the same format as its ENSEMBL rsync directoy path.
								    E.g.: rsync://ftp.ensembl.org/ensembl/pub/release-81/

		 target_ident 					to identify target sequences.
		    --output-dir=<path/to/dir>			- a working directory as unix path to directory.
		    --parameter-file=<path/to/dir>		- a parameter file in cld format as path to file.
		    --gene-list=<path/to/dir>			- a gene list file with ENSEMBL IDs new-line seprated as path to file.
			--scoring-module=<path/to/dir>		- the path and filename of a file defining a perl scoring function
			--scoring-weights=<path/to/dir>		- the path and filename of a file defining a weight matrix with multipliers for the different scoring aspects

		 library_assembly 				to format a library from an identification folder.
		    --output-dir=<path/to/dir>			- a working directory as unix path to directory.
		    --parameter-file=<path/to/dir>		- a parameter file in cld format as path to file.
		    --gene-list=<path/to/dir>			- a gene list file with ENSEMBL IDs new-line seprated as path to file. 
		    --cov=<int>						- Specify the minimum gene coverage as <int> default(15)
		    --lib-size=<int>				- Specify the maximum library size as <int> default(2000)
		    --lib-name=<string>				- Prefix for the final library as <string> default(test_lib).
		    --5-prime=<string>				- Define the adapter to be put in 5' before the target site.
												default(CTGAGCTCATAGAAGACCTCACC)
		    --3-prime=<string>				- Define the adapter to be put in 3' behind the target site.
								   				default(GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG)
		    --cor-5-prime=<string>			- Specify if the first 5' baspair should be corrected to a G
								   				 may be "true" or "false" default :true.	    
		    --input-folder=<path/to/dir>		- Specify the input folder for library assembly.
								    			this folder must be prepared by --task= target_ident
			--spread-over-transcripts=<string>	- should the designs be equally spread oer the different transcripts of the gene
													-an be : true or false (default:true)

		 end_to_end 							to perform and end-to-end analysis from target identification to library formatting
		    --output-dir=<path/to/dir>			- a working directory as unix path to directory.
		    --parameter-file=<path/to/dir>		- a parameter file in cld format as path to file.
		    --gene-list=<path/to/dir>			- a gene list file with ENSEMBL IDs new-line seprated as path to file. 
		    --cov=<int>						- Specify the minimum gene coverage as <int> default(15)
		    --lib-size=<int>				- Specify the maximum library size as <int> default(2000)
		    --lib-name=<string>				- Prefix for the final library as <string> default(test_lib).
		    --5-prime=<string>				- Define the adapter to be put in 5' before the target site.
								    			default(CTGAGCTCATAGAAGACCTCACC)
		    --3-prime=<string>				- Define the adapter to be put in 3' behind the target site.
								    			default(GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG)
		    --cor-5-prime=<string>			- Specify if the first 5' baspair should be corrected to a G.
			--spread-over-transcripts=<string>	- should the designs be equally spread oer the different transcripts of the gene
													-can be : true or false (default:true)
			--scoring-module=<path/to/dir>		- the path and filename of a file defining a perl scoring function
			--scoring-weights=<path/to/dir>		- the path and filename of a file defining a weight matrix with multipliers for the different scoring aspects

	    --version							- Show version.
	    --help								- Show this message.
	    
	    
