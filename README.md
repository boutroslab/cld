![]( https://cdn.rawgit.com/boutroslab/cld/master/logo.png )

[![DOI](https://zenodo.org/badge/20669/fheigwer/cld.svg)](https://zenodo.org/badge/latestdoi/20669/fheigwer/cld)

#CRISPR Library Designer (CLD): a software for the multispecies design of sgRNA libraries

CITATION

[F. Heigwer\*, T. Zhan\*, M. Breinig, J. Winter, D. Brügemann, S. Leible, M. Boutros, CRISPR library designer (CLD): software for multispecies design of single guide RNA libraries, Genome Biol., 2016, DOI:10.1186/s13059-016-0915-2](http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0915-2 "Access manuscript directly")

ABSTRACT

Here we describe CRISPR library designer (CLD), an integrated bioinformatics application for the design of custom single guide RNA (sgRNA) libraries for all organisms with annotated genomes. CLD is suitable for the design of libraries using modified CRISPR enzymes and targeting non-coding regions. To demonstrate its utility, we perform a pooled screen for modulators of the TNF-related apoptosis inducing ligand (TRAIL) pathway using a custom library of 12,471 sgRNAs.

**Quick-Start:**

On Mac/Linux\:
 1. If you haven't, install Xquartz from http://www.xquartz.org/
   
   1.1 When logging in remotely: log into your remote server by ssh -X
 2. Download CLD_GUI_Mac.zip or CLD_GUI_Ubuntu.zip according to your Operating system and unzip it.
 3. Double click on the application or open it by ./CLD in the Terminal.
 4. Download the database for your organism of interest.
 5. Enter its name in the reference organism field on the start page.
 6. Enter a list of gene identifiers in the "Gene List" tab and go to the "Design Parameter" tab to set your parameters.
 7. Go to the "Start Analysis" tab to start sgRNA search.
 8. The results will be created in the selected output directory ("Input/Output" tab).

**Command-Line-Start:**

cld can be called either with “--version”, printing its version number and copyrights, 
“--help” printing a more elusive help documentation and with “--task”. 

EXAMPLE to execute from the path containing all needed files:

cld --task=end_to_end --output-dir=. --parameter-file=./params.txt --gene-list=./gene_list.txt		    

cld can run 2 distinct tasks, database creation and 
library design.

Database creation is called using the “--task=make_database” command 
	giving the organism name of interest, as it is denoted in ENSEMBLs ftp folder structure
	e.g. homo_sapiens, and the rsync url to the current ftp server of ENSEMBL, examples 
 	can be found when cld  --help is called. After calling this function CLD will 
 	automatically download the latest toplevel FASTA, GFF and GTF files for the organism 
 	of interest and compile a database containing bowtie indexes, mygff files and 
 	reformatted sequence files. If not enough computing power is available to the user, 
 	these databases also might be downloaded from http://www.dkfz.de/signaling/crispr-downloads/. 

Library design can either be done in two steps: “cld 
	 --task=target_ident” and then “cld  --task=library_assembly” if the user wants 
 	to separate the two steps for example in order to only identify target sites without 
 	compiling a clonable library. 
 	Else “cld  --task=end_to_end” which automatically will perform the steps mentioned before 
 	and present the end-result in a user defined output folder. 
 	For reasons of manageability for high throughput design, output files are kept 
 	as simple and standardised as possible. However a genome wide library targeting 
 	the human genome quickly spans several GB depending on how strict the parameters 
 	are chosen. Since the end_to_end task takes most time we benchmarked its time 
 	consumption to be approximately 1 h wall-time for an 8-core cpu node.

 	
For running cld from the command line the syntax as outlined in the MANUAL must be used.

