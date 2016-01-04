cld can be called either with “--version”, printing its version number and copyrights, 
“--help” printing a more elusive help documentation and with “--task”. 

EXAMPLE to execute from the path containing all needed files:

cld -—task=end_to_end —-output-dir=. --parameter-file=./params.txt --gene-list=./gene_list.txt		    

cld can run 2 distinct tasks, database creation and 
library design.

Database creation 
	is called using the “--task=make_database” command 
	giving the organism name of interest, as it is denoted in ENSEMBLs ftp folder structure
	e.g. homo_sapiens, and the rsync url to the current ftp server of ENSEMBL, examples 
 	can be found when cld  --help is called. After calling this function CRISPR will 
 	automatically download the latest toplevel FASTA, GFF and GTF files for the organism 
 	of interest and compile a database containing bowtie indexes, mygff files and 
 	reformatted sequence files. If not enough computing power is available to the user, 
 	these databases also might be downloaded from http://www.dkfz.de/signaling/crispr-downloads/. 

Library design
	can either be done in two steps: “cld 
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
 	
In the following table every parameters as can be defined in the parameter-file is explained in more detail.

| parameter | explanations | type/value |

| ------------- | ------------- |

| purpose_exclusive	|	determins if the pupose exclusive choice should effect the design filtering criteria	|	boolean (true or false) |

| min_length	|	minimum length of the protospacer ( 5' sequence befor the PAM)	|	numeric |

| max_length	|	minimum length of the protospacer ( 5' sequence befor the PAM)	|	numeric |

| min_G	|	minimum total G content	|	numeric |

| max_G	|	maximum total G content	|	numeric |

| min_A	|	minimum total A content	|	numeric |
| max_A	|	maximum total A content	|	numeric |
| min_C	|	minimum total C content	|	numeric |
| max_C	|	maximum total C content	|	numeric |
| min_T	|	minimum total T content	|	numeric |
| max_T	|	maximum total T content	|	numeric |
| right_homology	|	if homology arms should be prosposed, this is the length of the arm 5' of the target site	|	numeric |
| left_homology	|	if homology arms should be prosposed, this is the length of the arm 3' of the target site	|	numeric |
| downstream_window	|	5' tolerance window of the proximity of a design to the start or stop codon for tagging analysis can be defined here 	|	numeric |
| upstream_window	|	3' tolerance window of the proximity of a design to the start or stop codon for tagging analysis can be defined here 	|	numeric |
| number_of_CDS	|	if knockout is chosen as the purpose criterium, than the allowed number of CDS downstream of the start codon can be defined here	|	numeric |
| minspacerlength	|	if paired design as for the double nickase, or FokI approach is chosen than the minimum spacer length can be defined here	|	numeric |
| maxspacerlength	|	if paired design as for the double nickase, or FokI approach is chosen than the maximum spacer length can be defined here	|	numeric |
| preceding	|	here can be defined if the protospacer should begin with a specific base, as for example when the U6 promotor would favour Guanine at this position	|	A, C, G, T or any |
| PAM	|	defines the Protospacer motif, which here ccan be NAG, NGG or any (allowing both)	|	NAG, NGG or any |
| ignore_intergenic	|	if off-targets, which are not in any gene should be ignored	|	true or false |
| purpose	|	Knockdown/-out as pupose requires designs to hit in coding sequences near the start codon of a gene , N-Terminal tagging requires that the start codon is targeted and C-Terminal tagging requires the stop codon to be span by the design	|	knockout, n-tagging, c-tagging, non-coding, CRISPRa or CRISPRi |
| gene_exclusive	|	define if desired need to target a region within a gene	|	boolean (true or false) |
| exon_exclusive	|	define if desired need to target a region within a exon/transcript/mRNA	|	boolean (true or false) |
| CDS_only	|	define if desired need to target a region within a coding region	|	boolean (true or false) |
| CpG_exclusive	|	define if desired are allowed to target only a region without a CpG island	|	boolean (true or false) |
| specific_exon	|	define if a specific exon number is to be targeted	|	numeric |
| retrieve_recomb_matrix	|	if the sequences for homology arms should be computed and reported	|	boolean (true or false) |
| bowtie_version	|	version of bowtie to be used for off-target analysis. Here bowtie is more sensitive to mismatches of single designs, and bowtie2 is optimized for paired alignments of sequences 	|	bowtie or bowtie2 |
| offtargetdb	|	define if offtargets should be searched in genomic sequences, c-DNA sequences or in gene models with introns but withoout intergenic regions	|	gDNA, cDNA |
| off-targets-allowed	|	define how many off-targets per design are tolerated before it is excluded from the report	|	numeric |
| unspecific_leading_bases	|	define the number of 5' base pairs of the target site to be ignored for the off-target mapping	|	numeric |
| edit_distance_allowed	|	define the edit distance (sum of all mismatch or INDEL positions) allowed in alignment to still count it as off-target	|	numeric |
| bowtie_mode	|	define the bowtie mode as referenced in the bowtie2 manual	|	sensitive, very sensitive, fast, very-fast |
| sec_off_target	|	it can be checked if designs do target off-target site, which are not nessecarily in the genome of interest. Those secondary sequences need to be provided in an extra fasta formatted file in the database path and named 'secondary_off_targets.fasta' 	|	boolean (true or false) |
| max_per_exon	|	maximum number of sgRNA allowed to be reported per exon	|	numeric |
| out_gff	|	should a gff ouput be generated	|	boolean (true or false) |
| specific_transcript	|	if only a specific transcript ( given by an ENSEMBL TR ID) should be targetedthis makes no sense if more than one gene is searched	|	ENSEMBL transcript ID or any |
| match_info	|	should a detailed alignment information be printed on any sgRNA mapping	|	boolean (true or false) |
| draw_html_report	|	should an html report be printed	|	boolean (true or false) |
| databasepath	|	a unix path to the cld formatted databases containing folder	|	e.g. /data/databases/ |
| ref_organism	|	the reference organism as in the name of the database and the sub directories e.g. if the database has the prefix homo_sapiens the organism should be named homo sapiens	|	e.g. homo_sapiens , dmel
| data_type	|	ensemble_acc or fasta	|	define if the input file contains fasta formatted sequences or ENSEMBL IDs, only if ENSEMBL IDs are given the entire complexity of cld can be used |
| ignore_missing_id	|	should the program die if ids are faced, which can not be found int the currently used database	|	boolean (true or false) |
| kind	|	should sgRNA target sites be found in single or paired mode suitable for the paired nickase or FokI paired nuclease approach	|	single or double |
| exclude_overlapping_genes	|	should designs targeting multiple overlapping genes/ antisense transcripts be targeted at all	|	boolean (true or false) |


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
								    E.g.: 
								    rsync://ftp.ensembl.org/ensembl/pub/release-81/
								    
								    rsync://ftp.ensemblgenomes.org/all/pub/protists/current/

									rsync://ftp.ensemblgenomes.org/all/pub/plants/current/

									rsync://ftp.ensemblgenomes.org/all/pub/fungi/current/
									
									rsync://ftp.ensemblgenomes.org/all/pub/metazoa/current/

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
	    	    
