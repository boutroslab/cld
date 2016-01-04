cld can be called either with “--version”, printing its version number and copyrights, 
“--help” printing a more elusive help documentation and with “--task”. 

EXAMPLE to execute from the path containing all needed files:

cld -—task=end\_to\_end —-output-dir=. --parameter-file=./params.txt --gene-list=./gene\_list.txt		    

cld can run 2 distinct tasks, database creation and 
library design.

Database creation 
	is called using the “--task=make\_database” command 
	giving the organism name of interest, as it is denoted in ENSEMBLs ftp folder structure
	e.g. homo\_sapiens, and the rsync url to the current ftp server of ENSEMBL, examples 
 	can be found when cld  --help is called. After calling this function CRISPR will 
 	automatically download the latest toplevel FASTA, GFF and GTF files for the organism 
 	of interest and compile a database containing bowtie indexes, mygff files and 
 	reformatted sequence files. If not enough computing power is available to the user, 
 	these databases also might be downloaded from http://www.dkfz.de/signaling/crispr-downloads/. 

Library design
	can either be done in two steps: “cld 
	 --task=target\_ident” and then “cld  --task=library\_assembly” if the user wants 
 	to separate the two steps for example in order to only identify target sites without 
 	compiling a clonable library. 
 	Else “cld  --task=end\_to\_end” which automatically will perform the steps mentioned before 
 	after another and present the end-result in a user defined output folder. 
 	For reasons of manageability for high throughput design, output files are kept 
 	as simple and standardised as possible. However a genome wide library targeting 
 	the human genome quickly spans several GB depending on how strict the parameters 
 	are chosen. Since the end\_to\_end task takes most time we benchmarked its time 
 	consumption to be approximately 1 h wall-time for an 8-core cpu node.
 	
In the following table every parameters as can be defined in the parameter-file is explained in more detail.

| parameter | explanations | type/value |
| ------------- | ------------- |
| purpose\_exclusive	|	determins if the pupose exclusive choice should effect the design filtering criteria	|	boolean (true or false) |
| min\_length	|	minimum length of the protospacer ( 5' sequence befor the PAM)	|	numeric |
| max\_length	|	minimum length of the protospacer ( 5' sequence befor the PAM)	|	numeric |
| min\_G	|	minimum total G content	|	numeric |
| max\_G	|	maximum total G content	|	numeric |
| min\_A	|	minimum total A content	|	numeric |
| max\_A	|	maximum total A content	|	numeric |
| min\_C	|	minimum total C content	|	numeric |
| max\_C	|	maximum total C content	|	numeric |
| min\_T	|	minimum total T content	|	numeric |
| max\_T	|	maximum total T content	|	numeric |
| right\_homology	|	if homology arms should be prosposed, this is the length of the arm 5' of the target site	|	numeric |
| left\_homology	|	if homology arms should be prosposed, this is the length of the arm 3' of the target site	|	numeric |
| downstream\_window	|	5' tolerance window of the proximity of a design to the start or stop codon for tagging analysis can be defined here 	|	numeric |
| upstream\_window	|	3' tolerance window of the proximity of a design to the start or stop codon for tagging analysis can be defined here 	|	numeric |
| number\_of\_CDS	|	if knockout is chosen as the purpose criterium, than the allowed number of CDS downstream of the start codon can be defined here	|	numeric |
| minspacerlength	|	if paired design as for the double nickase, or FokI approach is chosen than the minimum spacer length can be defined here	|	numeric |
| maxspacerlength	|	if paired design as for the double nickase, or FokI approach is chosen than the maximum spacer length can be defined here	|	numeric |
| preceding	|	here can be defined if the protospacer should begin with a specific base, as for example when the U6 promotor would favour Guanine at this position	|	A, C, G, T or any |
| PAM	|	defines the Protospacer motif, which here ccan be NAG, NGG or any (allowing both)	|	NAG, NGG or any |
| ignore\_intergenic	|	if off-targets, which are not in any gene should be ignored	|	true or false |
| purpose	|	Knockdown/-out as pupose requires designs to hit in coding sequences near the start codon of a gene , N-Terminal tagging requires that the start codon is targeted and C-Terminal tagging requires the stop codon to be span by the design	|	knockout, n-tagging, c-tagging, non-coding, CRISPRa or CRISPRi |
| gene\_exclusive	|	define if desired need to target a region within a gene	|	boolean (true or false) |
| exon\_exclusive	|	define if desired need to target a region within a exon/transcript/mRNA	|	boolean (true or false) |
| CDS\_only	|	define if desired need to target a region within a coding region	|	boolean (true or false) |
| CpG\_exclusive	|	define if desired are allowed to target only a region without a CpG island	|	boolean (true or false) |
| specific\_exon	|	define if a specific exon number is to be targeted	|	numeric |
| retrieve\_recomb\_matrix	|	if the sequences for homology arms should be computed and reported	|	boolean (true or false) |
| bowtie\_version	|	version of bowtie to be used for off-target analysis. Here bowtie is more sensitive to mismatches of single designs, and bowtie2 is optimized for paired alignments of sequences 	|	bowtie or bowtie2 |
| offtargetdb	|	define if offtargets should be searched in genomic sequences, c-DNA sequences or in gene models with introns but withoout intergenic regions	|	gDNA, cDNA |
| off-targets-allowed	|	define how many off-targets per design are tolerated before it is excluded from the report	|	numeric |
| unspecific\_leading\_bases	|	define the number of 5' base pairs of the target site to be ignored for the off-target mapping	|	numeric |
| edit\_distance\_allowed	|	define the edit distance (sum of all mismatch or INDEL positions) allowed in alignment to still count it as off-target	|	numeric |
| bowtie\_mode	|	define the bowtie mode as referenced in the bowtie2 manual	|	sensitive, very sensitive, fast, very-fast |
| sec\_off\_target	|	it can be checked if designs do target off-target site, which are not nessecarily in the genome of interest. Those secondary sequences need to be provided in an extra fasta formatted file in the database path and named 'secondary\_off\_targets.fasta' 	|	boolean (true or false) |
| max\_per\_exon	|	maximum number of sgRNA allowed to be reported per exon	|	numeric |
| out\_gff	|	should a gff ouput be generated	|	boolean (true or false) |
| specific\_transcript	|	if only a specific transcript ( given by an ENSEMBL TR ID) should be targetedthis makes no sense if more than one gene is searched	|	ENSEMBL transcript ID or any |
| match\_info	|	should a detailed alignment information be printed on any sgRNA mapping	|	boolean (true or false) |
| draw\_html\_report	|	should an html report be printed	|	boolean (true or false) |
| databasepath	|	a unix path to the cld formatted databases containing folder	|	e.g. /data/databases/ |
| ref\_organism	|	the reference organism as in the name of the database and the sub directories e.g. if the database has the prefix homo\_sapiens the organism should be named homo sapiens	|	e.g. homo\_sapiens , dmel
| data\_type	|	ensemble\_acc or fasta	|	define if the input file contains fasta formatted sequences or ENSEMBL IDs, only if ENSEMBL IDs are given the entire complexity of cld can be used |
| ignore\_missing\_id	|	should the program die if ids are faced, which can not be found int the currently used database	|	boolean (true or false) |
| kind	|	should sgRNA target sites be found in single or paired mode suitable for the paired nickase or FokI paired nuclease approach	|	single or double |
| exclude\_overlapping\_genes	|	should designs targeting multiple overlapping genes/ antisense transcripts be targeted at all	|	boolean (true or false) |


For running cld from the command line the following syntax must be used.

Usage: cld  --task=end\_to\_end [options=value] ...
Options:
	    --task=<task option>
		 make\_database 				to provide an cld ready data base.
		    --organism=<string>		Specify an organism to build the database for.
								    it must be one of the organisms available in ENSEMBLs ftp repository.
								    And in the same format as its ENSEMBL ftp directoy name.
								    E.g.: drosophila\_melanogaster or homo\_sapiens

		    --rsync-link=<rsync://path/to/dir>	         Specify an ftp repository to build the database from.
								    it must be one of the organisms available in ENSEMBLs ftp repository.
								    And in the same format as its ENSEMBL rsync directoy path.
								    E.g.: 
								    rsync://ftp.ensembl.org/ensembl/pub/release-81/
								    
								    rsync://ftp.ensemblgenomes.org/all/pub/protists/current/

									rsync://ftp.ensemblgenomes.org/all/pub/plants/current/

									rsync://ftp.ensemblgenomes.org/all/pub/fungi/current/
									
									rsync://ftp.ensemblgenomes.org/all/pub/metazoa/current/

		 target\_ident 					to identify target sequences.
		    --output-dir=<path/to/dir>			- a working directory as unix path to directory.
		    --parameter-file=<path/to/dir>		- a parameter file in cld format as path to file.
		    --gene-list=<path/to/dir>			- a gene list file with ENSEMBL IDs new-line seprated as path to file.
			--scoring-module=<path/to/dir>		- the path and filename of a file defining a perl scoring function
			--scoring-weights=<path/to/dir>		- the path and filename of a file defining a weight matrix with multipliers for the different scoring aspects

		 library\_assembly 				to format a library from an identification folder.
		    --output-dir=<path/to/dir>			- a working directory as unix path to directory.
		    --parameter-file=<path/to/dir>		- a parameter file in cld format as path to file.
		    --gene-list=<path/to/dir>			- a gene list file with ENSEMBL IDs new-line seprated as path to file. 
		    --cov=<int>						- Specify the minimum gene coverage as <int> default(15)
		    --lib-size=<int>				- Specify the maximum library size as <int> default(2000)
		    --lib-name=<string>				- Prefix for the final library as <string> default(test\_lib).
		    --5-prime=<string>				- Define the adapter to be put in 5' before the target site.
												default(CTGAGCTCATAGAAGACCTCACC)
		    --3-prime=<string>				- Define the adapter to be put in 3' behind the target site.
								   				default(GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG)
		    --cor-5-prime=<string>			- Specify if the first 5' baspair should be corrected to a G
								   				 may be "true" or "false" default :true.	    
		    --input-folder=<path/to/dir>		- Specify the input folder for library assembly.
								    			this folder must be prepared by --task= target\_ident
			--spread-over-transcripts=<string>	- should the designs be equally spread oer the different transcripts of the gene
													-an be : true or false (default:true)

		 end\_to\_end 							to perform and end-to-end analysis from target identification to library formatting
		    --output-dir=<path/to/dir>			- a working directory as unix path to directory.
		    --parameter-file=<path/to/dir>		- a parameter file in cld format as path to file.
		    --gene-list=<path/to/dir>			- a gene list file with ENSEMBL IDs new-line seprated as path to file. 
		    --cov=<int>						- Specify the minimum gene coverage as <int> default(15)
		    --lib-size=<int>				- Specify the maximum library size as <int> default(2000)
		    --lib-name=<string>				- Prefix for the final library as <string> default(test\_lib).
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
	    	    
