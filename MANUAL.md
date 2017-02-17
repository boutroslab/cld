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

	    --version							- Show version.
	    --help								- Show this message.
	    
	    

In the following table every parameters as can be defined in the parameter-file is explained in more detail.

| parameter | explanations | type/value |
| ------------- | ------------- | ------------- |
| purpose_exclusive | defines if the pupose exclusive choice should affect the design filtering criteria | boolean (1 or 0) |
| min_length | defines the minimum length of the protospacer (5' sequence before the PAM) | numeric |
| max_length | defines the maximum length of the protospacer (5' sequence before the PAM) | numeric |
| min_G | defines minimum total G content | numeric |
| max_G | defines maximum total G content | numeric |
| min_A | defines minimum total A content | numeric |
| max_A | defines maximum total A content | numeric |
| min_C | defines minimum total C content | numeric |
| max_C | defines maximum total C content | numeric |
| min_T | defines minimum total T content | numeric |
| max_T | defines maximum total T content | numeric |
| right_homology | if homology arms are chosen, this number defines the length of the arm 5' of the target site | numeric |
| left_homology | if homology arms are chosen, this number defines the length of the arm 3' of the target site | numeric |
| downstream_window | defines the 5' nucleotide distance window of a design to the start or stop codon for tagging analysis | numeric |
| upstream_window | defines the 3' nucleotide distance window of a design to the start or stop codon for tagging analysis | numeric |
| number_of_CDS | defines the allowed nucleotid number within CDS downstream of the start codon, if knockout is chosen as the purpose criterium | numeric |
| minspacerlength | defines the minimum spacer length, if paired design is chosen as the purpose criterium (as for the double nickase or FokI Cas9 approach) | numeric |
| maxspacerlength | defines the maximum spacer length, if paired design is chosen as the purpose criterium (as for the double nickase or FokI Cas9 approach) | numeric |
| preceding | defines if the protospacer should begin with a specific base (example: U6 promotor would favour G at this position) | IUPAC coded Nucleotide |
| PAM_location | defines if the PAM motif is 3' or 5' located with respect to the protospacer | 3_prime or 5_prime |
| PAM | defines the PAM sequence, which can be NAG, NGG or any (allowing both) | IUPAC coded Nucleotides |
| ignore_intergenic | defines if off-targets which are not in any gene should be ignored | boolean (1 or 0) |
| purpose | defines following purposes: knockdown/-out as pupose requires designs to hit in coding sequences near the start codon of a gene, N-terminal tagging requires the start codon to be targeted and C-terminal tagging requires the stop codon to be targeted by sgRNAs | knockout, n-tagging, c-tagging, non-coding, CRISPRa or CRISPRi |
| gene_exclusive | defines if the sgRNA needs to target a region within the targeted gene (for CRISPRa/i 500 before and after the gene are parsed too) | boolean (1 or 0) |
| exon_exclusive | defines if the sgRNA needs to target a region within an exon | boolean (1 or 0) |
| CDS_only | defines if the sgRNA needs to target a region within a coding region | boolean (1 or 0) |
| CpG_exclusive | defines if the sgRNA is allowed to target a region within a CpG island | boolean (1 or 0) |
| specific_exon | defines if a specific exon number is to be targeted | numeric |
| retrieve_recomb_matrix | defines if the sequences for homology arms should be computed and reported | boolean (1 or 0) |
| bowtie_version | defines which version of bowtie or blast should be used for off-target analysis. Bowtie is more sensitive to mismatches of single designs, and bowtie2 is optimized for paired alignments of sequences. Here Blast tends to be the most sensitive towards less homologous sequences. For all mapping algorithms only full-length alignments are counted. | bowtie, bowtie2 or blast |
| offtargetdb | defines if off-targets should be searched in genomic sequences, sequences of annotated genes or exons of protein coding sequences | genomeDNA, gDNA or cDNA |
| targets-allowed | defines how many targets per design are tolerated before it is excluded from the report. This should be grater or equal than 1. Else the target is excluded as well and there are no designs passing.| numeric |
| unspecific_leading_bases | defines the number of 5' base pairs of the target site to be ignored for the off-target mapping | numeric |
| edit_distance_allowed | defines the edit distance (sum of all mismatch or INDEL positions) allowed during alignment to be still counted as off-target | numeric |
| bowtie_mode | define the bowtie mode as referenced in the bowtie2 manual | sensitive, very sensitive, fast, very-fast |
| sec_off_target | defines if sgRNA targets sites that are not in the genome of interest (for example: GFP etc.). Those  sequences need to be provided in an extra fasta formatted file in the database path and named 'secondary_off_targets.fasta'  | boolean (1 or 0) |
| max_per_exon | defines the maximum number of sgRNA allowed to be reported per exon | numeric |
| out_gff | defines if a gff should be generated  | boolean (1 or 0) |
| specific_transcript | defines if only a specific transcript (provided as an ENSEMBL TR ID) should be targeted. This is not applicable if more than one gene is searched. | ENSEMBL transcript ID or any |
| working_path | defines if an unix path to the results should be used, else results are created in the current working directory (.) | e.g. /data/workdir/ |
| databasepath | defines if an unix path to the folder containing CLD formatted databases should be generated | e.g. /data/databases/ |
| ref_organism | defines the reference organism as given in the name of the database and the sub-directories e.g. if the organisms is homo_sapiens, the database needs to have the prefix homo_sapiens | e.g. homo_sapiens , drosophila_melanogaster |
| data_type | defines if the input file contains official gene symbols, ENSEMBL IDs or genomic coordinates. Coordinates need to be given as ID, chromosome (Ensembl_type), start, end. Coordinate data need to be tab separated and different entries need to be newline separated | should be either ensemble_acc, gene_symbol or coordinates |
| ignore_missing_id | defines if the program should die if IDs are faced, that can not be found in the currently used database | boolean (1 or 0) |
| kind | defines if sgRNA target sites should be found in a single or paired mode (suitable for the paired nickase or FokI paired nuclease approach) | single or double |
| exclude_overlapping_genes | defines if sgRNA designs targeting multiple overlapping genes/ antisense transcripts should be excluded | boolean (1 or 0) |
| sort_by_rank | defines if sgRNAs should be ranked additionallly by an on-target score | boolean (1 or 0) |
| scores | defines the on-target score to be used. The preset scores are derived from the algorithms proposed by  Xu et al. 2015 and Doench et al. 2014. However they are only defined for a 20 nt protospacer adjacent to a NGG PAM. | xu_score, doench_old or custom |
| custom_score | defines a custom scoring function in perl code. the function needs to be unnamed and dependent on sequence information of the 30mer described in Doench et al.. Results of the function need to be numeric. | string in perl language defining a anonymous funtion, which acts on the Doench 30mer |
| cover_many_transcripts | defines if priority in sgRNA choice for the final library should be given on maximum coverage of all transcripts of a gene. All other scores will be ignored but shown in the resulting tables. | boolean (1 or 0) |

In the following table all output columns of the \*.tab files are explained in more detail:

|Column|Heading|Meaning|
|---|------|--------|
|0|Name|target site ID|
|1|Length|target length|
|2|Chromosome|target chromosome|
|3|Start|target start relative to the chromosome|
|4|End|target end relative to the chromosome|
|5|Strand|strand it will target|
|6|Nucleotide sequence|target site nucleotide composition of the form target_PAM|
|7|Gene Name|ID::GENE::STRAND|
|8|Transcripts|ENSEMBL transcript Ids overlapping with the target site "_"-separated|
|9|Transcript:: Exon|ENSEMBL transcript::exon Ids overlapping with the target site|
|10|Number of Cpg Islands hit|Number of Cpg Islands overlapping with the target site|
|11|Sequence around the cutside|if it is chosen to save a recombination matrix ist sequence is here|
|12|%A %C %T %G|nucleotide compositions in per cent|
|13|S-Score|specificity score|
|14|A-Score|annotation score|
|15|Custom-Score|by default the score from Doench et al. 2014 else the score deviated from the custom Perl scoring script|
|16|Doench-Score|Efficacy score as introduced by Doench et al. 2014 Nat. Biotech.|
|17|Xu-Score|Efficacy score as introduced by Xu et al. 2015 Gen.Res.|
|18|percent of total transcripts hit|per cent of transcripts of the targeted gene being hit by that putative sgRNA|
|19|Match-Target|target genes found by remapping the target site|
|20|Match-Chromosome|target chromosome found by remapping the target site|
|21|Match-Start|alignment start with respect to the estimate target chromosome|
|22|Match-End|alignment end with respect to the estimate target chromosome|
|23|Matchstring|alignment representation "M" for match "X" for mismatch "I" for insertion "D" for deletion|
|24|Editdistance|estimated edit distance of the alignment (X+I+D, sum of all deletions, insertions or exchanges neccessary to edit the matched sequence to the target sequence)|
|25|Number of Hits|estimated number of target sites in the respective genome with the off-target parameters specified|
|26|Direction|strandedness of the target alignment|
|27|Start_rti|target start with respect to the estimate target gene|
|28|End_rti|target end with respect to the estimate target gene|