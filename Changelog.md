**Changes:**

| version | modifications |
| ------------- | ------------- |
|v1.1.0|+   added strandedness information for single exon genes |
| |+   fixed GUI rules for CRISPRa/i |
| |+   added option so TSS are only counted if belonging to the gene of interest |
| |+   changed bowtie off-target mode to -k ( offtargets allowed +1 ) so that bowtie always searches if designs have less or at least on off-target more than allowed. This way, sgRNAs are excluded which do have more than the allowed off-targets and included if less. Earlier one could allow more off-targets then bowtie was searching for. |
|v1.2.0|+   fixed overlapping criterion to beeing 3 bp downstream of the PAM (sgRNAend - 5bp) |
|v1.3.0|+   fixed error in version 1.2 where filter criteria like exons were not interpreted |
| |+   enabled automated deletion of single gene files in the server appliation so only the usmmary is kept|
| |+   fixed all binary builds |
|v1.3.1|+   fixed error in version 1.3.0 where gene names with underscores caused major trouble |
|v1.4.0|+   removed html output option |
| |+   fixed alignment match coordinates to be according to unspecific leading base pairs |
| |+   added the columns, Chromosome, Match-Chromosome, Start_rti, End_rti (relative to input) |
| |+   fixed Manual and Readme for the new cloumn headings |
| |+   fixed sgRNA ranking accordingly |
|v1.4.1|+   fixed a problem when purpose was set to non-coding so that not all possible sgRNA were considered for the library fitering |
|v1.4.3|+   Added TTTV as PAM |
| |+   fixed output fasta file when PAM is 5' |
|v1.4.4|+   fixed 5_prime PAM matchstring |
| |+   fixed 5_prime PAM off-target behaviour |
| |+   changed off-targets-allowed parameter to targets-allowed the number has to be greater than 1. Else all sgRNAs are excluded also those that only have one \"on-target\" |
|v1.4.5|+   fixed error handling with doubled gene ids. Now a warning is thrown and genes are processed in alphanumerical order |
| |+   fixed issues when Genes had similar names with numerics e.g. RPL3 and RPL30 |
| |+   fixed counting of missing genes|
|v1.4.6|+   fixed an issue with library output in fasta format, where only the PAM was in the output oligo |
