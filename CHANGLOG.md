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
