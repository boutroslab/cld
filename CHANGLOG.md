**Changes:**

v   1.1.0   +   added strandedness information for single exon genes
            +   fixed GUI rules for CRISPRa/i
            +   added option so TSS are only counted if belonging to the gene of interest
            +   changed bowtie off-target mode to -k ( offtargets allowed +1 ) so that bowtie
                always searches if designs have less or at least on off-target more than allowed
                this way sgRNAs are excluded which do have more than the allowed off-targets and included if less
                earlier one could allow more off-targets then bowtie was searching for