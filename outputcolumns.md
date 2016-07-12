|Name|target site ID|
|----|--------------|
|Length|target length|
|Start|target start with respect to the input sequence (gene -500 if a gene name was the input)|
|End|target start with respect to the input sequence (gene -500 if a gene name was the input)|
|Strand|strand it will target|
|Nucleotide sequence|target site nucleotide composition of the form target_PAM|
|Gene Name|ID::GENE|
|Transcripts|ENSEMBL transcript Ids overlapping with the target site|
|Transcript:: Exon|ENSEMBL transcript::exon Ids overlapping with the target site|
|Number of Cpg Islands hit|Number of Cpg Islands overlapping with the target site|
|Sequence around the cutside|if it is chosen to save a recombination matrix ist sequence is here|
|%A %C %T %G|nucleotide compositions in per cent|
|S-Score|specificity score|
|A-Score|annotation score|
|Custom-Score|by default the score from Doench et al. 2014 else the score deviated from the custom Perl scoring script|
|percent of total transcripts hit|per cent of transcripts of the targeted gene being hit by that putative sgRNA|
|Target|target genes by remapping the target site|
|Match-start|alignment start with respect to the estimate target gene|
|Match-end|alignment end with respect to the estimate target gene|
|Matchstring|alignment representation "M" for match "X" for mismatch "I" for insertion "D" for deletion|
|Editdistance|estimated edit distance of the alignment (X+I+D)|
|Number of Hits|estimate number of target sites in the respective genome with the off-target parameters specified|
|Direction|strandedness of the target alignment|
|Doench-Score|Efficacy score as introduced by Doench et al. 2014 Nat. Biotech.|
|Xu-Score|Efficacy score as introduced by Xu et al. 2015 Gen.Res.|