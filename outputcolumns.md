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