#!/usr/bin/perl
use strict;
#use warnings FATAL => 'all';
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Bio::SeqFeature::Generic; #important package to handle sequence formats and objects
use Bio::Location::Split; #library to make splitted location objects
use Set::IntervalTree; #library providing methods for using interval trees
use JSON::XS qw(encode_json decode_json); #library to encode perl objects in a file
use File::Slurp qw(read_file write_file); #library to write perl objects to a file
use List::MoreUtils qw{
      any all none notall true false
      firstidx first_index lastidx last_index
      insert_after insert_after_string
      apply indexes
      after after_incl before before_incl
      firstval first_value lastval last_value
      each_array each_arrayref
      pairwise natatime
      mesh zip uniq distinct minmax part}; #some math and list utilities, later needed to partition a list of entries
use List::Util qw(sum);
use Archive::Zip;
#use Bio::Graphics;
use Parallel::ForkManager; #important package to enable mutlithreading of the script
use Cwd;
use IPC::Cmd qw[can_run run];
use Getopt::Long qw(:config pass_through);	
use File::Grep qw( fgrep fmap fdo );
use Text::Wrap;
use Unix::Processors;

my $procs = new Unix::Processors;
my $max_parallel= my $parallel_number =$procs->max_online;
  
$| = 1;

my ($script_name,$script_version,$script_date,$script_years) = ('cld','0.1.7','2015-09-01','2013-2015');



my (
    $opt_task,
    $opt_param_file,
    $opt_gene_list,
    $opt_working_directory,
    $opt_library_name,
    $opt_coverage,
    $opt_total_lib_size,
    $opt_5_adapt,
    $opt_3_adapt,
    $correct_5_prime_G,
    $opt_input_folder,
    $opt_organism,
    $opt_rsync_link,
    $ver,
    $help,
	$cover_many_transcripts,
	$scoring_module,
	$scoreweight_file
    );
GetOptions(
	    'task=s'		=> \$opt_task,
	    'output-dir=s'	=> \$opt_working_directory,	
	    'parameter-file=s'	=> \$opt_param_file,
	    'gene-list=s'	=> \$opt_gene_list,	        
	    'cov=i'		=> 		\$opt_coverage,
	    'lib-size=i'	=> \$opt_total_lib_size,
	    'lib-name=s'	=> \$opt_library_name,
	    '5-prime=s'		=> \$opt_5_adapt,
	    '3-prime=s'		=> \$opt_3_adapt,
	    'cor-5-prime=s'	=> \$correct_5_prime_G,
	    'input-folder=s'	=> \$opt_input_folder,
	    'organism=s'	=> \$opt_organism,
	    'rsync-link=s'	=> \$opt_rsync_link,
	    'version'		=> \$ver,
	    'help'		=> \$help,
		'spread-over-transcripts=s'=> \$cover_many_transcripts,
		'scoring-module=s'=>\$scoring_module,
		'scoring-weights=s'=>\$scoreweight_file
	);
my $ver_str = "$script_name, version $script_version, $script_date\nAuthor $script_years Florian Heigwer\n";
my $help_str = qq{Usage: cld --task=end_to_end [options=value] ...
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
};

print (($ver ? $ver_str : ''), ($help ? $help_str : ''));
if(!can_run('bowtie') ){
	print "Cannot find bowtie. Please download and install from\n http://bowtie-bio.sourceforge.net/index.shtml\n";
	print $ver_str."\n".$help_str."\n";
	exit;
	}elsif(!can_run('bowtie2')){
	print "Cannot find bowtie2. Please download and install from\n http://bowtie-bio.sourceforge.net/bowtie2/index.shtml\n";
	print $ver_str."\n".$help_str."\n";
	exit;
}
for (my $i=0; $i<scalar(@ARGV); $i++)
{
    if (substr($ARGV[$i],0,1) eq '-' and $i < scalar(@ARGV)-1){
	if ($ARGV[$i] eq '-task'		){$opt_task		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-output-dir'		){$opt_working_directory= int($ARGV[++$i]); }	
	if ($ARGV[$i] eq '-parameter-file'	){$opt_param_file	= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-gene-list'		){$opt_gene_list	= int($ARGV[++$i]); }	        
	if ($ARGV[$i] eq '-cov'			){$opt_coverage		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-lib-size'		){$opt_total_lib_size	= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-lib-name'		){$opt_library_name	= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-5-prime'		){$opt_5_adapt		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-3-prime'		){$opt_3_adapt		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-cor-5-prime'		){$correct_5_prime_G	= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-organism'		){$opt_organism		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-rsync-link'		){$opt_rsync_link	= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-input-folder'	){$opt_input_folder		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-spread-over-transcripts'	){$cover_many_transcripts		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-scoring-module'	){$scoring_module		= int($ARGV[++$i]); }
	if ($ARGV[$i] eq '-scoring-weights'	){$scoreweight_file		= int($ARGV[++$i]); }
    }
}

if(!defined($opt_task)){
   print "\nThere must be any specific task specified.\n\n".$ver_str.$help_str;
   exit;
}



if($opt_task eq "make_database"){   
    make_database(
		  defined($opt_organism) ? $opt_organism : "drosophila_melanogaster",
		  defined($opt_rsync_link) ? $opt_rsync_link : "rsync://ftp.ensembl.org/ensembl/pub/release-77/"
		  );
}elsif($opt_task eq "target_ident"){
	if(!defined($opt_gene_list)){
		print "\nThere must be an gene list specified.\n\n".$ver_str.$help_str;
	   exit;
	}
    if(!defined($opt_gene_list) or !(-f $opt_gene_list)){ die "The gene list file $opt_gene_list could not be opened. Either the user has no rights the read it or the file does not exist." }
    if(!defined($opt_param_file) or !(-f $opt_param_file)){ die "The parameter list file $opt_param_file could not be opened. Either the user has no rights the read it or the file does not exist." }
    if(!defined($opt_working_directory) or !(-d $opt_working_directory)){ $opt_working_directory="~" }
    make_a_crispr_library(
			     $opt_param_file,
			     $opt_gene_list,
			     $opt_working_directory
			   ),
}elsif($opt_task eq "end_to_end"){
	if(!defined($opt_gene_list)){
		print "\nThere must be an gene list specified.\n\n".$ver_str.$help_str;
	   exit;
	}
    if(!defined($opt_gene_list) or !(-f $opt_gene_list)){ die "The gene list file $opt_gene_list could not be opened. Either the user has no rights the read it or the file does not exist." }
    if(!defined($opt_param_file) or !(-f $opt_param_file)){ die "The parameter list file $opt_param_file could not be opened. Either the user has no rights the read it or the file does not exist." }
    if(!defined($opt_working_directory) or !(-d $opt_working_directory)){ $opt_working_directory="~" }
    
    if(defined($opt_coverage)		and	$opt_coverage		=~m/([^\d]+)/g){ die "The gene coverage must be specified as an integer. Your input is not integer as it contains: $1." }
    if(defined($opt_total_lib_size)	and	$opt_total_lib_size	=~m/([^\d]+)/g){ die "The total library size must be specified as an integer. Your input is not integer as it contains: $1." }
    if(defined($opt_5_adapt)		and	$opt_5_adapt		=~m/([^ACGT]+)/g){ die "The 5' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
    if(defined($opt_3_adapt)		and	$opt_3_adapt		=~m/([^ACGT]+)/g){ die "The 3' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
    if(defined($correct_5_prime_G)	and	$correct_5_prime_G ne "true" and $correct_5_prime_G ne "false" ){ die "Correcting the 5' basepair to a G must be true or false" }    
    
    filter_library(		
		    make_a_crispr_library(
					    $opt_param_file,#"params.txt",
					    $opt_gene_list,# "still_left.tab"
					    $opt_working_directory #~
					),
		    defined($opt_library_name) ? $opt_library_name : "test_lib",
		    defined($opt_coverage) ? $opt_coverage : 15,
		    defined($opt_total_lib_size) ? $opt_total_lib_size : 2000, 
		    defined($opt_5_adapt) ? $opt_5_adapt : "CTGAGCTCATAGAAGACCTCACC", 
		    defined($opt_3_adapt) ? $opt_3_adapt : "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG",
		    defined($correct_5_prime_G) ? $correct_5_prime_G : "true", 
		    $opt_working_directory,
			$opt_gene_list,
			defined($cover_many_transcripts) ? $cover_many_transcripts : "true"
		);
}elsif($opt_task eq "library_assembly"){
	if(!defined($opt_gene_list)){
		print "\nThere must be an gene list specified.\n\n".$ver_str.$help_str;
	   exit;
	}
    if(!defined($opt_gene_list) or !(-f $opt_gene_list)){ die "The gene list file $opt_gene_list could not be opened. Either the user has no rights the read it or the file does not exist." }
    if(!defined($opt_param_file) or !(-f $opt_param_file)){ die "The parameter list file $opt_param_file could not be opened. Either the user has no rights the read it or the file does not exist." }
    if(!defined($opt_working_directory) or !(-d $opt_working_directory)){ $opt_working_directory="~" }
    if(!defined($opt_input_folder) or !(-d $opt_input_folder)){die "The input folder $opt_input_folder could not be opened. Either the user has no rights the read it or the file does not exist." }
    if(defined($opt_coverage)		and	$opt_coverage		=~m/([^\d]+)/g){ die "The gene coverage must be specified as an integer. Your input is not integer as it contains: $1." }
    if(defined($opt_total_lib_size)	and	$opt_total_lib_size	=~m/([^\d]+)/g){ die "The total library size must be specified as an integer. Your input is not integer as it contains: $1." }
    if(defined($opt_5_adapt)		and	$opt_5_adapt		=~m/([^ACGT]+)/g){ die "The 5' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
    if(defined($opt_3_adapt)		and	$opt_3_adapt		=~m/([^ACGT]+)/g){ die "The 3' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
    if(defined($correct_5_prime_G)	and	$correct_5_prime_G ne "true" and $correct_5_prime_G ne "false" ){ die "Correcting the 5' basepair to a G must be true or false" }    
    
    filter_library(
		    $opt_input_folder,#"Thu_Dec_4_17:28:18_20141417710498",
		    defined($opt_library_name) ? $opt_library_name : "test_lib",
		    defined($opt_coverage) ? $opt_coverage : 15,
		    defined($opt_total_lib_size) ? $opt_total_lib_size : 200000, 
		    defined($opt_5_adapt) ? $opt_5_adapt : "CTGAGCTCATAGAAGACCTCACC", 
		    defined($opt_3_adapt) ? $opt_3_adapt : "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG",
		    defined($correct_5_prime_G) ? $correct_5_prime_G : "true", 
		    $opt_working_directory,
			$opt_gene_list,
			defined($cover_many_transcripts) ? $cover_many_transcripts : "true"
		   );
}else{
	die "Some parameters are missing!\nOne of the four options make_database , target_ident , end_to_end or library assembly must be set.\n";
}


=head1 VERSION

Version 0.1.4

=cut

our $VERSION = '0.1.8';

=head1 SUBROUTINES/METHODS
=head2 calculate_CRISPR_score
#########################################################################################
#name:      calculate_CRISPR_score
#function:  helper-function for find_and_print_CRISPRS
#           creates the score for the given CRISPRS
#input:     (builded tree as referemce, something-Hashreference, start, end, chrom
#           1/0 for $score{"CDS"}++)
#output:    the calculated score as hash
#########################################################################################
=cut

sub calculate_CRISPR_score {
      my %score = ();
	  my $gene_name=$_[7];
      my @new_score=@{$_[6]};  
      my $expression="[";
      if (($_[1]->{"number_of_CDS"}>0)) {
         foreach my $number(1..$_[1]->{"number_of_CDS"}){
            $expression.=$number;
            }
         $expression.="]";
      }else{
             $expression="[.]";
      }
      $score{"CRISPRi"}=0;
      $score{"CRISPRa"}=0;
      my %transcripts=();
      if ( exists $_[0]->{$_[4]} ) { # check wethere the tree exists
            #search for annotations in the intervall from start (2) to end (3) and count them in score
            my $annotations = $_[0]->{$_[4]}->fetch( int($_[2]), int($_[3]) );
           
            foreach  my $anno ( @{$annotations} ) {
				
                  if ( $anno =~ m/gene_(\S+)_([0-9]+)_([0-9]+)/ ) {
                       my $gene_annotations = $_[0]->{$_[4]}->fetch( int($2), int($3) );                        
                        foreach  my $gene_anno ( @{$gene_annotations} ) {
                              if ($gene_anno=~m/exon::(\S+)::(\d+)::(\S+)\_(\d+)_(\d+)/) {    
                                    ${$transcripts{$1}}{"exon".$2}=$4."_".$5;
                              }
                        }
                        last;
                  }
            }
            foreach  my $anno ( @{$annotations} ) {
                  if ( $anno =~ m/gene_(\S+)_[0-9]+_[0-9]+/ ) {
						my $temp=$1;
                        $new_score[1]++;
						#print "$temp\t$gene_name\n";
						if ($temp=~m/$gene_name/) {
							${ $score{"gene"} }{$temp}++;
						}		
                  } elsif ( $anno =~ m/exon::(\S+)::(\d+)::(\S+)\_(\d+)_(\d+)/) {
                        ${ $score{"exon"} }{$2}++;
						${ $score{"gene_to_exon"} }{$3}++;
						$new_score[1]=$new_score[1]+5/$2;
                        if (exists $score{"transcripts"}) {
                              $score{"transcripts"}=$score{"transcripts"}."_".$1;
                        }else{
                              $score{"transcripts"}="_".$1;
                        }
                        ${$score{"transexons"}}{$1."::".$2}++;
                        
                  } elsif ( $anno =~ m/CpG/ ) {
                        $score{"CpG"}++;
			$new_score[1]--;
                  } elsif ( $anno =~ m/CDS::(\S+)::($expression)::(\S+)\_(\d+)_(\d+)$/ ) {
                        $new_score[1]=$new_score[1]+5/$2;
		       if($_[1]->{"specific_transcript"} ne "any"){
			   if ($1 eq $_[1]->{"specific_transcript"}) {
                                    $score{"CDS_1"}++;
                              }
                        } else{
                              $score{"CDS_1"}++;
                        }
                        if ($_[5] == 1) { # only for FORWARD needed @Flo ist das so gewollt?
                              $score{"CDS"}++;
                        }
                  } elsif ( $anno =~ m/CDS::(\S+)::(\d)::(\S+)\_(\d+)_(\d+)$/ ) {
                        $score{"CDS"}++;
						${ $score{"gene_to_CDS"}}{$3}++;
                        $new_score[1]++;
                  }
            }
            my $strand=1;
            my $first_start=0;
            my $first_end=0;
            my $second_start=0;
            foreach my $trans (keys %transcripts){
                  my %temp=%{$transcripts{$trans}};
				  if(exists $temp{"exon1"}){
					my @first_temp=split("_",$temp{"exon1"});
					$first_start=$first_temp[0];
					$first_end=$first_temp[1];
				  }
                  if(exists $temp{"exon2"}){
											my @second_temp=split("_",$temp{"exon2"});
											$second_start=$second_temp[0];
											if ($second_start<$first_start) {
												$strand=0;
										   }
					};
                  
                                 #print $_[1]->{"crispri_upstream"}."\t".$_[1]->{"crispri_downstream"}."\n";
                  if (
                        $strand==1
                        && ($first_start+$_[1]->{"crispri_upstream"})>=int($_[2])
                        && ($first_start-$_[1]->{"crispri_downstream"})<=int($_[2])
                        ) {
                          $score{"CRISPRi"}=1;
                     }elsif(
                           $strand==0
                         && ($first_end-$_[1]->{"crispri_upstream"})<=int($_[2])
                         && ($first_end+$_[1]->{"crispri_downstream"})>=int($_[2])
                     ){
                           $score{"CRISPRi"}=1;
                     }elsif(
                           $strand==1
                         && ($first_start-$_[1]->{"crispra_upstream"})<=int($_[2])
                         && ($first_start-$_[1]->{"crispra_downstream"})>=int($_[2])
                     ){
                           $score{"CRISPRa"}=1;
                     }elsif(
                           $strand==0
                         && ($first_end+$_[1]->{"crispra_upstream"})>=int($_[2])
                         && ($first_end+$_[1]->{"crispra_downstream"})<=int($_[2])
                     ){
                           $score{"CRISPRa"}=1;
                     }
				 # print  $first_start."\t".int($_[2])."\t".$score{"CRISPRi"}."\n";
            }
            
            #search for the start and stop coddon in the intervall from start (2) to end (3) with an up/downstream window
            $annotations = $_[0]->{$_[4]}->fetch( int( $_[2] - $_[1]->{"downstream_window"}), int( $_[3] + $_[1]->{"upstream_window"}) );
            foreach  my $anno ( @{$annotations} ) {
                  if ( $anno =~ m/start_codon/ ) {
                        $score{"start_codon"}++;
                        $new_score[1]++;
                  }elsif ( $anno =~ m/stop_codon/ ) {
                        $score{"stop_codon"}++;
                        $new_score[1]++;
                  }
            }            
      }
      $score{"new_score"}=\@new_score;
      return %score;
}
=head2 make_CRISPR_statistics
#########################################################################################
#name:      make_CRISPR_statistics
#function:  helper-function for find_and_print_CRISPRS
#           adjust statistics for the given CRISPRS
#input:     (something-Hashreference, score-Hashreference, dont_asses_context,
#           statistics-hashreference)
#output:    1 if delete and next in loop is needed, 0 if not
#########################################################################################
=cut
sub make_CRISPR_statistics {
      if ( $_[0]->{"gene_exclusive"} eq "true" && !exists $_[1]->{"gene"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they did not hit any gene"}++;
            return 1;
      }
      if ( $_[0]->{"exclude_overlapping_genes"} eq "true"  && scalar(keys(%{$_[1]->{"gene"}}))>1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they hit multiple genes"}++;
            return 1;
      }
      if (  $_[0]->{"exon_exclusive"} eq "true" && !exists $_[1]->{"exon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they did not hit any exon"}++;
            return 1;
      }
      if (  $_[0]->{"CpG_exclusive"} eq "true" && exists $_[1]->{"CpG"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were located in an CpG island"}++;
            return 1;
      }
      if (  $_[0]->{"purpose_exclusive"} eq "true" && ($_[0]->{"purpose"} eq "knockout") && !(exists $_[1]->{"CDS_1"}) && ($_[2] != 1) ) {
            $_[3]->{"Number of designs excluded because they were not directly behind the ATG of the specified transcript"}++;
            return 1;
      }
      if (  $_[0]->{"purpose_exclusive"} eq "true" && ($_[0]->{"purpose"} eq "n-tagging") && !exists $_[1]->{"start_codon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located at the start codon"}++;
            return 1;
      }
      if (  $_[0]->{"purpose_exclusive"} eq "true" && ($_[0]->{"purpose"} eq "c-tagging") && !exists $_[1]->{"stop_codon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located at the stop codon"}++;
            return 1;
      }
		if ( exists $_[0]->{"purpose_exclusive"} && ($_[0]->{"purpose"} eq "CRISPRa") && $_[1]->{"CRISPRa"}!=1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not amenable for CRISPRa"}++;
            return 1;
      }
      if ( exists $_[0]->{"purpose_exclusive"} && ($_[0]->{"purpose"} eq "CRISPRi") && $_[1]->{"CRISPRi"}!=1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not amenable for CRISPRi"}++;
            return 1;
      }
      if (  $_[0]->{"CDS_only"} eq "true" && !(exists $_[1]->{"CDS"}) && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located in a coding sequence"}++;
            return 1;
      }
      if (  $_[0]->{"exon_exclusive"} eq "true" && !( $_[0]->{"specific_exon"} eq "any") && $_[2] != 1  ) {
            my $counter = 0;
            foreach  my $exon ( keys (%{ $_[1]->{"exon"} }) ) {
                  if ( $exon == $_[0]->{"specific_exon"} ) {
                        $counter++;
                        last;
                  }
            }
            if ( $counter == 0 ) {
                  $_[3]->{"Number of designs excluded because they did not hit the specific exon"}++;
                  return 1;
            }
      }
      if ( !( $_[0]->{"specific_transcript"} eq "any") && $_[2] != 1 ) {
            my $counter = 0;
            foreach  my $transcript ( split("_",$_[1]->{"transcripts"} ) ) {
                  if ( $transcript eq $_[0]->{"specific_transcript"} ) {
                        $counter++;
                        last;
                  }
            }
            if ( $counter == 0 ) {
                  $_[3]->{"Number of designs excluded because they did not hit the specific transcript"}++;
                  return 1;
            }
      }
}
=head2 make_pos_index
#########################################################################################
#name:      make_pos_index
#function:  return an index of occurrences of a certain character in any given string
#input:     (string reference,character)
#output:    hash with character positions as keys
#########################################################################################
=cut
sub make_pos_index {
      my %pos     = ();
      my $result  = index( ${$_[0]}, $_[1], 0);
      while ( $result != -1 ) {
            $pos{$result}++;
            $result = index( ${$_[0]}, $_[1], ($result + 1) );
      }
      return %pos;
}
=head2 find_base_count
#########################################################################################
#name:      find_base_count
#function:  return an array containing the percentage each character in the string
#input:     (string)
#output:    array of percentages
#########################################################################################
=cut
sub find_base_count {
      my $seq  = $_[0];
      return      int( ($seq =~ tr/A/x/) * 100 / length($seq) ),
                  int( ($seq =~ tr/C/x/) * 100 / length($seq) ),
                  int( ($seq =~ tr/T/x/) * 100 / length($seq) ),
                  int( ($seq =~ tr/G/x/) * 100 / length($seq) );
}
=head2 reverse_comp
#########################################################################################
#name:      reverse_comp
#function:  return the reverse complement of a DNA sequence
#input:     (string)
#output:    reverse complement string
#########################################################################################
=cut
sub reverse_comp {
      ( my $rev = reverse $_[0] ) =~ tr/ACGTacgt/TGCAtgca/;
      return $rev;
}
=head2 mean
#########################################################################################
#name:      mean
#function:  return the mean of an array of number
#input:     (string)
#output:    mean as string
#########################################################################################
=cut
sub mean {
      return sum(@_) / @_;
}
=head2 variance
#########################################################################################
#name:      variance
#function:  return the variance of an array of number
#input:     (string)
#output:    variance as string
#########################################################################################
=cut
sub variance {
      return ( sum( map { ( $_ - mean( @{ $_[0] } ) )**2 } @{ $_[0] } ) / @{ $_[0] } );
}
=head2 gene_label
#########################################################################################
#name:      gene_label
#function:  return the last value in the features  primary tag (Bio::SeqFeature)
#input:     (Bio::SeqFeature::Generic->new)
#output:    label as string
#########################################################################################
=cut
sub gene_label {
      my $feature=shift;
      if ($feature->can("primary_tag") ) {
            my $notes=$feature->primary_tag;
            my @notes=split("::",$notes);
            $notes=pop(@notes);
            $notes;
      }else{
            my $notes="";
            $notes;
      }
}
=head2 make_my_id
#########################################################################################
#name:      make_my_id
#function:  return searchable ids for generation of the FASTA index
#input:     (FASTA header)
#output:    array of searchable ids (strings)
#########################################################################################
=cut
sub make_my_id {
      $_[0] =~m/^>(\S+) locus_tag= (\S+);/;
      return ( $1, $2);
}
=cut

=head2 build_tree
#########################################################################################
#name:      build_tree
#function:  if a pre.build tree exists, load this tree into memory if not build it 
#           from source and save to it to file and return the tree anyways
#input:     (path to data without file ending)
#output:    augmented black/red tree (Set::IntervalTree->new)
#########################################################################################
=cut

sub build_tree {
      my  $trees = Set::IntervalTree->new;
      if (-e $_[0].".otree") {
            $trees->LoadTree($_[0].".otree");
      }elsif(-e $_[0].".mygff")  {
            open (my $infile, "<", $_[0].".mygff");
                  foreach my $line (<$infile>) {
                        chomp $line;
                        my @line=split("\t",$line);
                        my $object= $line[0] . "_" . $line[1] . "_" . $line[2];
                        $trees->insert(  $object, $line[1], $line[2] );
                  }
            close($infile);
            $trees->SaveTree($_[0].".otree");
      }
      return $trees
}

=head2 round_digits
#########################################################################################
#name:      round_digits
#function:  cutoff after certain number of digits
#input:     (float number as string, digits after komma <int>) 
#output:    shortened float
#########################################################################################
=cut
sub round_digits{
      if($_[0]=~m/(\d+\.\d{$_[1]}).*/){
            return $1;
      }
      return $_[0];
}

=head2 make_mismatch_string
#########################################################################################
#name:      make_mismatch_string
#function:  convert a SAM file format mismatch string into an "Mismatch string" 
#           M for every match X for every mismatch D for every deletion I for every 
#           insertion
#input:     (Sequence as string-reference, $something{"unspecific_leading_bases"},
#           direction) 
#output:    shortened float
#########################################################################################
=cut
sub make_mismatch_string{
      my $mismatchstring      = "";
      my $pos                 = 0;
      my @stringarray=split("\t",${$_[0]});
      if(${$_[0]}=~m/MD:Z:(\S+)\s/){
            $mismatchstring=$1;
      }
      my @matchstring=split("",$stringarray[9]);
      my @matches=$stringarray[5]=~m/[0-9]+[MID]/g;
      foreach my $match (@matches){
            $match=~/([0-9]+)([MID])/;
            foreach (1..$1){
                  $matchstring[$pos]=$2;
                  $pos++;
            }
      }
      @matches    = $mismatchstring =~m/[0-9]+|[\^A-Z]+|[0-9]+$/g;
      $pos        = 0;
      foreach my $match (@matches){
            if($match=~/([0-9]+)/){
                  $pos += $1;
            }elsif($match=~/^[A-Z]$/){
                  $matchstring[$pos]="X";
                  $pos++;
            }else{
                  $pos++
            }
      }
      if ($_[2] eq "fw") {
            foreach (1..$_[1]){
                  unshift @matchstring , "n";
            }
      }
      else {
            foreach (1..$_[1]){
                  push @matchstring , "n";
            }
      }
      
      
      return(@matchstring);
}
=head2 print_offtarget_string
#########################################################################################
#name:      print_offtarget_string
#function:  convert a "Mismatch string" into html output with different color-highlighting
#           for M,X,D and I
#input:     (Mismatch string) 
#output:    converted string for html output
#########################################################################################
=cut
sub print_offtarget_string {
      my $string=$_[0];
      $string=~s/n/<span style="color: black;">n<\/span>/g;
      $string=~s/M/<span style="color: lightgreen;">M<\/span>/g;
      $string=~s/X/<span style="color: red;">X<\/span>/g;
      $string=~s/D/<span style="color: pink;">D<\/span>/g;
      $string=~s/I/<span style="color: pink;">I<\/span>/g;
      return $string;
}
=head2 score_micro_homolgy
#########################################################################################
#name:      score_micro_homolgy
#function:  calculate an microhomology score between 0 and 12 for a qiven position in a 
# indexed sequence
#input:     (indices hash of hashes, threshold int,length int, sequence ref to string ) 
#output:    micro homolgy score (int)
#########################################################################################
=cut
sub score_micro_homolgy {
	my $score=0;
	my $count=0;
	my $lengthhom=1;
	my $seq;
	my $stuff=0;
	my $outframe=1;
	my $inframe=1;
	my $limit_right=$_[2]+$_[1];
	if(($_[2]+$_[1]+$_[3])>length(${$_[4]})){
	    $limit_right=length(${$_[4]})-$_[3];
	}
	RIGTHSTART:foreach my $rightarmstart ($_[2]..($limit_right)){		
		LENGTH:foreach my $length (2..$_[3]){
			if($lengthhom>0){
				$lengthhom=0;
				my @right_seq=split("",substr(${$_[4]},$rightarmstart,$length));
				$stuff=$length+1;
				LEFTARM: while($stuff<$_[1]){
					$stuff++;
					$count=0;
					$seq="";
					foreach my $letter (@right_seq){
						if(exists(${${$_[0]}{$letter}}{($_[2]-$stuff+$count)})){
							$seq.=$letter;
							$count++;												
						}else{
							next LEFTARM;					
						}
					}
					if($count==scalar(@right_seq)){
						$lengthhom=1;
						my $gaplength=($stuff+($rightarmstart-$_[2]));
						if($gaplength%3 == 0){
							$inframe=$inframe+($count*exp(0.1*(-$gaplength)));
						}else{
							$outframe=$outframe+($count*exp(0.1*(-$gaplength)));
						}
						#print "(@right_seq): $seq ".scalar(@right_seq)." ".$count." $stuff $position $rightarmstart hit"."\n";
						#print substr($string,($position-20),40)."\n";
						#print substr($string,($position-20),(20-$stuff+$count));
						#print '-' x ($stuff+($rightarmstart-$position-$count)); 
						#print substr($string,($rightarmstart),20-($rightarmstart-$position))."\n";
					}else{
						next LEFTARM;
					}
				}	
			}else{
				$lengthhom=1;
				next RIGTHSTART;		
			}	
		}
	}
	return (10*log($outframe/$inframe));
}
=head2 hor_bar_chart
#########################################################################################
#name:      hor_bar_chart
#function:  print horizontal bar chart
#input:     (@array(score 1, score 2, score 3))
#output:    N/A
#########################################################################################
=cut
sub hor_bar_chart {
	my $outstring='
	 <table>
      <tbody>
        <tr>
          <td style="font-size: 0.6em;">S</td><td class="bar"><div style="width: '.$_[0].'px;" class="item1"></div></td>
           </tr><tr>
          <td style="font-size: 0.6em;">A</td><td class="bar"><div style="width: '.$_[1].'px;" class="item2"></div></td>
           </tr><tr>
          <td style="font-size: 0.6em;">E</td><td class="bar"><div style="width: '.$_[2].'px;" class="item3"></div></td>
        </tr>
      </tbody>
    </table> ';
    return($outstring);
}
=head2 filter_library
#########################################################################################
#name:      filter_library
#function:  filter library and print orderable fasta files
#input:     ()
#output:    N/A
#########################################################################################
=cut
sub filter_library{
    my $temp_dir=$_[0];
    my $lib_name=$_[1];
    my $coverage=$_[2];
    my $limit=$_[3];
    my $five_prime_extension=$_[4];
    my $three_prime_extension=$_[5];
    my $correct_five_prime=$_[6];
	my $gene_list_file=$_[8];
	my $many_transcripts=$_[9];
    my @info=();
    my @line=();
    my %designs=();
    my %genes=();
	my %genes_from_list=();
	
	my %id_for_lib;
    my $gene="";
    my $version=2;
	open (my $gene_list, "<", $gene_list_file) or die $!;
		while(<$gene_list>){
			chomp $_;
			$genes_from_list{$_}++;
		}
	close $gene_list;
	my %genes_avail;
    open (my $libgff, "<", $temp_dir . "/all_results_together.gff") or die $!;
       while(<$libgff>){
            if(!($_=~/\#/)){
				my $cline = $_;
                @line=split("\t",$cline);
                @info=split(";\ ",$line[8]);
				foreach my $key (keys(%genes_from_list)){
					if ($cline=~/$key/) {
						$gene=$key;
						$genes_avail{$gene}++;
					}
					
				}
                my $id="";
                foreach my $number (0..5) {
                    $info[$number]=~m/(.+?)=(.*)/;
                    if($number==0){
                        $id=$2;
                    }else{
                        ${$designs{$id}}{$1}=$2;
                        ${$designs{$id}}{"gene"}=$gene;
                    }
                }
            }
        }
    close $libgff;
    my $current_size=0;
    my $last_gene="";
    my @new_file=(); 
    my %ids=();
	my @ids=();
    my %count=();
    my %missing=();
	my %all_ids;
    my $general_coverage=0;
	my %id_for_lib;
	my %id_with_info;
	open (my $libtab, "<", $temp_dir . "/all_results_together.tab") or die $!;
	while(<$libtab>){		
        @line=split("\t",$_);
		my @transcripts=split("_",$line[7]);
		foreach my $trans (@transcripts){
			${${$ids{$line[6]}}{$trans}}{$line[0]}++;
		}
		${${$id_with_info{$line[6]}}{$line[0]}}{"anno_score"}=$line[16];
		${${$id_with_info{$line[6]}}{$line[0]}}{"spec_score"}=$line[15];
		${${$id_with_info{$line[6]}}{$line[0]}}{"eff_score"}=$line[17];
		${${$id_with_info{$line[6]}}{$line[0]}}{"seq"}=$line[5];
	}
	close $libtab;
	my %isgone;
	foreach my $key (keys %ids){
			my %total=();
				foreach my $subkey (keys %{$ids{$key}}){
					foreach my $subsubkey (keys %{${$ids{$key}}{$subkey}}){
						$total{$subsubkey}++
					}
				}
				if ((scalar (keys %total)) < $coverage) {
					delete $ids{$key};
					$isgone{$key}=(scalar (keys %total));
				}
				
		}
	if ($many_transcripts eq "true") {
		my $temp="";		
		foreach my $key (keys %ids){
			while ((scalar (keys %{$id_for_lib{$key}})) < $coverage) {
					foreach my $subkey (
					sort {scalar keys %{${$ids{$key}}{$a}} <=> scalar keys %{${$ids{$key}}{$b}} } keys %{$ids{$key}}
					){					
					if (defined((keys(%{${$ids{$key}}{$subkey}}))[0]) && ((scalar (keys %{$id_for_lib{$key}})) < $coverage)) {
						$temp= (keys(%{${$ids{$key}}{$subkey}}))[0];
						delete ${${$ids{$key}}{$subkey}}{$temp};
						${$id_for_lib{$key}}{$temp}++;
					}	
				}
			}			
		}
	}else{
		
		foreach my $key (keys %ids){
			while ($count{$key} < $coverage) {
				foreach my $element (
					sort { ${$id_with_info{$key}}{$b}->{"anno_score"} <=> ${$id_with_info{$key}}{$a}->{"anno_score"} }
					sort { ${$id_with_info{$key}}{$b}->{"spec_score"} <=> ${$id_with_info{$key}}{$a}->{"spec_score"} }
					sort { ${$id_with_info{$key}}{$b}->{"eff_score"} <=> ${$id_with_info{$key}}{$a}->{"eff_score"} }
					keys %{$id_with_info{$key}}
					){
					if( 		$count{$key} < $coverage &&
									!(${${$id_with_info{$key}}{$element}}{"seq"}=~m/GAAGAC/) &&
								   !(${${$id_with_info{$key}}{$element}}{"seq"}=~m/GTCTTC/) &&
								   !(${${$id_with_info{$key}}{$element}}{"seq"}=~m/GAATTC/) &&
								   !(${${$id_with_info{$key}}{$element}}{"seq"}=~m/CTTAAG/) &&
								   !(${${$id_with_info{$key}}{$element}}{"seq"}=~m/CAATTG/) &&
								   !(${${$id_with_info{$key}}{$element}}{"seq"}=~m/GTTAAC/) &&
								   !(${${$id_with_info{$key}}{$element}}{"seq"}=~m/CTCGAG/) &&
								   !(${${$id_with_info{$key}}{$element}}{"seq"}=~m/GAGCTC/) 
						){
							$count{$key}++;
							${$id_for_lib{$key}}{$element}++;
					}else{
						delete ${$id_with_info{$key}}{$element};
					}
			   }
			}
		}
	}
	
	foreach my $key (keys %id_for_lib){		
		if ($general_coverage>$limit) {
			delete $id_for_lib{$key};
		}else{
			$general_coverage+=scalar keys %{$id_for_lib{$key}}
		}		
	}
	
	#foreach my $key (keys %id_for_lib){
	#	print $key."\t".join("||",keys %{$id_for_lib{$key}})."\n";
	#}

    $lib_name=$lib_name.".$version."."$coverage.".scalar(keys(%id_for_lib));
    open (my $libtab_out, ">", $lib_name.".tab") or die $!;
    open (my $libfa_out, ">", $lib_name.".fasta") or die $!;
    open (my $libtab, "<", $temp_dir . "/all_results_together.tab") or die $!;
    my %fasta=();
        while (<$libtab>) {
            @line=split("\t",$_);
            if (${$id_for_lib{$line[6]}}{$line[0]}) {
				$all_ids{$line[0]}++;
				if (!exists($fasta{$line[0]})) {
					print $libtab_out $_;
					if ($correct_five_prime eq "true") {
						if ($line[5]=~m/\w(\w+)\s[NACGT]+_\w(\w+)\s[NACGT]+/) {
							print $libfa_out ">".$line[0]."_left\n".$five_prime_extension."G".$1."$three_prime_extension\n";
							print $libfa_out ">".$line[0]."_right\n".$five_prime_extension."G".$2."$three_prime_extension\n";
						}else{
							$line[5]=~m/\w(\w+)\s\w+/;
							print $libfa_out ">".$line[0]."\n".$five_prime_extension."G".$1."$three_prime_extension\n";
						}
					}else{
						if ($line[5]=~m/(\w+)\s[NACGT]+_(\w+)\s[NACGT]+/) {
							print $libfa_out ">".$line[0]."_left\n".$five_prime_extension.$1."$three_prime_extension\n";
							print $libfa_out ">".$line[0]."_right\n".$five_prime_extension.$2."$three_prime_extension\n";
						}else{
							$line[5]=~m/(\w+)\s\w+/;
							print $libfa_out ">".$line[0]."\n"."$five_prime_extension".$1."$three_prime_extension\n";
						}
					
					}
					
					$fasta{$line[0]}++;
				}
            }
            
        }
    close ($libtab);
    close ($libfa_out);
    close ($libtab_out);    
    open ($libgff, "<", $temp_dir . "/all_results_together.gff") or die $!;
    open (my $libgff_out, ">", $lib_name.".gff") or die $!;
        foreach my $line (<$libgff>) {
            if($line=~m/id\=(.+?)\;/){
                if ($all_ids{$1}) {
                    print $libgff_out $line;     
                }
            }else{
                print $libgff_out $line;  
            }
        }
    close $libgff;
    close $libgff_out;
    open ($libtab_out, ">", $lib_name.".large.tab") or die $!;
    open ($libtab, "<", $temp_dir . "/all_results_together.tab") or die $!;
        while (<$libtab>) {
                print $libtab_out $_;
        }
    close $libtab_out;
    close $libtab;
    open ($libgff_out, ">", $lib_name.".large.gff") or die $!;
    open ($libgff, "<", $temp_dir . "/all_results_together.gff") or die $!;
        while (<$libgff>) {
                print $libgff_out $_;
        }
    close $libgff_out;
    close $libgff;
    open (my $coverage_file, ">", $lib_name.".coverage.tab") or die $!;   
        foreach my $key (keys %id_for_lib) {
            print $coverage_file $key."\t".$count{$key}."\n";
			delete $genes_from_list{$key};
        }
    close $coverage_file;
    open (my $mis, ">", $lib_name.".missing.tab") or die $!;   
        foreach my $otherkey (keys %id_for_lib) {
			my $temp=0;
			foreach my $key (keys %genes_from_list) {
				if ($otherkey=~m/$key/) {
					$temp=1;
				}
			}
			if ($temp==0) {
				print $mis $otherkey."is missing from the library. It was covered by ".$isgone{$otherkey}." designs. Maybe it was covered two low or not found in the cld database.\n";
			}
			
        }
		print $mis (scalar keys %genes_from_list) - (scalar keys %id_for_lib )."are missing ebcause of two harsh design citeria\n.";
    close $mis;    
}


=head2 create_popup
#########################################################################################
#name:      create_popup
#function:  creates popup with. the target-, match- and querry-string adjusted to each
#input:     (line of the outfile.tab as string, databasepath as string,
#           $something{"unspecific_leading_bases"}, int overflow before and after Sequence)
#output:    HTML formated string
#########################################################################################
=cut
sub create_popup {
      my @line          = @{$_[0]};
      my $report        = "";
      my $target        = "";
      my $match         = $line[17];
      my $querry        = "";
      my $adjust        = 0;
      my $adjustleft    = 0;
      my $start         = 0;
      my $end           = 0;
      my $insert        = 0;
      
      if ($line[20] eq "rc") { #turn and translate the match and querry string if direction is reverse complementary
            $querry = reverse_comp($line[5]);
            $querry =~s/\s+//ig;
            $querry =~s/^(\w{2})\w(\w+)/$1N$2/ig;
      }
      else {
            $querry = $line[5];
            $querry =~s/\s+//ig;
             $querry =~s/(\w+)\w(\w{2})$/$1N$2/ig;
            $adjust = $_[2];
      }
      
      #search the target-string in database - substr workaround used, because the Bioperl function does not work proper with large numbers
      my $db            = Bio::DB::Fasta->new( $_[1] . '.all.dna.fa', -makeid => \&make_my_id );
      my @target_name   = split("::", $line[14]);
      my $obj           = $db->get_Seq_by_id($target_name[0]);
      if (!defined $obj) { #have to search in whole chromosom
            $db         = Bio::DB::Fasta->new( $_[1] . '.dna.toplevel.fa', -makeid => \&make_my_id );
            $obj        = $db->get_Seq_by_id($line[14]);
      }
      $start   = $line[15]-$adjust-$_[3]-1;
      if ($start < 0 ) {
            $adjustleft = abs($start);
            $start = 0;
      }
      
      $end     = ($line[16]-$adjust+$_[3]);
      $target  = substr $obj->seq, $start, $end - $start - 1;
      #count insertions and adjust the end property to display
      $insert  = $line[17] =~ tr/I/x/;
      
      $report .= '<table>';
      $report .= '<tr><td>Target:</td><td>|'.$start.'*|</td>'.create_popup_string($target, $_[3], 0, 0, $match, 0, 0, $adjustleft).'<td>|'.($end-$insert).'*|</td></tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match, $_[3], 1, 1, $match, 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry, $_[3], 1, 0, $match, 0, 0, $adjustleft).'</tr>';
      $report .= '</table>';
      $report .= '<hr /><p style="font-size:50%">*Start and End-point in the target</p>';
      
      return $report;
}
=head2 create_popupD
#########################################################################################
#name:      create_popupD
#function:  for double sequence - does the same as create_popup
#input:     (line of the outfile.tab as string, databasepath as string,
#           $something{"unspecific_leading_bases"}, int overflow before and after Sequence)
#output:    HTML formated string
#########################################################################################
=cut
sub create_popupD {
      my @line          = @{$_[0]};
      my $report        = "";
      my $target        = "";
      my @match         = split("-", $line[17]);
      my $tmp           = $line[5];
      $tmp              =~s/\s+//ig;
      my @querry        = split("_", $tmp );
      my $leadingbases  = $_[2];
      my $adjustleft    = 0;
      my $start         = 0;
      my $end           = 0;
      my $insert0       = 0;
      my $insert1       = 0;
      my $deletion0     = 0;
      my $deletion1     = 0;
      my $spacer        = "";
      my $match0        = "";
      my $match1        = "";
      
      #adjust the first/second querry-string 
      if ($line[20] eq "rc") {
            $querry[0]   = reverse_comp($querry[0]);
            $match0 = $match[0];
            $match1 = $match[1];
            #count insertions and deletions to adjust the end/start property
            $deletion0   = $match0 =~ tr/D/o/;
            $deletion1   = $match1 =~ tr/D/o/;
            #calc start and end of the target sequence
            $start       = $line[15] - $_[3] - 1 - $deletion1; # start of the match - the wanted overflow on the left side
            $end         = $line[16] + $line[21] + length($match[0]) + $_[3] - ($leadingbases*2) - ($deletion0*2); # end of the match + spacer size + length of the matchstring + the wanted overflow on the right side - the count of unspec leadingbases*2 (1st and 2nd)
      }
      else {
            $querry[1]   = reverse_comp($querry[1]);
            #swap both query and matchstring and change the start and end points ($line[15] is no longer the start of the target cause the first query swaped)
            @querry[0,1] = @querry[1,0];
            @match[0,1]  = @match[1,0];
            $match0 = $match[0];
            $match1 = $match[1];
            #count deletions to adjust the end property
            $deletion0   = $match0 =~ tr/D/o/;
            $deletion1   = $match1 =~ tr/D/o/;
            #calc start and end of the target sequence
            $start       = $line[15] - $line[21] - length($match[0]) - $_[3] + $leadingbases - 1 - $deletion1; # swaped of the if-clause and reverted (-/+), but only the leadingbases of the 2nd string count here
            $end         = $line[16] + $_[3] - $leadingbases - ($deletion0*2); # swaped of the if-clause and reverted (-/+), but only the leadingbases of the 1st string count here
      }
      
      if ($start < 0 ) {
            $adjustleft = abs($start);
            $start = 0;
      }
      
      #count insertions and deletions to adjust the end/start property
      $insert0   = $match0 =~ tr/I/o/;
      $insert1   = $match1 =~ tr/I/o/;
      for (my $i = 0; $i < ($line[21] - $leadingbases*2 - $deletion0 + $deletion1); $i++){
            $spacer .= " ";
      }
      
      #search the target-string in database - substr workaround used, because the Bioperl function does not work proper with large numbers
      my $db             = Bio::DB::Fasta->new( $_[1] . '.all.dna.fa', -makeid => \&make_my_id );
      my @target_name    = split("::", $line[14]);
      my $obj            = $db->get_Seq_by_id($target_name[0]);
      if (!defined $obj) { #have to search in whole chromosom
            $db          = Bio::DB::Fasta->new( $_[1] . '.dna.toplevel.fa', -makeid => \&make_my_id );
            $obj         = $db->get_Seq_by_id($line[14]);
      }
      
      $target  = substr $obj->seq, $start, $end - $start - 1;
      
      #build the table for the popup
      $report .= '<table style="font-size:75%">';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry[0], $_[3], 1, 0, $match[0], 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match[0], $_[3], 1, 1, $match[0], 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Target:</td><td>|'.$start.'*|</td>'.create_popup_string($target, $_[3], 0, 0, $match[0].$spacer.$match[1], 0, 1, $adjustleft).'<td>|'.($end - 1 - $insert1).'*|</td></tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match[1], $_[3], 1, 1, $match[1], (length($match[0]) + int($line[21]) + $insert0 - ($leadingbases*2) - $deletion0 + $deletion1), 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry[1], $_[3], 1, 0, $match[1], (length($match[0]) + int($line[21]) + $insert0 - ($leadingbases*2) - $deletion0 + $deletion1), 0, $adjustleft).'</tr>';
      $report .= '</table>';
      $report .= '<hr /><p style="font-size:50%">*start- and endpoint in the target sequence</p>';
      
      return $report;
}
=head2 create_popup_string
#########################################################################################
#name:      create_popup_string
#function:  adjusts a sequence for the popup-output
#input:     (string sequence, int overflow from create_popup, boolean if it should be used,
#           boolean if colored or not, string matchstring, 2nd overflow for double seq,
#           boolean if the last letter should  only be popped for 2nd half - 1 = yes,
#           int additional left re-adjustment)
#output:    HTML formated string
#########################################################################################
=cut

sub create_popup_string {
      my $string  = "";
      my @letters = split("",$_[0]);
      my @match   = split("",$_[4]);
      my $indel   = 0;
      my $letter  = "";
      my $length  = length($_[4]);
      
      if ($_[2] == 1) {
            for (my $i = 0; $i < $_[1] + $_[5] - $_[7]; $i++){
                  $string .= '<td></td>';
            }
      }
       
      
      for(my $i = 0; $i < scalar(@letters); $i++){
            if ($_[3] == 1) { # matchstring needs color
                  $string .= '<td align="center" valign="middle">'.print_offtarget_string($letters[$i]).'</td>';
            }
            else{
                  #adjust the letters for indel
                  $letter = $letters[$i];
                  if ($_[2] == 1 && $indel < $length) { # for querystring (same size as matchstring)
                        if ($match[$indel] eq "D") {
                              $letter = "_";
                              $i--;
                        }
                  }
                  elsif ($indel - $_[1] >= 0 && ($indel - $_[1]) < $length) { # for targetstring (needs adjustment)
                        if ($match[$indel - $_[1]] eq "I") {
                              $letter = "_";
                              $i--;
                              if ($_[6] == 0) {
                                    pop (@letters); #delete last Element
                              }
                              
                              elsif ($indel - $_[1] > $length/2) {
                                    pop (@letters); #delete last Element only for I's in the 2nd matchstring
                              }
                        }
                  }
                  
                  $string .= '<td align="center" valign="middle">'.$letter.'</td>';
                  $indel++;
            }
      }
      
      return $string;
}

=head2 make_temp_fasta_file
#########################################################################################
#name:      make_temp_fasta_file
#function:  creates a temporary fasta file for the bowtie index and builds a trees
#input:     (given id-Array, tree as referemce, something-Hashreference,
#           enzyme db, temp_dir, 1/0 if file or not)
#output:    N/A
#########################################################################################
=cut

sub make_temp_fasta_file {
      my %something=%{$_[2]};
      if ( !( $_[2]->{"specific_transcript"} eq "any") && scalar(@{$_[0]}) >1) {
            die "Transcript specificity is only defined for single gene analyses.\n";
      }
      open (my $tempfile, ">", $_[4] . "/tempfile.fasta");
            foreach my $id (@{$_[0]}) { 
                  $id =~ s/\s//ig;
                  my $seq_obj = $_[3]->get_Seq_by_id($id); # get a PrimarySeq obj
                  if ($seq_obj) {
                        my $header = $_[3]->header($id); # get the header, or description line
                        $header =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig;
                        my $chrom = $1;
                        my $location_offset = $2;
                        my $location_end = $3;
                        if ( !exists $_[1]->{$chrom} ) {
                              $_[1]->{$chrom} = build_tree( $something{"databasepath"} . $_[2]->{"ref_organism"} . "/" . $chrom . "_indexed" );
                        }
                        print $tempfile ">", $seq_obj->display_id(), " ", $header, "\n", $seq_obj->seq(), "\n";
                  } else {
                        open (my $failfile, ">>", $_[4] . "/failfile.tab");
                              print $failfile substr($id,0)."\n";
                        close($failfile);
                  }
            }
      close $tempfile;
      if ( open (my $failfile, "<", $_[4] . "/failfile.tab") && !( $_[2]->{"ignore_missing_id"} eq "true")) {
            my $error="";
            while (<$failfile>) {
                  $error.=$_."\n";
            }
            close($failfile);
            die "No Database entry found for <br> \"". $error."\" in the \" ".$_[2]->{"ref_organism"}."\" genome.\n Please enter a valid ensembl ID or gene symbol (case sensitive) or change the input option above to FASTA sequence.\n";;
      }
      
}
=head2 make_a_crispr_library
#########################################################################################
#name:      make_a_crispr_library
#function:  makes a large scale crispr library file from a list of gene names and a parameter file
#input:     (gene_list <string>, parameterfile <string>, database <string>)
#output:    table of CRISPRS and hash for statistics
#########################################################################################
=cut

sub make_a_crispr_library{
      #define default starting variables
      my $temp_dir            = "";
      my %trees               = ();
      my $seqio_obj           = "";
      my %something           = ();
	  my %weights=();
	  if (defined $scoreweight_file) {
		open(my $scorefile, "<", $scoreweight_file) or die $!;
			while (<$scorefile>) {
				if ($_=~m/$(.+)=(\d+)/) {
					$weights{$1}=$2;
				}				
			}			
		close($scoreweight_file);
	  }
	  
      my $parallel_number     = 2;
      $something{"input_file"}=$_[1];
      #create a time stamped output folder
      $temp_dir = scalar localtime();
      $temp_dir =~ s/\s+/\_/ig;
	  $temp_dir =~ s/\W+/\_/ig;
      $temp_dir=$_[2].$temp_dir;
      mkdir( "$temp_dir") or die $!;
      #read in the parameter file, must contain each parameter with name and value '=' separated, without any quotes
      open(my $parameterfile, "<", $_[0]) or die $_[0]."could not be opened. No such file or directory.";
            foreach my $element (<$parameterfile>) {
                  chomp $element;
                  my @line=split("=",$element);
                  $something{$line[0]} = $line[1];
            }
      close $parameterfile;
      
      my $databasepath = $something{"databasepath"} . $something{"ref_organism"} . "/" . $something{"ref_organism"};
      #################################################################################################################################################################################
      # upload a file and save it in the temp directory for bowtie
      #################################################################################################################################################################################
      if ($something{"sec_off_target"} eq "true") {
			if (-e $something{"databasepath"} .'secondary_off_targets.fasta') {
				system('bowtie2-build '.$something{"databasepath"}.'/secondary_off_targets.fasta '.$temp_dir.'/temp_sec ;');
			}else{
				die $something{"databasepath"} .'secondary_off_targets.fasta'."could not be opened. No such file or directory.";
			}
      }
      #################################################################################################################################################################################
      # For ENSEMBLE: define the path to the bowtie index and do checks if an database entry is found for the ensemble accesion number - if all checks passed, create the $seqio_obj
      #################################################################################################################################################################################
      if ( $something{"data_type"} eq "ensemble_acc" ) {
            my $db = Bio::DB::Fasta->new( $databasepath . '.all.dna.fa', -makeid => \&make_my_id );
            my @ids = ();
            open (my $infile, "<", $something{"input_file"});
                  while (<$infile>) {
                        my $line = $_;
                        chomp $line;
                        push @ids, $line;
                  }
            close $infile;
            make_temp_fasta_file(\@ids, \%trees, \%something, $db, $temp_dir, 1);
            $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ); #read the temporary fasta file            
      } else { 
            #################################################################################################################################################################################
            # For FASTA: define the path to the bowtie index and do checks if the fasta sequence/file is in the right format - if all checks passed, create the $seqio_obj
            #################################################################################################################################################################################
            my $count=0;
            my $temp="";
            open(my $infile, "<",$something{"input_file"});
                        while (my $line = <$infile>){
                              if ($line=~m/^(>.+)/) {
                                    $count++;
                              }elsif ($line=~m/([^ACGTUN\s]+)/){
                                    die $something{"input_file"}." is not a FASTA format file because it contains \"$1\" as bases\n" ;
                              }
                              $temp=$temp.$line;
                        }
            close $infile;
            if ($temp=~m/^(>[^\n]+)$/) {
                      die "\"$1\" is not a FASTA format sequence.\nFASTA sequences need a header starting with \">\" which is new line separated from the sequence.\n e.g.: >some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG\n";
                  
            }
            if($count==0){
                        die "A FASTA format sequence needs to have a header starting with \">\"\n e.g.: >some Sequence\nAGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG\n";
                        
            }
            $seqio_obj = Bio::SeqIO->new( -file => $something{"input_file"}, -format => "fasta" ) or die "Your input wasn't FASTA format \n e.g.: >some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG\n"; #if neither online or pasted sequences are used it will use the input file as sequence input
      }
      
      #################################################################################################################################################################################
      # Start the creation of the report (index.html)
      #################################################################################################################################################################################
      #define empty object to be filled in the process
      my $fname               = "";
      my $dont_asses_context  = 0;
      if ($something{"draw_html_report"} eq "true" ){
            open( my $report, ">", "$temp_dir/index.html" ) or die "can not open report file";
                    print $report '<!DOCTYPE html>
										<html lang="enc">
										<head>
										<meta http-equiv="Content-Type" content="text/html;charset=utf-8">
										<title>cld Results</title><tr><td><br><br>
                                 <a href="all_results_together.tab" target="_blank" download="all_results_together.tab"><input type="button" align="center" value="Download a tabular report for all query sequences together"></a>
                                 <br><br>';
		  if(-e "$temp_dir/index.html"){ print $report '<a href="failfile.tab" target="_blank" download="failfile.tab"><input type="button" align="center" value="Download a table of failed ids"></a>';}
		  print $report '<br><br><a href="libraryfile.tab" target="_blank" download="libraryfile.tab"><input type="button" align="center" value="Download the filtered library table"></a><br><br><a href="library.fasta" target="_blank" download="library.fasta"><input type="button" align="center" value="Download the ready for ordering fasta file"></a></td></tr>';
		  print $report '<tr><td colspan="2" rowspan="1" style="vertical-align: top;"><embed style="width: 100%; height: 10px;" alt="there should appear a  line" src="/cld/gradientline.svg" type="image/svg+xml"><br></td></tr>';
                  close $report;
                  chmod 0755, $temp_dir . "/index.html";
      }
      #################################################################################################################################################################################
      #start the main loop which is looping through the sequence found in the SeqIO object, may be one or many more
      #################################################################################################################################################################################
      my %statistics    = ();
      my %CRISPR_hash   = ();
      my %CRISPR_cnt    = ();
      my @seq_array     = ();
      while( my $seq = $seqio_obj->next_seq() ) {
            push(@seq_array,\$seq);
      }
      my $seqcount=0;
      my $maxcount=scalar(@seq_array);
      my @fname_array=();
      foreach my $seq_obj ( @seq_array ) {
            $seq_obj=$$seq_obj;
            my %statistics    = ();
            my %CRISPR_hash   = ();
            #find the genomic coordinates of the sequence either directly from the fasta sequence or from other smyces
            my $chrom               = ""; #create an empty name chromosome
            my $location_offset     = 0;
            my $location_end        = 0; #set the location offset value to 0
            if ( ( $seq_obj->description ) =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig) { #if the descrition of the sequence contains information of the form loc=chr1:0..13000000
                  #print the result of the pattern matching
                  $chrom = $1; #save the location information in the certain objects
                  $location_offset = $2;
                  $location_end = $3; #the location offset is the start of the sequence on that certain chromosome, for the whole chromosome this is of cause 0
                  if ( !exists $trees{$chrom} ) {
                        $trees{$chrom} = build_tree( $something{"databasepath"} . $something{"ref_organism"} . "/" . $chrom . "_indexed" );
                  }
            } else {
                  $dont_asses_context = 1;
            }
            #define the name for the current sequence as its display id
            $fname                                                = ( $seq_obj->display_id );#current name is the sequence' id
            $statistics{$fname}{"seq_name"}                       = $fname; 
            $statistics{$fname}{"seq_length"}                     = $seq_obj->length;
            $statistics{$fname}{"seq_location"}                   = $chrom."::".$location_offset."::".$location_end;
            $statistics{$fname}{"Number of successful designs"}   = 0;
            if ( $chrom eq "" ) { $chrom = $fname; } #if $chrom is still empty fill it with the sequence' id
            my $whole_seq = $seq_obj->seq; #deduce the complete nucleotide sequence as a alphanumerical string from the SeqIo object
            if ($statistics{$fname}{"seq_length"}!=0) {
                  if (int(log($statistics{$fname}{"seq_length"})) <= 8) {
                        if (int(log($statistics{$fname}{"seq_length"}))>=2) {
                              $parallel_number=int(log($statistics{$fname}{"seq_length"}));
                        }else{
                             $parallel_number=2;
                        }
                  }else{
                        $parallel_number=$max_parallel;
                  }
            }else{
                  $parallel_number=2;
            }
            ###########################################################################################################################################################################
            #create the hashes for the CRISPR and the statistics
            ###########################################################################################################################################################################
            if (!(-e $temp_dir . "/" .$fname . '.json')) {
                  if (!exists $CRISPR_hash{$fname}) { #only execute following calucaltions if it is not already done for the sequence
                        my @findings = find_and_print_CRISPRS(    \$seq_obj,
                                                                  $chrom,
                                                                  $location_offset,
                                                                  \%trees,
                                                                  $dont_asses_context,
                                                                  $temp_dir,
                                                                  $parallel_number,
                                                                  \%something,
                                                                  $location_end);
                        %{$CRISPR_hash{$fname}} = %{$findings[0]};
                        %{$statistics{$fname}} = (%{$statistics{$fname}},%{$findings[1]});
                  }
                  my $json = JSON::XS::encode_json(\%CRISPR_hash);
                  write_file( $temp_dir . "/" .$fname . '.json', { binmode => ':raw' }, $json );
                  $json = JSON::XS::encode_json(\%statistics);
                  write_file( $temp_dir . "/" . $fname . 'stats.json', { binmode => ':raw' }, $json );
            }
            undef $CRISPR_hash{$fname};
            $CRISPR_cnt{$fname} = 0;
            $seqcount++;
            push(@fname_array,$fname);
            print "$fname has been searched for designs. Search is done ".($seqcount/$maxcount*100)." %.\n";
      } # end the loop after Hash and statistic were created, so the Bowtie Magic will be executed only once
      %statistics    = ();
      %CRISPR_hash   = ();
      foreach my $fname ( @fname_array ) {
            my $json = read_file( $temp_dir . "/" .$fname. '.json', { binmode => ':raw' } );
            %CRISPR_hash = ( %CRISPR_hash, %{ decode_json $json } );
            unlink $temp_dir . "/" . $fname. ".json";
            $json = read_file( $temp_dir . "/" . $fname. 'stats.json', { binmode => ':raw' } );
            %statistics=( %statistics, %{ decode_json $json });
            unlink $temp_dir . "/" . $fname. "stats.json";
      }          
      ###########################################################################################################################################################################
      #Bowtie for single sequence
      ###########################################################################################################################################################################
                        if ($something{"kind"} eq "single") {
                              open (my $crisprs, ">", $temp_dir . "/temp_CRISPRS.fasta");
                                    foreach my $seq ( sort( keys(%CRISPR_hash) ) ) {
                                          foreach my $key ( sort( keys(%{$CRISPR_hash{$seq}}) ) ) {
                                                foreach my $letter ("A","C","G","T"){
                                                           print $crisprs "\>" . $key . "\n";
                                                            if(${ ${ $CRISPR_hash{$seq} } {$key} }{"strand"} eq "minus"){
                                                                  print $crisprs substr(reverse_comp(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}),$something{"unspecific_leading_bases"},(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-3-$something{"unspecific_leading_bases"})).$letter.substr(reverse_comp(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}),(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-2)). "\n"; #$whole_CRISPR_seq
                                                            }else{
                                                                 print $crisprs substr(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"},$something{"unspecific_leading_bases"},(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-3-$something{"unspecific_leading_bases"})).$letter.substr(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"},(length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-2)). "\n"; #$whole_CRISPR_seq
                                                            }
                                                }
                                          }
                                    }
                              close $crisprs;
                              
                              #####################################################################################################################################################################
                              #teemp_sec
                              if ($something{"bowtie_version"} eq "bowtie2") {
                                    if ($something{"offtargetdb"} eq "gDNA") {
                                          system( 'bowtie2 -p '.$parallel_number.' -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath .".dna". ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }elsif($something{"offtargetdb"} eq "cDNA"){
                                          system( 'bowtie2 -p '.$parallel_number.' -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath.".cdna" . ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }else{
                                          system( 'bowtie2 -p '.$parallel_number.' -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath.".genome" . ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }
                              }else{
                                    if ($something{"offtargetdb"} eq "gDNA") {
                                          system( 'bowtie ' . $databasepath .".dna". ' ' . $temp_dir . "/" . 'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir . '/temp_out.bwt' );
                                    }elsif($something{"offtargetdb"} eq "cDNA"){
                                          system( 'bowtie ' . $databasepath .".cdna". ' ' . $temp_dir . "/" . 'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir . '/temp_out.bwt' );
                                    }else{
                                          if (-e $databasepath.".genome.1.ebwtl") {
                                                system( 'bowtie ' . $databasepath.".genome" . ' ' . $temp_dir . "/" .'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --large-index  --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir .'/temp_out.bwt' );
                                          }else{
                                                system( 'bowtie ' . $databasepath.".genome" . ' ' . $temp_dir . "/" .'temp_CRISPRS.fasta -f -v 3 -y -k 30 -S --sam-nohead --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir .'/temp_out.bwt' );
                                          }
                                    }
                              }
                              open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                                    while (my $line = <$bowtie>) {
                                          chomp $line;
                                          my @line = split( "\t", $line );
                                          if ($line[2] eq "*") {
                                                $line[0] =~m/(\S+)_(\S+)_/;
                                                my $seq = $1;
                                                ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."¤¤0¤¤0¤¤NA¤¤0¤¤NA";
                                          }else{
                                                $line=~m/NM:i:(\d+)/;
                                                my $edit_distance=$1;
                                                if (($edit_distance <= $something{"edit_distance_allowed"})) {
                                                      my @line = split( "\t", $line );
                                                      #decide if it is forward (fw) or backward (bw) query-sequence
                                                      my $direction = "rc";
                                                      if ($line[1] == 0 || $line[1] == 256) {
                                                            $direction = "fw";
                                                      }
                                                      if ( $line[0] =~ m/([^_]+_[^_]+_[^_]+)/ig ) {
                                                            $line[0] =~m/(\S+)_(\S+)_/;
                                                            my $seq = $1;
                                                            my @matchstringo=make_mismatch_string (\$line,$something{"unspecific_leading_bases"}, $direction);
                                                            if ( ($direction eq "fw" && $matchstringo[scalar(@matchstringo)-1] ne "X" && $matchstringo[scalar(@matchstringo)-2] ne "X" && $matchstringo[scalar(@matchstringo)-3] ne "X")
                                                                || ($direction eq "rc" && $matchstringo[0] ne "X" && $matchstringo[1] ne "X" && $matchstringo[2] ne "X")
                                                                ) {
                                                                  my $startcoordinate=0;
                                                                  if ($something{"offtargetdb"} eq "genomicDNA") {
                                                                        my $namestuff="";
                                                                        if ( !exists $trees{$line[2]} ) { #TODO if und else fast identisch, aber relativ kurzer Part und daher Funktion performancetechnisch nachteilig
                                                                              $trees{$line[2]} = build_tree( $something{"databasepath"}  . $something{"ref_organism"} . "/" . $line[2] . "_indexed" );
                                                                        }
                                                                        my $annotations = $trees{$line[2]}->fetch( int($line[3]), int(($line[3])) );
                                                                        foreach  my $anno ( @{$annotations} ) {
                                                                              if ( $anno =~ m/gene_(\S+)_([0-9]+)_([0-9]+)/ ) {
                                                                                    $namestuff=$1;
                                                                                    #$startcoordinate=$2;
                                                                              }
                                                                        }
                                                                       if ($namestuff eq "" && $something{"ignore_intergenic"} eq "true") {                                                                              
                                                                        }elsif($namestuff ne ""){
                                                                              ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$namestuff."¤¤".($line[3]-$startcoordinate)."¤¤".($line[3]+@matchstringo-$startcoordinate)."¤¤".join("",@matchstringo)."¤¤".$edit_distance."¤¤".$direction;
                                                                        }elsif($namestuff eq "" && $something{"ignore_intergenic"} eq "false"){
                                                                              ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."¤¤".($line[3]-$startcoordinate)."¤¤".($line[3]+@matchstringo-$startcoordinate)."¤¤".join("",@matchstringo)."¤¤".$edit_distance."¤¤".$direction;
                                                                        }
                                                                  }else{
                                                                        ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."¤¤".($line[3])."¤¤".($line[3]+@matchstringo)."¤¤".join("",@matchstringo)."¤¤".$edit_distance."¤¤".$direction;
                                                                  }
                                                            }
                                                      }    
                                                }
                                          }
                                    }
                              close $bowtie;
                              unlink $temp_dir . "/temp_CRISPRS.fasta";
                              unlink $temp_dir . "/temp_out.bwt";
                              if ( $something{"sec_off_target"} eq "true" ) { #ckeck if this is wanted
                                    open (my $crisprs, ">", $temp_dir . "/temp_CRISPRS.fasta");
                                          foreach my $seq ( sort( keys(%CRISPR_hash) ) ) {
                                                foreach my $key ( sort( keys(%{$CRISPR_hash{$seq}}) ) ) {
                                                      print $crisprs "\>" . $key . "\n";
                                                      print $crisprs ${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"} . "\n"; #$whole_CRISPR_seq
                                                }
                                          }
                                    close $crisprs;
                                    
                                    ###############################################################################################################################################################
                                    
                                    #do send a bowtie2 job fot the two temporary written fasta files as if they were paired seqencing reads and save the result in a ~out.bwt file
                                    system( 'bowtie2 -p '.$parallel_number.'  -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $temp_dir .'/temp_sec -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                                          my $edit_distance=0;
                                          while (my $line = <$bowtie>) {
                                                chomp $line;
                                                my @line = split( "\t", $line );
                                                if($line=~m/NM:i:(\d+)/){
                                                      $edit_distance=$1;
                                                }
                                                if (($edit_distance <= $something{"edit_distance_allowed"}) && ($line[2] ne "*")) {
                                                      if ( $line[0] =~ m/([^_]+_[^_]+_[^_]+)/ig ) {
                                                            my $key = $1;
                                                            $line[0] =~m/(\S+)_(\S+)_/;
                                                            my $seq = $1;
                                                            push @{ ${ ${ $CRISPR_hash{$seq} } {$key} }{"sec_hits"} }, $line[2];
                                                      }
                                                }
                                          }
                                    close $bowtie;
                                    unlink $temp_dir . "/temp_CRISPRS.fasta";
                                    unlink $temp_dir . "/temp_out.bwt";
                              }
                        }else{
                              
                              #####################################################################################################################################################################
                              #Bowtie for double sequence
                              #####################################################################################################################################################################
                              
                              open (my $leftcrisprs, ">", $temp_dir . "/temp_LEFTCRISPRS.fasta");
                                    open (my $rightcrisprs, ">", $temp_dir . "/temp_RIGHTCRISPRS.fasta");
                                          foreach my $seq ( sort( keys(%CRISPR_hash) ) ) {
                                                foreach my $key ( sort( keys(%{$CRISPR_hash{$seq}}) ) ) {
                                                      my $stringleft = @{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[0];
                                                      my $stringright = @{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[1];
                                                      foreach my $left_PAM ("CCA","CCC","CCT","CCG","CTA","CTC","CTT","CTG"){
                                                            foreach my $right_PAM ("CCA","CCC","CCT","CCG","CTA","CTC","CTT","CTG"){
                                                                  print $leftcrisprs "\>" . $key . "\n";
                                                                  print $leftcrisprs $left_PAM.substr(@{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[0], 3, (length($stringleft) -3 - $something{"unspecific_leading_bases"})) . "\n"; #TODO am Ende entfernen (...,0,leading)
                                                                  print $rightcrisprs "\>" . $key . "\n";
                                                                  print $rightcrisprs $right_PAM.substr(reverse_comp(@{${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}}[1]), 3, (length($stringright) -3 - $something{"unspecific_leading_bases"})) . "\n"; #TODO am Ende entfernen (...,0,leading)
                                                            }
                                                      }
                                                }
                                          }
                                    close ($leftcrisprs);
                              close ($rightcrisprs);
                              
                              #####################################################################################################################################################################
                              #teemp_sec
                              if ($something{"offtargetdb"} eq "gDNA") {
                                    system( 'bowtie2 -p '.$parallel_number.' -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath .".dna". ' -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }elsif($something{"offtargetdb"} eq "cDNA"){
                                    system( 'bowtie2 -p '.$parallel_number.' -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath.".cdna". ' -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }else{
                                    system( 'bowtie2 -p '.$parallel_number.' -f -k 30 --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath.".genome".'  -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }
                              open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                              my $edit_distance=0;
                              my $was_hit=0;
                                    while (my $line = <$bowtie>) {
                                          chomp $line;
                                          my @line = split( "\t", $line );
                                          if ($line[2] ne "*") {
                                                #decide if it is forward (fw) or backward (rc) query-sequence - only for the first pair, the 2nd will be the oposite - Workaround: oposite interpretation of the flag-numbers
                                                my $direction           = "fw";
                                                if ($line[1]==97 || $line[1]==99 || $line[1]==355 || $line[1]==161 || $line[1]==163 || $line[1]==419)  {
                                                      $direction        = "rc";
                                                }
                                                if ( $line[0] =~ m/([^_]+_[^_]+_[^_]+)/ig ) {
                                                      $line[0] =~m/(\S+)_(\S+)_/;
                                                      my $seq = $1;
                                                      my @matchstringo=make_mismatch_string (\$line,$something{"unspecific_leading_bases"}, $direction);
                                                      my $startcoordinate=0;
                                                      my $spacer = abs($line[8]) - ((abs($line[8]) - abs($line[3] - $line[7])) * 2);
                                                      if ($something{"offtargetdb"} eq "genomicDNA") {                                                      
                                                            if ($line[1]==81 || $line[1]==83 || $line[1]==97 || $line[1]==99 || $line[1]==339 || $line[1]==355) {
                                                                  $line=~m/NM:i:(\d+)/;
                                                                  $edit_distance=$1;
                                                                  if (($edit_distance <= $something{"edit_distance_allowed"})) {
                                                                        if ( (($direction eq "fw" && $matchstringo[scalar(@matchstringo)-1] ne "X" && $matchstringo[scalar(@matchstringo)-2] ne "X" && $matchstringo[scalar(@matchstringo)-3] ne "X")
                                                                              || ($direction eq "rc" && $matchstringo[0] ne "X" && $matchstringo[1] ne "X" && $matchstringo[2] ne "X"))
                                                                           ) {
                                                                              my $namestuff="";
                                                                              if ( !exists $trees{$line[2]} ) { #TODO if und else fast identisch, aber relativ kurzer Part und daher Funktion performancetechnisch nachteilig
                                                                                    $trees{$line[2]} = build_tree( $something{"databasepath"} . $something{"ref_organism"} . "/" . $line[2] . "_indexed" );
                                                                              }
                                                                              my $annotations = $trees{$line[2]}->fetch( int($line[3]), int(($line[3])) );
                                                                              foreach  my $anno ( @{$annotations} ) {
                                                                                    if ( $anno =~ m/gene_(\S+)_([0-9]+)_([0-9]+)/ ) {
                                                                                          $namestuff=$1;
                                                                                          #$startcoordinate=$2;
                                                                                    } 
                                                                              }
																			  if ($namestuff eq "" && $something{"ignore_intergenic"}=="true") {                                                                              
																				}elsif($namestuff ne ""){
																					  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$namestuff."¤¤".($line[3]-$startcoordinate)."¤¤".($line[3]+@matchstringo-$startcoordinate)."¤¤".join("",@matchstringo)."¤¤".$edit_distance."¤¤".$direction."¤¤".$spacer;
																				}elsif($namestuff eq "" && $something{"ignore_intergenic"}=="false"){
																					  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."¤¤".($line[3]-$startcoordinate)."¤¤".($line[3]+@matchstringo-$startcoordinate)."¤¤".join("",@matchstringo)."¤¤".$edit_distance."¤¤".$direction."¤¤".$spacer;
																				}
                                                                              $was_hit=1;
                                                                        }
                                                                  }
                                                            }elsif($line[1]==161 || $line[1]==163 || $line[1]==145 || $line[1]==147 || $line[1]==403 || $line[1]==419){
                                                                  if ($was_hit==1) {                                                                  
                                                                        $line=~m/NM:i:(\d+)/;
                                                                        $edit_distance=$1;
                                                                        if ( (($direction eq "fw" && $matchstringo[scalar(@matchstringo)-1] ne "X" && $matchstringo[scalar(@matchstringo)-2] ne "X" && $matchstringo[scalar(@matchstringo)-3] ne "X")
                                                                              || ($direction eq "rc" && $matchstringo[0] ne "X" && $matchstringo[1] ne "X" && $matchstringo[2] ne "X"))
                                                                              && ($edit_distance <= $something{"edit_distance_allowed"})
                                                                        ) {
                                                                                    my @hitarray=split(";;",${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"});
                                                                                    my @lasthitarray=split("¤¤",$hitarray[-1]);
                                                                                    $lasthitarray[3].="-".join("",@matchstringo);
                                                                                    $hitarray[-1]=join("¤¤",@lasthitarray);
                                                                                    ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}=join(";;",@hitarray);
                                                                        }else{
                                                                                    my @hitarray=split(";;",${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"});
                                                                                    pop @hitarray;
                                                                                    ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}=join(";;",@hitarray);
                                                                        }
                                                                        $was_hit=0;
                                                                  }
                                                            }                                                  
                                                      }else{
                                                            if ($line[1]<147) {
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."¤¤".($line[3]-$startcoordinate)."¤¤".($line[3]+@matchstringo-$startcoordinate)."¤¤".join("",@matchstringo)."¤¤".$edit_distance."¤¤".$direction."¤¤".$spacer;
                                                            }elsif($line[1]==147){
                                                                  my @hitarray=split(";;",${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"});
                                                                  my @lasthitarray=split("¤¤",$hitarray[-1]);
                                                                  $lasthitarray[3].="-".join("",@matchstringo);
                                                                  $hitarray[-1]=join("¤¤",@lasthitarray);
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}=join(";;",@hitarray);
                                                            }else{
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."¤¤".($line[3]-$startcoordinate)."¤¤".($line[3]+@matchstringo-$startcoordinate)."¤¤".join("",@matchstringo)."¤¤".$edit_distance."¤¤".$direction."¤¤".$spacer;
                                                            }
                                                            $was_hit=1;
                                                      }
                                                }                                                     
                                          }
                                    }
                              close $bowtie;
                              unlink $temp_dir . "/temp_LEFTCRISPRS.fasta";
                              unlink $temp_dir . "/temp_RIGHTCRISPRS.fasta";
                              unlink $temp_dir . "/temp_out.bwt";
                              #####################################################################################################################################################################
                              #  if ( exists $something{"sec_off_target"} ) needed for double! @Florian
                              #####################################################################################################################################################################
                        }
                  
                  
                  ###########################################################################################################################################################################
                  #Evaluate the infos of the CRISPR-Hash (additional information, changes)
                  ###########################################################################################################################################################################
                  
                  my $number_of_hits = 1;
                  SEQLOOP: foreach my $fname ( keys(%CRISPR_hash) ) {
                        CRISPRHASHLOOP: foreach my $key ( keys(%{$CRISPR_hash{$fname}}) ) {
                              if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} ) {
                                    $number_of_hits = scalar(split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"}));
                                    ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"} = $number_of_hits-1;
                                    if (( $number_of_hits > $something{"off-targets-allowed"} + 2 || $number_of_hits < 1 ) ) {
                                          delete $CRISPR_hash{$fname}{$key};
                                          $statistics{$fname}{"Number of designs excluded because they hit multiple targets or none"}++;
                                          next CRISPRHASHLOOP;
                                    }else{
                                          $statistics{$fname}{"Number of designs that hit a specific target"}++;
                                    }
                              }else{
                                    delete $CRISPR_hash{$fname}{$key};
                                    next CRISPRHASHLOOP;
                              }
                              
                              
                              if ($something{"sec_off_target"} eq "true" ){
                                    if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"sec_hits"} ) {
                                          $number_of_hits = @{ ${ ${ $CRISPR_hash{$fname} } {$key} }{"sec_hits"} };
                                          if (( $number_of_hits > 0 ) ) {
                                                delete $CRISPR_hash{$fname}{$key};
                                                $statistics{$fname}{"Number of designs excluded because they hit a secondary offtarget"}++;
                                                next CRISPRHASHLOOP;
                                          }else{
                                                $statistics{$fname}{"Number of designs that do not hit a secondary target"}++;
                                          }
                                    }else{
                                          $statistics{$fname}{"Number of designs that do not hit a secondary target"}++;
                                    }
                              }
                              
                              my $whole_crisp_seq = "";
                              if($something{"kind"} ne "single"){
                                    $whole_crisp_seq = join( "", @{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}});
                              }else{
                                    $whole_crisp_seq = ${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"};
                              }
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"} = 120;
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"exon"} = join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} } {"exon"} } ) );
                        }
                        
                        my %exon_count = ();
                        EXONLOOP: foreach my $key (keys(%{$CRISPR_hash{$fname}}) ) {
                              my $exon = ${ ${ $CRISPR_hash{$fname} } {$key} }{"exon"};
                              $exon_count{$exon}++;
                              if ( $exon_count{$exon} > $something{"max_per_exon"} ) {
                                    delete $CRISPR_hash{$fname}{$key};
                                    $statistics{$fname}{"Number of designs excluded because the maximum of designs per exon was exceeded"}++;
                                    next EXONLOOP;
                              }
                        }
                  }
                  
                  #################################################################################################################################################################################
                  #reopen the main loop which is looping through the sequence found in the SeqIO object, may be one or many more
                  #################################################################################################################################################################################
                  if ( $something{"data_type"} eq "ensemble_acc" ) {
                        $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ); #read the temporary fasta file            
                  } else { 
                        $seqio_obj = Bio::SeqIO->new( -file => $something{"input_file"}, -format => "fasta" ) or die "Your input wasn't FASTA format \n e.g.: >some Sequence<br>AGCTGATCGATCTAGCTAGCTGCTAGCTAGTCGATCGATGCTAGCTAGCTAGCTGCTAG\n"; #if neither online or pasted sequences are used it will use the input file as sequence input
                  }
                  @seq_array=();
                  while( my $seq = $seqio_obj->next_seq() ) {
                        my $seq_obj=$seq;
                        #reinitiate values...
                        my $chrom               = ""; #create an empty name chromosome
                        my $location_offset     = 0;
                        my $location_end        = 0; #set the location offset value to 0
                        if ( ( $seq_obj->description ) =~ m/chrom:([\w\.]+):(\d+)..(\d+)/ig) { #if the descrition of the sequence contains information of the form loc=chr1:0..13000000
                              #print the result of the pattern matching
                              $chrom = $1; #save the location information in the certain objects
                              $location_offset = $2;
                              $location_end = $3; #the location offset is the start of the sequence on that certain chromosome, for the whole chromosome this is of cause 0
                        } else {
                              $dont_asses_context = 1;
                        }
                        #define the name for the current sequence as its display id
                        $fname = ( $seq_obj->display_id ); #current name is the sequence' id
                        if ( $chrom eq "" ) { $chrom = $fname; } #if $chrom is still empty fill it with the sequence' id
                        my $whole_seq = $seq_obj->seq; #deduce the complete nucleotide sequence as a alphanumerical string from the SeqIo object
                        
                        ###########################################################################################################################################################################
                        #Build te image of the gene
                        ###########################################################################################################################################################################
                        my %transcripts_hash=();
                        my %CDS_hash=();
                        if ($CRISPR_cnt{$fname} == 0) { #only execute following calucaltions if it is not already done for the sequence
                              if ( exists $trees{$chrom} ) {
                                    my $annotations   = $trees{$chrom}->fetch( int($location_offset), int($location_end) );
                                    foreach my $anno ( sort( @{$annotations} ) ) {
                                          if ( $anno =~ m/exon::(\S+)::(.+)_(\d+)_(\d+)$/ig ) {
                                                my @pair = ( $3, $4 );
                                                ${$transcripts_hash{$1}}{$2}=\@pair;
                                          }
                                          if ( $anno =~ m/CDS::(\S+)/ig ) {
                                                $CDS_hash{$1}++;
                                          }
                                    }
                              }
                              my %strand=();
                              foreach my $key ( sort( keys(%transcripts_hash) ) ) {
                                          my $curstart = 0;
                                          my $diff = 0;
                                          $strand{$key} = 1;
                                          foreach my $pair_ref (sort(keys %{$transcripts_hash{$key}})) {
                                                $diff = ${$transcripts_hash{$key}}{$pair_ref}->[0]-$curstart;
                                                $curstart = ${$transcripts_hash{$key}}{$pair_ref}->[0];
                                          }
                                          if ($diff<0) {
                                                $strand{$key} = -1;
                                          }
                              }
                        }
                        foreach my $key ( keys(%{$CRISPR_hash{$fname}}) ) {
                              my @targets = split(";;", ${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} );
                              my $number_off_targets=1;
                              foreach my $hit (@targets){
                                    if ($hit ne "") {
                                          #print $hit."\n";
                                          my @splithit=split("¤¤",$hit);
                                          if (${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}-((20-(100/${ ${ $CRISPR_hash{$fname} } {$key} }{"length"}*$splithit[4]))/$number_off_targets)>0) {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}-((20-(100/${ ${ $CRISPR_hash{$fname} } {$key} }{"length"}*$splithit[4]))/$number_off_targets);
                                                $number_off_targets++;
                                          }else{
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=0;
                                          }
                                    }
                              }
                              if(exists(${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"})){
                                    @{${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"}}[0]=${ ${ $CRISPR_hash{$fname} } {$key} }{"score"};
                                    ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"};
                                    @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[1]=@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[1]*100/((5*(scalar(keys(%CDS_hash))))+(scalar(keys(%CDS_hash)))+(5*(scalar(keys(%transcripts_hash))))+1);
                                    #@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[2]=100*(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[2]-(-21))/(40-(-21));
                              }
                              foreach my $i (0..2){
                                    @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[$i]=round_digits(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[$i],4);
                              }
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"spec_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[0];
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"anno_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[1];
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"eff_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[2];
                        }
                       
                        #####################################################################################################################################################################
                        # Create the Table for the findings
                        #####################################################################################################################################################################
                              
                              open (my $outfiletab, ">", $temp_dir . "/" . $fname . "_" . "table.tab");
                                    if ($something{"kind"} eq "single") {
                                          print $outfiletab "Name\tLength\tStart\tEnd\tStrand\tNucleotide sequence\tGene Name\tTranscripts\tTranscript:: Exon\tNumber of Cpg Islands hit\tSequence around the cutside\t%A %C %T %G\tS-Score\tA-Score\tE-Score\tEXTRA_Score\tpercent of total transcripts hit\tTarget\tMatch-start\tMatch-end\tMatchstring\tEditdistance\tNumber of Hits\tDirection\n";
                                    }
                                    else {
                                          print $outfiletab "Name\tLength\tStart\tEnd\tStrand\tNucleotide sequence\tGene Name\tTranscripts\tTranscript:: Exon\tNumber of Cpg Islands hit\tSequence around the cutside\t%A %C %T %G\tS-Score\tA-Score\tE-Score\tEXTRA_Score\tpercent of total transcripts hit\tTarget\tMatch-start\tMatch-end\tMatchstring\tEditdistance\tNumber of Hits\tDirection\tSpacer\n";
                                    }
									if ($something{"purpose"} eq "non-coding") {
											PRINTLOOP: foreach my $key ( sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} cmp $CRISPR_hash{$fname}{$a}->{"spec_score"} } keys(%{$CRISPR_hash{$fname}}) ) {
											  $statistics{$fname}{"Number of successful designs"}++;
											  my @targets=split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} );
											  #write the tab-delimited file
											  HITLOOP: foreach my $hit (@targets){
													if ($hit ne "") {
														  #print the candidates name
														  print $outfiletab "$key\t";
														  my @splithit = split("¤¤",$hit);
														  #print its length on this whole sequence these are not genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} . "\t";
														  }
														  #print its start on this whole sequence these are not genomic coordinates
														   my @locus=split("::",$statistics{$fname}{"seq_location"});
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1] . "\t";
														  }
														  #print its end on this whole sequence these are  genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1] . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"} ) {
																if ($something{"kind"} eq "single") {
																	  if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
																			print $outfiletab reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,2)."N ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},3)) . "\t";
																	  }else{
																			print $outfiletab substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-3)." N".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-2) . "\t";
																	  }
																}else{
																	  print $outfiletab reverse_comp(
																									 substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],0,2)."N ".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],3))
																	  ."_" .
																								   substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],0,length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-3)." N".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-2). "\t";
																}
														  }
														  #print the gene name it overlaped with if there was any
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"gene"} ) {
																print $outfiletab join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"gene"} } ) ) . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  my $percent_of_transcripts_hit = "NA";
														  my $number_of_transcripts = scalar(keys(%transcripts_hash));
														  my $number_of_target_transcripts = scalar(keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"}}));
														  if ($number_of_transcripts != 0) {
																$percent_of_transcripts_hit = $number_of_target_transcripts*100/$number_of_transcripts;
														  }
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ) {
																print $outfiletab ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"}."\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  #print the exons it overlaped with if there was any
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"} ) {
																print $outfiletab join( " ", keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"}}))  . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  #print the number of CpG islands it overlapped with if there was any
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} ) {
																print $outfiletab ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  #write the CRISPR as annotation to the original sequence
														  my $whole_crisp_seq = substr( $whole_seq, ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}, ( ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} + 2 ) );
														  #print homology matrix also sub divided into left and right sequence
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} ) {
																print $outfiletab "left::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} . "::right::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"right"} . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
														  
														  print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"}, "\t", join("\t",@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}});
														  print $outfiletab "\t$percent_of_transcripts_hit\t";
														  print $outfiletab $splithit[0]."\t".$splithit[1]."\t".$splithit[2]."\t".$splithit[3]."\t".$splithit[4]."\t";
														  print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"}, "\t";
														  print $outfiletab $splithit[5]."\t";
														  if (!($something{"kind"} eq "single")) {
																print $outfiletab $splithit[6]."\t";
														  }
														  
														  
														  #print the end of the line as a new line
														  print $outfiletab "\n";
														  #make a featureanntotaion for that CRISPR
													}
											  }
											  #print the candidates name
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} }{"gene"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} = join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} } {"gene"} } ) );
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} = "NA";
											  }
											  my $percent_of_transcripts_hit="NA";
											  my $number_of_transcripts = scalar(keys(%transcripts_hash));
											  my $number_of_target_transcripts = scalar(keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} }{"transexons"}}));
											  if ($number_of_transcripts != 0) {
													$percent_of_transcripts_hit=$number_of_target_transcripts*100/$number_of_transcripts;
											  }
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transcripts"} = ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ;
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transcripts"} = "NA";
											  }
											  #print the exons it overlaped with if there was any
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transexons"} = join( " ", keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} } {"transexons"} }));
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transexon"} = "NA";
											  }
											  #print the number of CpG islands it overlapped with if there was any
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"CpG"} = ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"};
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"CpG"} = "NA";
											  }
											  #write the CRISPR as annotation to the original sequence
											  my $whole_crisp_seq = substr( $whole_seq, ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}, ( ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} + 2 ) );
											  #print homology matrix also sub divided into left and right sequence
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} = "left::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} . "::right::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"right"};
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} = "NA";
											  }
											  undef ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"}; #delete the "context" hash of the crispr hash it is not longer needed
											  ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
										}
									}else{									
										PRINTLOOP: foreach my $key ( sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} } sort { $CRISPR_hash{$fname}{$b}->{"anno_score"} <=> $CRISPR_hash{$fname}{$a}->{"anno_score"} } sort { $CRISPR_hash{$fname}{$b}->{"eff_score"} <=> $CRISPR_hash{$fname}{$a}->{"eff_score"} } sort { $CRISPR_hash{$fname}{$b}->{"exon"} cmp $CRISPR_hash{$fname}{$a}->{"exon"} } keys(%{$CRISPR_hash{$fname}}) ) {
											  $statistics{$fname}{"Number of successful designs"}++;
											  my @targets=split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} );
											  #write the tab-delimited file
											  HITLOOP: foreach my $hit (@targets){
													if ($hit ne "") {
														  #print the candidates name
														  print $outfiletab "$key\t";
														  my @splithit = split("¤¤",$hit);
														  #print its length on this whole sequence these are genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} . "\t";
														  }
														  #print its start on this whole sequence these are  genomic coordinates
														  my @locus=split("::",$statistics{$fname}{"seq_location"});
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1] . "\t";
														  }
														  #print its end on this whole sequence these are  genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1] . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"} ) {
																if ($something{"kind"} eq "single") {
																	  if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
																			print $outfiletab reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,2)."N ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},3)) . "\t";
																	  }else{
																			print $outfiletab substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-3)." N".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-2) . "\t";
																	  }
																}else{
																	  print $outfiletab reverse_comp(
																									 substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],0,2)."N ".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],3))
																	  ."_" .
																								   substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],0,length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-3)." N".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-2). "\t";
																}
														  }
														  #print the gene name it overlaped with if there was any
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"gene"} ) {
																print $outfiletab join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"gene"} } ) ) . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  my $percent_of_transcripts_hit = "NA";
														  my $number_of_transcripts = scalar(keys(%transcripts_hash));
														  my $number_of_target_transcripts = scalar(keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"}}));
														  if ($number_of_transcripts != 0) {
																$percent_of_transcripts_hit = $number_of_target_transcripts*100/$number_of_transcripts;
														  }
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ) {
																print $outfiletab ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"}."\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  #print the exons it overlaped with if there was any
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"} ) {
																print $outfiletab join( " ", keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"}}))  . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  #print the number of CpG islands it overlapped with if there was any
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} ) {
																print $outfiletab ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  #write the CRISPR as annotation to the original sequence
														  my $whole_crisp_seq = substr( $whole_seq, ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}, ( ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} + 2 ) );
														  #print homology matrix also sub divided into left and right sequence
														  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} ) {
																print $outfiletab "left::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} . "::right::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"right"} . "\t";
														  } else {
																print $outfiletab "NA\t";
														  }
														  ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
														  
														  print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"}, "\t", join("\t",@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}});
														  print $outfiletab "\t$percent_of_transcripts_hit\t";
														  print $outfiletab $splithit[0]."\t".$splithit[1]."\t".$splithit[2]."\t".$splithit[3]."\t".$splithit[4]."\t";
														  print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"}, "\t";
														  print $outfiletab $splithit[5]."\t";
														  if (!($something{"kind"} eq "single")) {
																print $outfiletab $splithit[6]."\t";
														  }
														  
														  
														  #print the end of the line as a new line
														  print $outfiletab "\n";
														  #make a featureanntotaion for that CRISPR
													}
											  }
											  #print the candidates name
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} }{"gene"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} = join( "_", keys( %{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} } {"gene"} } ) );
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"gene"} = "NA";
											  }
											  my $percent_of_transcripts_hit="NA";
											  my $number_of_transcripts = scalar(keys(%transcripts_hash));
											  my $number_of_target_transcripts = scalar(keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} }{"context"} }{"transexons"}}));
											  if ($number_of_transcripts != 0) {
													$percent_of_transcripts_hit=$number_of_target_transcripts*100/$number_of_transcripts;
											  }
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transcripts"} = ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transcripts"} ;
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transcripts"} = "NA";
											  }
											  #print the exons it overlaped with if there was any
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"transexons"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transexons"} = join( " ", keys(%{ ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} } {"transexons"} }));
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"transexon"} = "NA";
											  }
											  #print the number of CpG islands it overlapped with if there was any
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"CpG"} = ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"CpG"};
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} }{"CpG"} = "NA";
											  }
											  #write the CRISPR as annotation to the original sequence
											  my $whole_crisp_seq = substr( $whole_seq, ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}, ( ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} + 2 ) );
											  #print homology matrix also sub divided into left and right sequence
											  if ( exists ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} ) {
													${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} = "left::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"left"} . "::right::" . ${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} }{"right"};
											  } else {
													${ ${ $CRISPR_hash{$fname} } {$key} } {"homology"} = "NA";
											  }
											  undef ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"}; #delete the "context" hash of the crispr hash it is not longer needed
											  ${ ${ $CRISPR_hash{$fname} } {$key} }{"Nuc_comp"} = join( " ", my @basecomp = find_base_count( $whole_crisp_seq ) ); #store the nucleotide composition as a string object in the hash
										} #end of printloop
									}
                              close $outfiletab;
                              #####################################################################################################################################################################
                              #If asked, write all the stuff to a gff format file
                              #####################################################################################################################################################################
                              
                              if ($something{"out_gff"} eq "true") { 
                                    open(my $gfffile, ">",$temp_dir . "/" . $fname . ".gff" ) or die $!;
                                     print $gfffile "##gff-version 3\n";
                                          PRINTLOOP: foreach my $key ( sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} } sort { $CRISPR_hash{$fname}{$b}->{"anno_score"} <=> $CRISPR_hash{$fname}{$a}->{"anno_score"} } sort { $CRISPR_hash{$fname}{$b}->{"eff_score"} <=> $CRISPR_hash{$fname}{$a}->{"eff_score"} } sort { $CRISPR_hash{$fname}{$b}->{"exon"} cmp $CRISPR_hash{$fname}{$a}->{"exon"} } keys(%{$CRISPR_hash{$fname}}) ) {
                                                my @locus=split("::",$statistics{$fname}{"seq_location"});
                                                print $gfffile $locus[0]."\tcld\tCRISPRtarget\t".(${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1])."\t";
                                                print $gfffile (${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1])."\t".sum(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}})."\t";
                                                if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                      print $gfffile "-"."\t."."\t";
                                                }else{
                                                      print $gfffile "+"."\t."."\t";
                                                }
                                                print $gfffile "id=".$key."; ";
                                                print $gfffile "spec_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"spec_score"}."; ";
                                                print $gfffile "anno_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"anno_score"}."; ";
                                                print $gfffile "eff_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"eff_score"}."; ";
                                                if ($something{"kind"} eq "single") {
                                                      if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                  print $gfffile "seq=".reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,2)."N_".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},3)) . "; ";
                                                      }else{
                                                                  print $gfffile "seq=".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-3)."_N".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-2) . "; ";
                                                      }
                                                }else{
                                                      print $gfffile "seq=".reverse_comp(
                                                            substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],0,2)."N_".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[0],3))
                                                                  ."_" .
                                                            substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],0,length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-3)."_N".substr(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1],length(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"}}[1])-2). "; ";
                                                }
                                                print $gfffile "offtargetcount=".${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"}.";\n";
                                          }
                                          close $gfffile;
                              }
                              
                            #####################################################################################################################################################################
                              # ZIP the stuff
                              #####################################################################################################################################################################
                              
                              unlink $temp_dir . "/tempfile.fasta";
                              my $zip    = Archive::Zip->new();
                              my $member = "";
                              if ( $something{"out_gff"} eq "true") { $member = $zip->addFile( $temp_dir . "/" . $fname . ".gff", $fname . "_CRISPR.gff" ); }
                              $member = $zip->addFile( $temp_dir . "/" . $fname . "_" . "table.tab", $fname . "_CRISPR.tab" );
                              $zip->writeToFileNamed( $temp_dir . '/' . $fname . '.zip' );
                              
                              #####################################################################################################################################################################
                              # Print the report site with the results (header was created earlier in the loop)
                              #####################################################################################################################################################################
                            if ($something{"draw_html_report"} eq "true" ){
                              open( my $report, ">>", "$temp_dir/index.html" ) or die "can not open report file";
                                    #print the statistics
                                    print $report '   <tr><td>
                                                            <span class="main3"><B>Query name: '.$statistics{$fname}{"seq_name"}.'   Query length: '.$statistics{$fname}{"seq_length"}.'   Query location: '.$statistics{$fname}{"seq_location"}.'</B><br><br>';
                                    foreach my $statistic(sort {$b cmp $a} keys(%{$statistics{$fname}})){
                                          if ($statistic =~ m/Number/ig) {
                                                print $report $statistic.' = '.$statistics{$fname}{$statistic}.'<br>';
                                          }
                                    }
                                    print $report '         </span><br><br>
                                                      </td></tr>';
                                    #print the table
                                    print $report '   <tr><td>
                                                            <table class="talenHits" style="text-align: center; width: 800px; table-layout: fixed; border-bottom: 2px solid black;border-top: 1px solid black; padding: 2px; background-color: white;">';
                                    open (my $tabfile, "<", $temp_dir . "/" . $fname . "_" . "table.tab");
                                         my $count = 0;
                                                while (my $line = <$tabfile>){
                                                      chomp $line;
                                                      my @line=split("\t",$line);
                                                      print $report '<tr class="entry" id="'.$line[0].'">';
                                                      foreach my $element (0,5,12,16,19,21){
                                                            if ($count == 0){
                                                                  if ($element==12) {
                                                                       print $report '<th style="text-align: center; border-bottom: 1px solid black; color: white; padding: 2px; background-color: #2662C3;"> SAE-Score </th>';
                                                                  }else{
                                                                        print $report '<th style="text-align: center; border-bottom: 1px solid black; color: white; padding: 2px; background-color: #2662C3;"> '.$line[$element].' </th>';
                                                                  }
                                                            }else{
                                                                  my $hittype = "bad";
                                                                  if ($line[21] == 1) {
                                                                        $hittype = "good";
                                                                  }
                                                                  if ($element == 19) {
                                                                        #create popup for the matchstring-info if the option is chosen, otherwise only the matchstring itself
                                                                        if ($something{"match_info"} eq "true") {
                                                                              my $popup = "";
                                                                              if ($something{"kind"} eq "single") {
                                                                                    $popup = create_popup(\@line, $databasepath, $something{"unspecific_leading_bases"}, 10);
                                                                              }else {
                                                                                    $popup = create_popupD(\@line, $databasepath, $something{"unspecific_leading_bases"}, 5);
                                                                              }
                                                                              $line[0] =~ s/\W/_/ig; #for html
                                                                              $line[14] =~ s/\W/_/ig; #for html
                                                                              print $report '<td class="'.$hittype.' main-button" style="text-align: center; word-wrap: break-word; font-family: arial, monospace; color: black; font-size:12px; background-color: white;">'
                                                                                                .'<div id="dialog_'.$line[0].'__'.$line[16].'__'.$line[19].'" class="dialog" title="Matchstring Info for '.$line[0].' on Target '.$line[14].'">'.$popup.'</div>';
                                                                              print $report     '<button id="'.$line[0].'__'.$line[16].'__'.$line[19].'" class="opener">Matchstring Info</button>'
                                                                                          .'</td>';
                                                                        }else {
                                                                              print $report '   <td class="'.$hittype.' main" style="text-align: center; word-wrap: break-word; font-family: arial, monospace; color: black; font-size:12px; background-color: white;">'
                                                                                                      . print_offtarget_string($line[$element]) .
                                                                                                '</td>';
                                                                        }
                                                                  }elsif( $element == 12 ){
                                                                        print $report '<td class="'.$hittype.'" >';
                                                                        print $report hor_bar_chart($line[12],$line[13],$line[14]);
                                                                        print $report '</td>';
                                                                  }else{
                                                                        print $report '   <td class="'.$hittype.' main" style="text-align: center; word-wrap: break-word; font-family: arial, monospace; color: black; font-size:12px; background-color: white;"> '
                                                                                                .$line[$element].
                                                                                          '</td>';
                                                                  }
                                                            }
                                                      }
                                                      $count++;
                                                      print $report '</tr>';
                                                }                                          
                                    close $tabfile;
                                    print $report '         </table>
                                                      </td></tr>
                                                      <tr><td>
                                                            <br><br>
                                                      </td></tr>';
                                    #print the Image, the ZIP-Download and the bottom-line
                                    print $report '   <table colspan="2" rowspan="2" style="vertical-align: top;border: 1px dashed black;">
                                                            <tr><td><br><br>';
                                    print $report '               <br><br>
                                                            </td></tr>
                                                            <tr><td>
                                                                  <a href="' . $fname . '.zip' . '" target="_blank" download="'.$fname.'_CRISPRS.zip" class="main"><input type="button" value="Download  results as zipped folder"></a><br><br><br>
                                                            </td></tr>
                                                      </table>';
                                    print $report '   <tr><td>
                                                            <br><br>
                                                      </td></tr>
                                                      <tr><td colspan="2" rowspan="1" style="vertical-align: top;">
                                                            <embed style="width: 100%; height: 10px;" alt="there should appear a  line" src="/cld/gradientline.svg" type="image/svg+xml"><br>
                                                      </td></tr>
                                                      <tr><td>
                                                            <br><br>
                                                      </td></tr>';
                              close $report;                              
                              $CRISPR_cnt{$fname}++;
                        }
                  print "$fname is completed 100%\n";
                  } #end Sequence loop
                  
                  if ($something{"draw_html_report"} eq "true" ){
                	  open( my $report, ">>", "$temp_dir/index.html" ) or die "can not open report file";
                        	
									print $report ' </tr>
													</table>
													</div>
													</body>
													</html>';
					close $report;
                  }
                  opendir(TEMPDIR,$temp_dir);
                        open (my $outfile, ">", $temp_dir."/all_results_together.tab");
                              open (my $outgff, ">", $temp_dir."/all_results_together.gff");
                              my $head_line=0;
                              my $head_line_gff=0;
                              foreach my $filename (readdir(TEMPDIR)){
                                    if($filename=~m/table\.tab/){
                                          open (my $file, "<", $temp_dir."/".$filename);
                                                while (<$file>){
                                                      if (($_=~m/^Name/) && ($head_line==0)) {
                                                            print $outfile $_;
                                                            $head_line++;
                                                      }elsif(($_=~m/^Name/) ){
                                                      }else{                                                      
                                                            print $outfile $_;
                                                      }
                                                }
                                          close $file;
                                          unlink $temp_dir."/".$filename;
                                    }elsif(($filename=~m/\.gff/) && !($filename=~m/all_results_together\.gff/)){
                                          open (my $file, "<", $temp_dir."/".$filename);
                                                while (<$file>){
                                                      if (($_=~m/^.*gff.*/) && ($head_line_gff==0)) {
                                                            print $outgff $_;
                                                            $head_line_gff++;
                                                      }elsif(($_=~m/^.*gff.*/) ){
                                                      }else{                                                      
                                                            print $outgff $_;
                                                      }
                                                }
                                          close $file;
                                          unlink $temp_dir."/".$filename;
                                    }
                              }
                        close $outgff;     
                  close $outfile;
                  closedir TEMPDIR;
                  return $temp_dir;
}



=head2 find_and_print_CRISPRS
#########################################################################################
#name:      find_and_print_CRISPRS
#function:  find the CRISPRS and build a hash for CRISPRS and the statistics...
#input:     (chrom, location_offset, builded tree as reference, dont_asses_context,
#           global temp_dir, parallel_number, something-Hashreference)
#output:    hash for CRISPRS and hash for statistics
#########################################################################################
=cut

sub find_and_print_CRISPRS {
      my $seq_obj_ref               = $_[0];
      my $chrom                     = $_[1];
      my $location_offset           = $_[2];
      my %trees                     = %{ $_[3] };
      my $dont_asses_context        = $_[4];
      my $temp_dir                  = $_[5];
      my $parallel_number           = $_[6];
      my %something                 = %{ $_[7] };
      my $seq_obj                   = $$seq_obj_ref;
	  my $gene_id					= $seq_obj->id;
      my $whole_seq                 = $seq_obj->seq;
      my $count                     = 0;
      my %finished_CRISPR_hash      = ();
      my $pm                        = Parallel::ForkManager->new($parallel_number);
      my $cutnumber                 = int( length($whole_seq) / int( $parallel_number - 1 ) );
      my $cut                       = 0;
      my @cuts                      = ();
      my %tempstatistics            = ();
      my $start_of_start            = 0 ;
      if ($dont_asses_context==0) {
            my $annotations   = $trees{$chrom}->fetch( int($_[2]), int(int($_[8])) );
            foreach my $anno ( sort( @{$annotations} ) ) {
                  if ( $anno =~ m/start_codon::(\S+)::(.+)_(\d+)_(\d+)$/ig ) {
                        if( $start_of_start==0){
                              $start_of_start=$3-$location_offset;
                        }else{
                              $start_of_start=($start_of_start+($3-$location_offset))/2;
                        }
                  }
            }
      }    
            
      while ( $cut <= length($whole_seq) ) {
            push @cuts, $cut;
            $cut = $cut + $cutnumber;
      }
      #################################################################################################################################################################################
      # cut the sequence into equal peaces, so that each forked child can work on one part (paralell!)
      #################################################################################################################################################################################
      
      foreach $cut (@cuts) {
            $pm->start and next;
            my $seq = substr( $whole_seq, $cut, $cutnumber );
            my %CRISPR_hash = ();
            my %Gpos = make_pos_index( \$seq, "G" );
            my %Cpos = make_pos_index( \$seq, "C" );
            my %Apos = make_pos_index( \$seq, "A" );
            my %Tpos = make_pos_index( \$seq, "T" );
			my %combined;
			% {$combined{"G"}}=%Gpos;
			% {$combined{"A"}}=%Apos;
			% {$combined{"C"}}=%Cpos;
			% {$combined{"T"}}=%Tpos;
			my %dont_care_ind=();
			my %dont_care_ind_right=();
			my %PAMindex = ();
            
            ###########################################################################################################################################################################
            # Single Sequence
            ###########################################################################################################################################################################
            
            if ($something{"kind"} eq "single") {
                  
                  #####################################################################################################################################################################
                  # Foward Sequence Calculations
                  #####################################################################################################################################################################
                  
                  if ($something{"preceding"} eq "A") {
                        %dont_care_ind=%Apos;
                  }elsif($something{"preceding"} eq "G"){
                        %dont_care_ind=%Gpos;
                  }elsif($something{"preceding"} eq "C"){
                        %dont_care_ind=%Cpos;
                  }elsif($something{"preceding"} eq "T"){
                        %dont_care_ind=%Tpos;
                  }else{
                        %dont_care_ind=(%Gpos,%Tpos,%Cpos,%Apos);
                  }
                  if ($something{"PAM"} eq "NAG") {
                        %PAMindex=%Apos;
                  } elsif ($something{"PAM"} eq "NGG") {
                        %PAMindex=%Gpos;
                  } else{
                        %PAMindex=(%Apos,%Gpos);
                  }
                  POSLOOP: foreach my $Gposind ( sort( keys(%dont_care_ind) ) ) {
                        LENGTHLOOP: foreach my $length ( ($something{"min_length"}+1) .. ($something{"max_length"}+1) ) {
                              if ( exists $PAMindex{ ( $Gposind + $length ) } && exists $Gpos{ ( $Gposind + $length + 1 ) } ) {
                                    my $taleseq = substr( $seq, $Gposind, $length + 2 );
                                    my @flank_array = find_base_count( $taleseq );
                                    $tempstatistics{"Total number of possible designs"}++;
                                    if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                          $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                          $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                          $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                          !($taleseq=~m/TTTTT/) 
                                    ) {
                                         my $name = ($seq_obj->display_id)."_" . $count . "_" . $cut. "." .(int(abs($Gposind + $cut-$start_of_start)/3));										
										my @new_score=(0,0,0,0);
										if (defined $scoring_module) {
											require $scoring_module;
											$new_score[3]=calc_score(substr( $seq, ($Gposind-4), 30));
											if ($weights{"custom_efficacy"}) {												
												$new_score[2]=$new_score[2]+$weights{"custom_efficacy"}*$new_score[3];
											}else{
												$new_score[2]=$new_score[2]+$new_score[3];
											}	
										}
                                          if(exists($Gpos{$Gposind}) && exists($Gpos{$Gposind+1})){
											if ($weights{"preceedingGG_efficacy"}) {												
												$new_score[2]=$new_score[2]+$weights{"preceedingGG_efficacy"};
											}else{
												$new_score[2]=$new_score[2]+1
											}	
                                          }
                                          if(exists($Gpos{$Gposind})){
                                              if ($weights{"startingG_efficacy"}) {												
												$new_score[2]=$new_score[2]+$weights{"startingG_efficacy"};
											}else{
												$new_score[2]=$new_score[2]+1
											}
                                          }
                                          if(($flank_array[3]+$flank_array[1])>50){
												if ($weights{"totalgc_efficacy"}) {												
													$new_score[2]=$new_score[2]+$weights{"totalgc_efficacy"};
												}else{
													$new_score[2]=$new_score[2]+1
												}
                                          }
                                          @flank_array = find_base_count( substr( $taleseq, $length-7,6) );
                                          if(($flank_array[3]+$flank_array[1])>70){
                                              if ($weights{"seedgc_efficacy"}) {												
													$new_score[2]=$new_score[2]+$weights{"seedgc_efficacy"};
												}else{
													$new_score[2]=$new_score[2]+1
												}
                                          }
                                           if ($weights{"microhom_efficacy"}) {
												$new_score[2]=$new_score[2]+$weights{"microhom_efficacy"}*score_micro_homolgy(\%combined,30,( $Gposind + $length +2 -5 ),5,\$seq);
										   }else{
												$new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Gposind + $length +2 -5 ),5,\$seq);
										   }
                                          ${ $CRISPR_hash{$name} }{"start"} = ($Gposind) + $cut;
                                          ${ $CRISPR_hash{$name} }{"end"} = ( $Gposind + $length + 2 ) + $cut;
                                          ${ $CRISPR_hash{$name} }{"length"} = $length + 2;
                                          my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset - 500;
                                          my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset - 500;
                                          my %score = calculate_CRISPR_score(\%trees, \%something, ($end-5), ($end-5), $chrom, 1, \@new_score , $gene_id);
                                          
                                          #############################################################################################################################################
                                          #Statistics
                                          #############################################################################################################################################
                                          
                                          if (make_CRISPR_statistics(\%something, \%score, $dont_asses_context, \%tempstatistics) == 1){
                                                delete $CRISPR_hash{$name};
                                                next LENGTHLOOP;
                                          }
                                          
                                          if ($something{"retrieve_recomb_matrix"} eq "true") {
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                          }
                                          
                                          %{ ${ $CRISPR_hash{$name} }{"context"} } = %score;
                                          ${ $CRISPR_hash{$name} }{"nucseq"} = $taleseq;
                                          ${ $CRISPR_hash{$name} }{"strand"} = "plus";
                                          $count++;
                                          #############################################################################################################################################
                                          
                                    } else {
                                          $tempstatistics{"Number of designs excluded because their nucleotide composition was too invariable or contained TTTTT"}++;
                                          next LENGTHLOOP;
                                    }
                              }
                        }
                  }
                  
                  #####################################################################################################################################################################
                  # Backward Sequence Calculations
                  #####################################################################################################################################################################
                  
                  if ($something{"preceding"} eq "A") {
                        %dont_care_ind=%Tpos;
                  } elsif($something{"preceding"} eq "G"){
                        %dont_care_ind=%Cpos;
                  } elsif($something{"preceding"} eq "C"){
                        %dont_care_ind=%Gpos;
                  } elsif($something{"preceding"} eq "T"){
                        %dont_care_ind=%Apos;
                  } else{
                        %dont_care_ind=(%Gpos,%Tpos,%Cpos,%Apos);
                  }
                  if ($something{"PAM"} eq "NAG") {
                        %PAMindex=%Tpos;
                  } elsif ($something{"PAM"} eq "NGG") {
                        %PAMindex=%Cpos;
                  } else{
                        %PAMindex=(%Tpos,%Cpos);
                  }
                  POSLOOP: foreach my $Cposind ( sort( keys(%Cpos) ) ) {
                        LENGTHLOOP: foreach my $length ( ($something{"min_length"}+1) .. ($something{"max_length"}+1) ) {
                              if ( exists $PAMindex{ ( $Cposind + 1 ) } && exists $dont_care_ind{ ( $Cposind + $length + 1 ) } ) {
                                    my $taleseq = substr( $seq, $Cposind, $length + 2 );
                                    my @flank_array = find_base_count( $taleseq );
                                    $tempstatistics{"Total number of possible designs"}++;
                                    if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                          $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                          $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                          $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                          !($taleseq=~m/AAAAA/) 
                                    ){
                                         my $name = ($seq_obj->display_id)."_" . $count . "_" . $cut. "." .(int(abs($Cposind + $cut-$start_of_start)/3));
										my @new_score=(0,0,0,0);
										if (defined $scoring_module) {
											require $scoring_module;
											$new_score[3]=calc_score(reverse_comp(substr( $seq, $Cposind-3, 30)));
											if ($weights{"custom_efficacy"}) {												
												$new_score[2]=$new_score[2]+$weights{"custom_efficacy"}*$new_score[3];
											}else{
												$new_score[2]=$new_score[2]+$new_score[3];
											}	
										}
                                          if(exists($Cpos{($Cposind + $length + 1)}) && exists($Cpos{($Cposind + $length)})){
                                               if ($weights{"preceedingGG_efficacy"}) {												
													$new_score[2]=$new_score[2]+$weights{"preceedingGG_efficacy"};
												}else{
													$new_score[2]=$new_score[2]+1
												}
                                          }
                                          if(exists( $Cpos{($Cposind + $length + 1)}) ){
                                                if ($weights{"startingG_efficacy"}) {												
												$new_score[2]=$new_score[2]+$weights{"startingG_efficacy"};
											}else{
												$new_score[2]=$new_score[2]+1
											}
                                          }
                                          if(($flank_array[3]+$flank_array[1])>50){
                                                if ($weights{"totalgc_efficacy"}) {												
													$new_score[2]=$new_score[2]+$weights{"totalgc_efficacy"};
												}else{
													$new_score[2]=$new_score[2]+1
												}
                                          }
                                          @flank_array = find_base_count( substr( $taleseq, 0,6) );
                                          if(($flank_array[3]+$flank_array[1])>70){
                                               if ($weights{"seedgc_efficacy"}) {												
													$new_score[2]=$new_score[2]+$weights{"seedgc_efficacy"};
												}else{
													$new_score[2]=$new_score[2]+1
												}
                                          }
										  if ($weights{"microhom_efficacy"}) {
											$new_score[2]=$new_score[2]+$weights{"microhom_efficacy"}*score_micro_homolgy(\%combined,30,( $Cposind + 5 ),5,\$seq);
										  }else{
											$new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Cposind + 5 ),5,\$seq);
										  }
                                          ${ $CRISPR_hash{$name} }{"start"} = ($Cposind) + $cut;
                                          ${ $CRISPR_hash{$name} }{"end"} = ( $Cposind + $length + 2 ) + $cut;
                                          ${ $CRISPR_hash{$name} }{"length"} = $length + 2;
                                          my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset - 500;
                                          my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset - 500;
                                          my %score = calculate_CRISPR_score(\%trees, \%something, ($end-5), ($end-5), $chrom, 0,\@new_score,$gene_id);
                                          
                                          #############################################################################################################################################
                                          #Statistics
                                          #############################################################################################################################################
                                          
                                          if (make_CRISPR_statistics(\%something, \%score, $dont_asses_context, \%tempstatistics) == 1){
                                                delete $CRISPR_hash{$name};
                                                next LENGTHLOOP;
                                          }
                                          if ($something{"retrieve_recomb_matrix"} eq "true") {
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                          }
                                          
                                          %{ ${ $CRISPR_hash{$name} }{"context"} } = %score;
                                          ${ $CRISPR_hash{$name} }{"nucseq"} = $taleseq;
                                          ${ $CRISPR_hash{$name} }{"strand"} = "minus";
                                          $count++;
                                          
                                          #############################################################################################################################################
                                          
                                    } else {
                                          $tempstatistics{"Number of designs excluded because their nucleotide composition was too invariable or contained TTTTT"}++;
                                          next LENGTHLOOP;
                                    }
                              }
                        }
                  }
            } else{
                  
                  #####################################################################################################################################################################
                  # Double Sequence - only forward calculations needed
                  #####################################################################################################################################################################
                  
                  if ($something{"preceding"} eq "A") {
                        %dont_care_ind=%Tpos;
                        %dont_care_ind_right=%Apos;
                  }elsif($something{"preceding"} eq "G"){
                        %dont_care_ind=%Cpos;
                        %dont_care_ind_right=%Gpos;
                  }elsif($something{"preceding"} eq "C"){
                        %dont_care_ind=%Gpos;
                        %dont_care_ind_right=%Cpos;
                  }elsif($something{"preceding"} eq "T"){
                        %dont_care_ind=%Apos;
                        %dont_care_ind_right=%Tpos;
                  }else{
                        %dont_care_ind=(%Gpos,%Tpos,%Cpos,%Apos);
                        %dont_care_ind_right=(%Gpos,%Tpos,%Cpos,%Apos);
                  }
                  my %PAMindex_right=();
                  if ($something{"PAM"} eq "NAG") {
                        %PAMindex=%Tpos;
                        my %PAMindex_right=%Apos;
                  } elsif ($something{"PAM"} eq "NGG") {
                        %PAMindex=%Cpos;
                        %PAMindex_right=%Gpos;
                  } else{
                        %PAMindex=(%Tpos,%Cpos);                        
                        %PAMindex_right=(%Gpos,%Apos);
                  }
                  POSLOOP: foreach  my $Cposind ( sort( keys(%Cpos) ) ) {
                        SPACERLOOP: foreach my $spacerlength ( ($something{"minspacerlength"}) .. ($something{"maxspacerlength"}) ) {
                              LENGTHLOOP: foreach  my $length ( ($something{"min_length"}+1) .. ($something{"max_length"}+1) ) {
                                    if ( exists $PAMindex{ ( $Cposind + 1 ) } && exists $dont_care_ind{ ( $Cposind + $length + 1 ) } && exists $dont_care_ind_right{ ( $Cposind + $length + 1 + $spacerlength ) } && exists $PAMindex_right{ ( $Cposind + $length +$length + 1 + $spacerlength ) } && exists $Gpos{ ( $Cposind + $length +$length + 2 + $spacerlength ) } ) {
                                          my $left_taleseq = substr( $seq, $Cposind, $length + 2 );
                                          my $right_taleseq = substr( $seq, ( $Cposind + $length + 1 + $spacerlength ), ( $length + 2));
                                          my $completeseq=$left_taleseq.$right_taleseq;
                                          my @flank_array = find_base_count( ($left_taleseq.$right_taleseq) );
                                          $tempstatistics{"Total number of possible designs"}++;
                                          if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                                $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                                $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                                $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                                !($completeseq=~/TTTTT/)
                                          ) {
                                          my $name = ($seq_obj->display_id)."_" . $count . "_" . $cut. "." .(int(abs($Cposind + $cut-$start_of_start)/3));
											 my @new_score=(0,0,0,0);
												if (defined $scoring_module) {
													require $scoring_module;
													$new_score[3]=calc_score(reverse_comp(substr( $seq, $Cposind-3, 30)))+calc_score(substr( $seq, ( $Cposind + $length + 1 + $spacerlength-4),30));
													if ($weights{"custom_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"custom_efficacy"}*$new_score[3];
													}else{
														$new_score[2]=$new_score[2]+$new_score[3];
													}	
												}
                                                if(exists($Cpos{( $Cposind + $length + 1 )}) && exists($Cpos{( $Cposind + $length)})){
                                                    if ($weights{"preceedingGG_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"preceedingGG_efficacy"};
													}else{
														$new_score[2]=$new_score[2]+1
													}
                                                }
                                                if(exists($Cpos{( $Cposind + $length + 1 )})){
                                                   if ($weights{"startingG_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"startingG_efficacy"};
													}else{
														$new_score[2]=$new_score[2]+1
													}
                                                }
                                                if(exists($Gpos{( $Cposind + $length + 1 + $spacerlength )}) && exists($Gpos{( $Cposind + $length + 2 + $spacerlength )})){
                                                    if ($weights{"preceedingGG_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"preceedingGG_efficacy"};
													}else{
														$new_score[2]=$new_score[2]+1
													}
                                                }
                                                if(exists($Gpos{( $Cposind + $length + 1 + $spacerlength )})){
                                                   if ($weights{"startingG_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"startingG_efficacy"};
													}else{
														$new_score[2]=$new_score[2]+1
													}
                                                }
                                                if(($flank_array[3]+$flank_array[1])>50){
                                                    if ($weights{"totalgc_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"totalgc_efficacy"};
													}else{
														$new_score[2]=$new_score[2]+1
													}
                                                }
                                                
                                                @flank_array = find_base_count(substr( $left_taleseq, 0,6));
                                                if(($flank_array[3]+$flank_array[1])>70){
                                                    if ($weights{"seedgc_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"seedgc_efficacy"};
													}else{
														$new_score[2]=$new_score[2]+1
													}
                                                }
                                                
                                                @flank_array = find_base_count(substr( $right_taleseq, $length-7,6));
                                                if(($flank_array[3]+$flank_array[1])>70){
                                                    if ($weights{"seedgc_efficacy"}) {												
														$new_score[2]=$new_score[2]+$weights{"seedgc_efficacy"};
													}else{
														$new_score[2]=$new_score[2]+1
													}
                                                }
												if ($weights{"microhom_efficacy"}) {
													$new_score[2]=$new_score[2]+$weights{"microhom_efficacy"}*score_micro_homolgy(\%combined,30,( $Cposind + $length +5 ),5,\$seq); 
													$new_score[2]=$new_score[2]+$weights{"microhom_efficacy"}*score_micro_homolgy(\%combined,30,( $Cposind + $length +2 + $spacerlength -5 ),5,\$seq);
												}else{
													$new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Cposind + $length +5 ),5,\$seq); 
													$new_score[2]=$new_score[2]+score_micro_homolgy(\%combined,30,( $Cposind + $length +2 + $spacerlength -5 ),5,\$seq);
												}
                                                @{${ $CRISPR_hash{$name} }{"lengthcombo"}}=($length,$spacerlength);
                                                ${ $CRISPR_hash{$name} }{"start"} = ($Cposind) + $cut;
                                                ${ $CRISPR_hash{$name} }{"end"} = ( $Cposind + $length+$spacerlength+$length+2 + 2 ) + $cut;
                                                ${ $CRISPR_hash{$name} }{"length"} =  $length+$spacerlength+$length+2 + 2;
                                                my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset - 500;
                                                my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset - 500;
                                                my %score = calculate_CRISPR_score(\%trees, \%something, ($end-5), ($end-5), $chrom, 0, \@new_score, $gene_id);
                                                
                                                #######################################################################################################################################
                                                #Statistics
                                                #######################################################################################################################################
                                                
                                                if (make_CRISPR_statistics(\%something, \%score, $dont_asses_context, \%tempstatistics) == 1){
                                                      delete $CRISPR_hash{$name};
                                                      next LENGTHLOOP;
                                                }
                                                
                                                if ($something{"retrieve_recomb_matrix"} eq "true" ) {
                                                      ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                      ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                                }
                                                
                                                %{ ${ $CRISPR_hash{$name} }{"context"} } = %score;
                                                @{${ $CRISPR_hash{$name} }{"nucseq"}} = ($left_taleseq,$right_taleseq);
                                                ${ $CRISPR_hash{$name} }{"strand"} = "plus";
                                                $count++;
                                                
                                                #######################################################################################################################################
                                                
                                          } else {
                                                $tempstatistics{"Number of designs excluded because their nucleotide composition was too invariable or contained TTTTT"}++;
                                                next LENGTHLOOP;
                                          }
                                    }
                              }
                        }
                  }
            }
            
            ##########################################################################################################################################################################
            #store the CRISPR and the Statistics in temporary files - for each child process and end the fork
            ##########################################################################################################################################################################
            
            {
                  my $json = JSON::XS::encode_json(\%CRISPR_hash);
                  write_file( $temp_dir . "/" .$seq_obj->display_id . $cut . '.json', { binmode => ':raw' }, $json );
                  $json = JSON::XS::encode_json(\%tempstatistics);
                  write_file( $temp_dir . "/" . $seq_obj->display_id . $cut . 'stats.json', { binmode => ':raw' }, $json );
            }
            
            $pm->finish();
      }
      
      ##########################################################################################################################################################################
      #parent wait till all children are done and then rebuild the CRISPR and the Statistics out of the temporary files
      ##########################################################################################################################################################################
      
      $pm->wait_all_children();
      foreach  my $cut (@cuts) {
            my $json = read_file( $temp_dir . "/" .$seq_obj->display_id . $cut . '.json', { binmode => ':raw' } );
            %finished_CRISPR_hash = ( %finished_CRISPR_hash, %{ decode_json $json } );
            unlink $temp_dir . "/" . $seq_obj->display_id . $cut . ".json";
            $json = read_file( $temp_dir . "/" . $seq_obj->display_id . $cut . 'stats.json', { binmode => ':raw' } );
            my %sechash=%{ decode_json $json };
            foreach  my $seckey (keys(%sechash)){
                  if ( $tempstatistics{$seckey}) {
                        $tempstatistics{$seckey}=$tempstatistics{$seckey}+$sechash{$seckey};
                  }else{$tempstatistics{$seckey}=0}
            }
            unlink $temp_dir . "/" . $seq_obj->display_id . $cut . "stats.json";
      }
      return (\%finished_CRISPR_hash,\%tempstatistics);
}



sub make_database{
        if (can_run('wget') && (can_run('bowtie-build') || can_run('bowtie2-build'))) {                
	    print $_[0]."\n"; #read organism from command-line
	    if(!(-d $_[0])){
			mkdir $_[0];
		}
	    chdir $_[0];		
		system('rsync '.$_[1].'fasta/'.$_[0].'/dna/ > temp.log;');
		# Void context
		if ( fgrep { /README/ } "temp.log" ) {
			print "rsync could connect\nand your files are downloaded";
		}else{
			die "\n\nThere was some problem with the rsync connection to ensembl.\nMaybe you have a typo in the server address or some proxy is hindering the access.\n"
		}
		unlink('temp.log');
	    system('
			rsync -av --progress '.$_[1].'gtf/'.$_[0].'/ .;
			rsync -av --progress --exclude "*abinitio*" '.$_[1].'fasta/'.$_[0].'/cdna/ .;
			rsync -av --progress --exclude "*primary_assembly*" --exclude "*dna_rm*" --exclude "*dna.chromosome*" --exclude "*dna_sm*" '.$_[1].'fasta/'.$_[0].'/dna/ .;
		');
	    print "All files were dowloaded\n";
			system('for f in *.gz ; do gunzip $f; done ;');
			system('for f in *.gtf ; do cat $f > '.$_[0].'.new.gta ; done ;');
			system('for f in *.cdna.all.fa  ; do cat $f > '.$_[0].'.cdna.all.new.fu ; done ;');	    
			system('for f in *.dna.toplevel.fa ; do cat $f > '.$_[0].'.dna.toplevel.new.fu ; done ;');
			system('for f in *.gtf ; do rm $f;  done ;');
			system('for f in *.cdna.all.fa ; do rm $f;  done ;');
			system('for f in *.dna.toplevel.fa ; do rm $f ; done ;');
			system('for f in *.new.gta ; do cat $f > '.$_[0].'.gtf ; done ;');
			system('for f in *.cdna.all.new.fu ; do cat $f > '.$_[0].'.cdna.all.fa ; done ;');
			system('for f in *.dna.toplevel.new.fu; do cat $f > '.$_[0].'.dna.toplevel.fa ; done ;');
			system('for f in *.gta ; do rm $f ; done ;');
			system('for f in *.cdna.all.new.fu ; do rm $f ; done ;');
			system('for f in *.dna.toplevel.new.fu ; do rm $f ; done ;');
			print "All files were unzipped\n";			
			create_mygff($_[0].'.dna.toplevel.fa',$_[0].'.gtf',$_[0].'.all.dna.fa');
			print "All files were converted to gff\n";
			wrap_sequences($_[0].'.cdna.all.fa');
			cpg_for_all($_[0].'.dna.toplevel.fa'); #store CpG-island information in csv-files
			system('for f in *.fasta ; do rm $f ;done ;');
			print "The entire genome was checked for CPG islands\n";
			system('rm *all*.fasta;');
			system('rm *.flat;');
			system('rm *.gdx');
			system('rm \#*');
			system('rm *README*');
			system('rm *CHECKSUMS*');
			correct_cdna($_[0].".cdna.all.fa"); #change header of [organsim].cdna.all.fa in [organsim].cdna.all.facorrected.fa
			system('for f in *.cdna.all.fa ; do rm $f ;done ;');
			print "CDNA files were corrected\n";
			system("mv ".$_[0].".cdna.all.facorrected.fa ".$_[0].".cdna.all.fa;"); #rename [organsim].cdna.all.facorrected.fa
			system('for f in *.gff ; do rm $f ;done ;');
			 print "Annotation were formatted\n";
			include_cpg("."); #add CpG-island information to mygff-files
			system('for f in *.csv ; do rm $f ;done ;');
			system('for f in *.gtf ; do rm $f ;done ;');
			print "All prerequisites for building alignment indeces were built correctly.\nNow Bowtie indeces will be build, depending on the size of the target genome, this can take a while.\n";
			if(can_run('bowtie2-build')){system("bowtie2-build ".$_[0].".cdna.all.fa ".$_[0].".cdna & bowtie2-build ".$_[0].".all.dna.fa ".$_[0].".dna & bowtie2-build ".$_[0].".dna.toplevel.fa ".$_[0].".genome;");} #create files with bowtie2-indices
			if(can_run('bowtie-build')){system("bowtie-build ".$_[0].".cdna.all.fa ".$_[0].".cdna & bowtie-build ".$_[0].".all.dna.fa ".$_[0].".dna & bowtie-build ".$_[0].".dna.toplevel.fa ".$_[0].".genome;"); }#create files with bowtie-indices
			opendir my $curr_dir , ".";
			while (readdir($curr_dir)) {
				if (-z $_) {
					unlink($_);
				}				
			}
			closedir($curr_dir);			
			print "The database for the organism ".$_[0]." has been built in the following path:\n";
			system('pwd');
			}else{
                print "wget and bowtie or bowtie2 need to be installed and executable from the \$PATH variable.\n You can test this by running \"which wget\" and \"which bowtie\" from your terminal."
			}
}
sub create_mygff{
	print "creating Fasta Database....";
      my $db       = Bio::DB::Fasta->new($_[0]);
	  print "done\n";
      my @line;
      my $seqstr  ;
	  $Text::Wrap::columns = 61;
      my $chom;
      my $locus_tag;
      my $note;
      my $id;
      my $chrom_file_fasta;
      my $chrom_file_gff;
      my $temp;
	  my $curr_chrom;
      open(my $in_gtf,"<", $_[1]) or die $!;
      open($chrom_file_fasta, ">",$_[2]) or die $!;
      while (<$in_gtf>) {
       if($_!~/^\#/){
         @line=split("\t",$_);
         if ($line[0] ne $curr_chrom) {    
           if(defined $chrom_file_gff ){close($chrom_file_gff) };   
           open($chrom_file_gff, ">", $line[0]."_indexed.mygff") or die $!;
         }
         $curr_chrom=$line[0];
         if ($line[2]=~/gene/) {
           if (($line[3]-500)>0) {
             $temp=$line[3]-500;
           }    
           $seqstr   = wrap('', '', $db->seq($line[0], ($temp) => ($line[4]+500)));
           $line[8]=~m/gene_id \"(.+?)\"/;
           $id=$1;
           $line[8]=~m/gene_name \"(.+?)\"/;
           $locus_tag=$1;
           $line[8]=~m/gene_biotype \"(.+?)\"/;
           $note=$1;
           print $chrom_file_fasta ">$id locus_tag= $locus_tag;note= $note;chrom:".$line[0].":".$line[3]."..".$line[4]."\n$seqstr\n";
           print $chrom_file_gff "gene_".$id."::".$locus_tag."\t".$line[3]."\t".$line[4]."\n";
         }elsif($line[2]=~m/exon/){ #if feature is a start_codon, stop_codon, CDS or exon
           if($line[8]=~m/gene_id \"(.+?)\".*?transcript_id \"(\S+)\"\;.+?exon_number \"(\S+)\"\;/){ #store transcriptID and exonnumber
              print $chrom_file_gff $line[2]."::".$2."::".$3."::".$1."\t".$line[3]."\t".$line[4]."\n"; #print "feature::transcriptID::exonnumber\t[start]\t[end]
              if ($3==1) {
               if ($line[6] eq "+") {
                 print $chrom_file_gff "TSS::".$2."::".$3."::".$1."\t".$line[3]."\t".$line[3]."\n"; #print "feature::transcriptID::exonnumber\t[start]\t[end]
               }else{
                 print $chrom_file_gff "TSS::".$2."::".$3."::".$1."\t".$line[4]."\t".$line[4]."\n"; #print "feature::transcriptID::exonnumber\t[start]\t[end]
               }       
              }       
           }
         }elsif($line[2]!~m/transcript/){
           if($line[8]=~m/gene_id \"(.+?)\".*?transcript_id \"(\S+)\"\;.+?exon_number \"(\S+)\"\;/){ #store transcriptID and exonnumber
             print $chrom_file_gff $line[2]."::".$2."::".$3."::".$1."\t".$line[3]."\t".$line[4]."\n"; #print "feature::transcriptID::exonnumber\t[start]\t[end]
           }
         }
       }  
      }
      close($in_gtf);
      close($chrom_file_fasta);
	  if(defined $chrom_file_gff ){close($chrom_file_gff) };   
}

sub wrap_sequences{
	$Text::Wrap::columns = 61;
	my $db_file_name=$_[0];
	
	open my ($file_to_wrap),$db_file_name ;
	open my $out_file , ">temp" ;
	while (<$file_to_wrap>) {
		if ($_!~m/>/){
			print $out_file wrap('', '', $_);
		}else{
			print $out_file $_;
		}    
	}
	close $file_to_wrap;
	close $out_file;
	rename("temp",$db_file_name);
}

sub build_tree_new{
        my $file=$_[0]; 
        my $dir = getcwd;
        opendir INDIR, $dir;
        
        foreach my $infile (readdir(INDIR)){
            if($infile=~m/(\S+)\.gff/){	#for every gff-file
                my $outfile=$1."_indexed.mygff";	
                my %chroms=();
                my %trees=();
                my $file=$infile;
                open INFILE ,$file; #open infile [seqname].gff
                open OUTFILE ,">".$outfile; #open outfile [seqname]_indexed.mygff
                #my $count="";
                while(my $line = <INFILE>) {
                        chomp $line;
                        my @line= split("\t", $line); #split the line after tab stops
                        my $chrom=$line[0];
                        $chroms{$chrom}++;
                        my $gene_name="";						
                        my $start=$line[3]; #start position of the feature, with sequence numbering starting at 1
                        my $end=$line[4]; #end position of the feature, with sequence numbering starting at 1
                        if($line[2]=~/gene/ ){ #if feature is a gene
                                if($line=~/.*\s(\S+)\slocus_tag= (\S+);.*/){ #if there is an annotation for the locus_tag of the gene
                                        my $temp="gene_".$1."::".$2;
                                        print OUTFILE $temp,"\t",$start,"\t",$end,"\n"; #print "gene_[geneID]::[locus_tag]\t[start]\t[end]
                                } elsif($line=~/.*\s(\S+)\snote= .*/){ #if there is a note for the gene
                                        my $temp="gene_".$1;
                                        print OUTFILE $temp,"\t",$start,"\t",$end,"\n"; #print "gene_[geneID]\t[start]\t[end]
                                } elsif($line[8]=~/.+/) { #if there is no additional information for the gene
                                        my $temp="gene_".$line[8];
                                        print OUTFILE $temp,"\t",$start,"\t",$end,"\n";	#print "gene_[geneID]\t[start]\t[end]
                                }
                        }
                }
                close(INFILE);
                close(OUTFILE);
                }
        }
        close INDIR;
        
        open INFILE, $file.".gtf";
        my %chrom=();
        my $count=0;
        while(my $line = <INFILE>){
            chomp $line;
            my @line=split("\t",$line);	#split the line after tab stops
            $chrom{$line[0]}++;
            
            if($chrom{$line[0]}==1){
                if ($count==0) { #if this is the first annotation for this seqname
                    open OUTFILE, ">>".$line[0]."_indexed.mygff"; #open outfile [seqname]_indexed.mygff
                }else{
                    close OUTFILE;
                    open OUTFILE, ">>".$line[0]."_indexed.mygff";
                }
            }
            if($line=~m/start_codon/ || $line=~m/stop_codon/ || $line=~m/CDS/ || $line=~m/exon/){ #if feature is a start_codon, stop_codon, CDS or exon
                if($line=~m/transcript_id \"(\S+)\"\;.+?exon_number \"(\S+)\"\;/){ #store transcriptID and exonnumber
                   print OUTFILE $line[2]."::".$1."::".$2."\t".$line[3]."\t".$line[4]."\n"; #print "feature::transcriptID::exonnumber\t[start]\t[end]
                }
            }  
            $count++;
        }
        close OUTFILE;
}

sub correct_cdna{
        my $filename=$_[0];
        open INFILE, $filename;
        open OUTFILE, ">".$filename."corrected.fa";
        
        foreach my $line (<INFILE>){
                $line=~s/^>(\S+).*gene:(\S+).*$/>$2 transcript:$1/;
                print OUTFILE $line;
        }
        close INFILE;
        close OUTFILE;
}
sub cpg_for_all{
        my $pm = Parallel::ForkManager->new($max_parallel);
        my $db       = Bio::DB::Fasta->new($_[0]);
		my @ids      = $db->get_all_primary_ids;
        foreach my $id (@ids){
                my $pid = $pm->start and next;				
					my $seq=$db->seq($id);
					open(my $temp, ">", $id.".fasta") or die $!;
					print $temp $seq;
					close($temp);
					predict_cpg_islands($id.".fasta"); #prints CpG-islands to csv-file
                $pm->finish;             
        }
        $pm->wait_all_children;
}

sub include_cpg{
        opendir INDIR, $_[0];
        foreach my $file (readdir(INDIR)){
            if ($file=~m/(.*)\.csv/) { #for every csv-file
				if (-e $1."_indexed.mygff") {
					open OUTFILE, ">>".$1."_indexed.mygff";
                open INFILE, $file;
                while(my $line = <INFILE>){
                    my @line=split(",",$line); #split the line after commas
                    print OUTFILE "CpG_".$line[1]."_".$line[2]."\t".$line[1]."\t".$line[2]."\n"; #add CpG_[start]_[stop]\t[start]\t[stop] to [seqname]_indexed.mygff
                }
                close INFILE;
                close OUTFILE; 
			}else{
				unlink($file);
			}
				
                
            }
        }
        closedir INDIR;
}

sub convert_to_gff{
        my ($seqfile) = $_[0];
        die("must define a valid seqfile to read") unless ( defined $seqfile && -r $seqfile);
        
        my $pm = Parallel::ForkManager->new($max_parallel);
        my $seqio = new Bio::SeqIO(-format => 'genbank', -file => $seqfile);
        while( my $seq = $seqio->next_seq ) { #for every sequence
			my $fname = $seq->display_id; #store seqname
			if ($fname=~/\#/) {
				next;
			}	
            my $pid = $pm->start and next;
            		
                open SEQOUT, ">$fname.fasta";
                open GFFOUT, ">$fname.gff";
                    foreach my $feature ( $seq->top_SeqFeatures() ) {		
                        if($feature->primary_tag()=~m/gene/){ #if feature is a gene		
                                my $chrom=$fname.":".$feature->start()."..".$feature->end(); #store seqname:[start]..[stop]
                                my $name="";
                                for my $tag ($feature->get_all_tags) {         
                                         for my $value ($feature->get_tag_values($tag)) { #add all tag-value-pairs
                                                if($tag eq "gene"){                
                                                        $name.="$value ";
                                                }else{  
                                                        $name.="$tag= $value;";
                                                }           
                                            }		          
                                }
                                print GFFOUT "$fname\tENSEMBL\tgene\t".$feature->start()."\t".$feature->end()."\t.\t+\t.\t$name\n"; #print genes with attributes to gff-file
                                $name.="chrom:$chrom";	
                                my $currseq=substr($seq->seq(),$feature->start(),($feature->end()-$feature->start()));
                                #$currseq=~s/(\w{5000})/$1\n/g;
                                chomp $currseq;
                                print SEQOUT ">$name\n".$currseq."\n";	#print sequence to fasta-file				
                        }
                    }
                close SEQOUT;
                close GFFOUT;
            $pm->finish;
        }
        $pm->wait_all_children;
}

sub predict_cpg_islands{
        ## LICENCE AGREEMENT
        ##
        ## The intellectual property right of CpG island search script belongs 
        ## to Daiya Takai and Peter A. Jones.  Use of the this script is free for 
        ## academic users. Daiya Takai is responsible for the scientific content 
        ## of the script.
        ##
        ## Integration of this script into another script/program  will require
        ## explicit permission of authors. Such permission will only be granted 
        ## if there is a valid scientific or technical reason to encapsulate the 
        ## entire content of the script into a new resource. 
        ##
        ## Integration of any part of the script into any commercial product is only 
        ## permitted by the agreement with the authors.
        ##
        ## Publication of results obtained with the aid of the script should cite:
        ##
        ## Takai D and Jones PA. Comprehensive analysis of CpG islands in human 
        ## chromoseme 21 and 22. PNAS 2002 99(6):3740-5
        ##
        ## If you have questions, please send an email to the following address: 
        ## Daiya Takai : takai_d@ccnt.hsc.usc.edu
        ##
        ## This script is reviced on May 21st, 2002.
        
        my $GCC=55;
        my $OE=0.65;
        my $LENGTH=500;
        
        our $filename= shift @_;
        while(@_){
                my $commandlineoption=shift @_;
                my @commandlineoption=split ("=",$commandlineoption);
                        if ($commandlineoption[0] =~ /GCC/i){
                                if (($commandlineoption[1] >=50) and ($commandlineoption[1] <=70)){
                                        $GCC=$commandlineoption[1];
                                }else{
                                        $GCC=55;
                                }
                        }elsif($commandlineoption[0] =~ /OE/i){
                                if (($commandlineoption[1] >=0.60) and ($commandlineoption[1] <=1.0)){
                                        $OE=$commandlineoption[1];
                                }else{
                                        $OE=0.65;
                                }
                        }elsif($commandlineoption[0] =~ /LENGTH/i){
                                if (($commandlineoption[1] >=200) and ($commandlineoption[1] <=1500)){
                                        $LENGTH=$commandlineoption[1];
                                }else{
                                        $LENGTH=500;
                                }
                        }
        }
        &extract($filename,$GCC,$OE,$LENGTH);
        
        sub extract($$$$) {
                my $sequence=undef;
                open (IN,$_[0]) || die "$_[0]: $!";
                        while(<IN>){
                                chomp;
                                s/\r\n//g; # remove return
                                s/\r//g;
                                s/\n//g;
                                s/\s//g;
                                s/\d//g;
                                $sequence=$sequence.$_;
                        }
                close(IN);
        my $exportfile="";
        if ($filename=~m/(\S+)\.fasta/) {	#save locusname
                 $exportfile=$1.".csv";		#create csv file [seqname].csv
        } else {
                my @filename_array=split ('\.',$_[0]);
                 $exportfile=$filename_array[0].'.csv';
        }
        
        open(OUT, ">$exportfile") || die "$exportfile: $!";
        
        
        my $resolution=$_[3];													#definition of resolution
        my $window_cg=7;
        my $window_cgr=$_[1];
        my $window_cgs=$_[2];
        my($criteria_cgr)=$_[1];
        my($criteria_cgs)=$_[2];
        my($criteria_length)=$_[3];
        
        
        my $seqlength=length $sequence;
		if (!defined($seqlength)) {
			$seqlength=0;
		}
		
        my $CpGno=1;
        my $CpGstart=0;														#CpGstart+1 is equal to the 
                                                                                                                                                #start position of sequences.
        
                while($CpGstart<=$seqlength-$resolution){						#determine 5' border of CpG island
                        my $sample1=substr $sequence, $CpGstart, $resolution;			#get a sample
                        my ($cgr,$cgs,$cg)=&get_parameter($sample1);
        #		print "$CpGstart\n";
        
                                if (($cgr >= $window_cgr) and ($cgs >= $window_cgs) and ($cg >= $window_cg)){
                                        
                                       my  $start=$CpGstart+1;
					my $faultstart=$start;									#temporal 5' border
                                       my  $faultend=$start+$resolution-1;						#temporal 3' border
                                                if ($CpGstart<=$seqlength-$resolution*2){
                                                        $CpGstart=$CpGstart+$resolution;			#shift the window to next portion
                                                }else{
                                                        $CpGstart=$seqlength-$resolution;
                                                }
        enddefinition:
                                                while($CpGstart<=$seqlength-$resolution){
                                                       my  $sample2=substr $sequence, $CpGstart, $resolution;
                                                       my  ($cgr,$cgs,$cg)=&get_parameter($sample2);
                                                                if( (($cgr < $window_cgr) or ($cgs < $window_cgs) or ($cg < $window_cg)) or (($CpGstart == $seqlength-$resolution) and ($cgr >= $window_cgr) and ($cgs >= $window_cgs) and ($cg >= $window_cg))){
                                                                      my   $tempCpGstart=$CpGstart;			#when the window does not meet the criteria
                                                                                                                                                #then shift the window back to 5' side
                                                                                                                                                #until the window meets the cirteria
                                                                                                                                                
                                                                                while ($tempCpGstart >= $CpGstart-$resolution){
                                                                                     my    $sample3=substr $sequence, $tempCpGstart, $resolution;	
                                                                                      my   ($cgr,$cgs,$cg)=&get_parameter($sample3);
                                                                                                if (($cgr >= $window_cgr) and ($cgs >= $window_cgs) and ($cg >= $window_cg)){
                                                                                                      my   $end=$tempCpGstart+$resolution;
                                                                                                       my  $cpglength=$end-$start+1;
                                                                                                      my   $flag=0;
                                                                                                        
                                                                                                                while($cpglength>=$resolution){
                                                                                                                      my   $sample4=substr $sequence, $start-1, $cpglength;
                                                                                                                      my   ($cgr,$cgs,$cg)=&get_parameter($sample4);
        
                                                                                                                                if (($cgr>=$criteria_cgr) and ($cgs>=$criteria_cgs) and ($cpglength>=$criteria_length)){
                                                                                                                                        $cgr=(int($cgr*10))/10;
                                                                                                                                        $cgs=(int($cgs*1000))/1000;
                                                                                                                                #	print "<FONT SIZE=-1>CpG island $CpGno start=$start,end=$end,\%GC=$cgr,ObsCpG/ExpCpG=$cgs,Length=$cpglength</FONT><BR>\n";
                                                                                                                                        print OUT "$start,$end,$cgr,$cgs,$cpglength\n";
                                                                                                                                        $CpGstart=$end-1;
                                                                                                                                        $CpGno++;
                                                                                                                                        last enddefinition;
                                                                                                                                }elsif (($end<=$faultend) or ($cpglength<=$resolution)){
                                                                                                                                        $cpglength=$faultend-$faultstart+1;
                                                                                                                                     my    $sample5=substr $sequence, $faultstart-1, $cpglength;
                                                                                                                                     my    ($cgr,$cgs,$cg)=&get_parameter($sample5);
                                                                                                                                        $CpGstart=$faultend-1;
                                                                                                                                        $cgr=(int($cgr*10))/10;
                                                                                                                                        $cgs=(int($cgs*1000))/1000;
                                                                                                        #					if ($cpglength<$criteria_length){
                                                                                                        #							last enddefinition;
                                                                                                        #						}elsif(($cgr<$criteria_cgr) or ($cgs<$criteria_cgs)){
                                                                                                        #							$faultstart=$faultstart+$resolution/2;
                                                                                                        #							$faultend=$faultend-$resolution/2;
                                                                                                        #						}else{
                                                                                                                                                                print OUT "$faultstart,$faultend,$cgr,$cgs,$cpglength\n";
                                                                                                                                                                $CpGno++;
                                                                                                                                                                last enddefinition;
                                                                                                        #						}
                                                                                                                                }elsif ($flag==0) {
                                                                                                                                        $end--;		#decrement 3' border
                                                                                                                                        $cpglength=$end-$start+1;
                                                                                                                                        $flag=1;
                                                                                                                                }elsif ($flag==1) {
                                                                                                                                        $start++;	#increment 5' border
                                                                                                                                        $cpglength=$end-$start+1;
                                                                                                                                        $flag=0;
                                                                                                                                }
                                                                                                                }
                                                                                                }
                                                                                }continue{
                                                                                        $tempCpGstart--;
                                                                                }
                                                                }
                                                                
                                                                
                                                                
                                                                
                                                }continue{
                                                        if ($CpGstart <= $seqlength-$resolution*2){
                                                                $CpGstart=$CpGstart+$resolution;
                                                                $faultend=$faultend+$resolution;
                                                        }else{
                                                                $CpGstart=$seqlength-$resolution;
                                                                $faultend=$seqlength;
                                                        }
                                                }
                                }
                }continue{
                        $CpGstart=$CpGstart+1;
                }
                                
        
        close(OUT);
        
        open(IN,$exportfile)|| die "$exportfile: $!";
	my @cpginf=();
       my  $j=0;
                while(<IN>){
                        chomp;
                       my  @temp=split(",",$_);
                                for (my $k=0;$k<=4;$k++){
                                       $cpginf[$j]->[$k]=$temp[$k];
                                }
                        $j++;
                }
        close(IN);
                for (my $l=1;$l<=$j-1;$l++){
                        if ($cpginf[$l]->[0]-$cpginf[$l-1]->[1]<=100){
                                my $length=$cpginf[$l]->[1]-$cpginf[$l-1]->[0]+1;
                                my $samplex=substr $sequence,$cpginf[$l-1]->[0]-1,$length;
                                my ($cgr,$cgs,$cg)=&get_parameter($samplex);
                                        if (($cgr>=$criteria_cgr) and ($cgs>=$criteria_cgs) and ($length>=$criteria_length)){
                
                                                $cpginf[$l-1]->[1]=$cpginf[$l]->[1];
                                                $cpginf[$l-1]->[2]=(int($cgr*10))/10;
                                                $cpginf[$l-1]->[3]=(int($cgs*1000))/1000;
                                                $cpginf[$l-1]->[4]=$length;
                                                
                                                        for (my $p=$l;$p<=$j-1;$p++){
                                                                for (my $q=0;$q<=4;$q++){
                                                                        $cpginf[$p]->[$q]=$cpginf[$p+1]->[$q];
                                                                }
                                                        }
                                                        
                                                        for (my $q=0;$q<=4;$q++){
                                                                $cpginf[$j-1]->[$q]="";
                                                        }
                                                $j--;
                                                $l--;
                                        }
                                
                        }
                }
        
        #print "Selected lower limits: \%GC=$_[1], ObsCpG/ExpCpG=$_[2], Length=$_[3]\n";
        my $contigname=substr $exportfile, 0, 8;
        open(OUT,">$exportfile") || die "$exportfile: $!";
        my $q=1;
                for (my $m=1;$m<=$j;$m++){
                                if ($cpginf[$m-1]->[4]>=$criteria_length){
                                        print OUT "CpG island "."$q,$cpginf[$m-1]->[0],$cpginf[$m-1]->[1],$cpginf[$m-1]->[2],$cpginf[$m-1]->[3],$cpginf[$m-1]->[4]\n";
                                        #print "$contigname, CpG island $q, start=$cpginf[$m-1]->[0], end=$cpginf[$m-1]->[1], \%GC=$cpginf[$m-1]->[2], ObsCpG/ExpCpG=$cpginf[$m-1]->[3], Length=$cpginf[$m-1]->[4]\n";
        #				print "$contigname, CpG island $q, $cpginf[$m-1]->[0], $cpginf[$m-1]->[1], $cpginf[$m-1]->[2], $cpginf[$m-1]->[3], $cpginf[$m-1]->[4]\n";
                                        $q++;
                                }
                }
        if ($q==1){
                #print "No CpG islands was found.\n";
        }
        
        close(OUT);
        #unlink "$exportfile";
        
        sub get_parameter ($){
                my $cgs;
                my($sequence)=$_[0];
                my($seqlength)=length $sequence;
                my($c)=($sequence =~s/c/c/ig);
                my($g)=($sequence =~s/g/g/ig);
                my($cg)=($sequence =~s/cg/cg/ig);
                my($cgr)=($c+$g)/$seqlength*100;
                        if (($c==0) or ($g==0)) {
                                $cgs=0;
                        }else{
                                $cgs=$cg*$seqlength/$c/$g;
                        }
                ($cgr,$cgs,$cg)	
        }
        
        }
		unlink($filename)
}


sub mock{};
=cut
=head1 AUTHOR

Florian Heigwer, C<< <f.heigwer at dkfz.de> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-etools-ecrisp at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=ETOOLS-ECRISP>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.



GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG
=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc ETOOLS::ECRISP


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=ETOOLS-ECRISP>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/ETOOLS-ECRISP>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/ETOOLS-ECRISP>

=item * Search CPAN

L<http://search.cpan.org/dist/ETOOLS-ECRISP/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Florian Heigwer.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

1; # End of ETOOLS::ECRISP
