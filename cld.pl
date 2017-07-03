#!/usr/bin/perl
###################################################################################################################################################################################################
# initialize geneneral cld settings and import perl dependencies
###################################################################################################################################################################################################

#use strict;
#use warnings FATAL => 'all';
use Bio::DB::Fasta;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Scalar::Util qw(looks_like_number);
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
use Parallel::ForkManager; #important package to enable mutlithreading of the script
use Cwd;
use IPC::Cmd qw[can_run run];
use Getopt::Long qw(:config pass_through);	
use File::Grep qw( fgrep fmap fdo );
use Text::Wrap;
use Unix::Processors;
require Tk;							#loads the Tk library, neccessary for producing a simple Tk GUI
Tk->import;
require Tk::PathEntry;					#loads the interactive pathentry widget
require Tk::Dialog;						#loads the Tk Fileopener dialog widget
require Tk::Dressing;
require Tk::BrowseEntry;
require Tk::Optionmenu;
require Tk::Widget;
require Tk::Frame;
require Tk::Entry;
require Tk::Label;
require Tk::Button;
require Tk::Scrollbar;
require Tk::Checkbutton;
require Tk::MainWindow;
require Tk::NoteBook;
require Tk::Text::SuperText;

#open(my $debug_log, ">", "debug_log.txt") or die $!;

my $procs = new Unix::Processors;
my $max_parallel= my $parallel_number =$procs->max_online;
my $aligner_path="";
if(-d $ENV{PAR_TEMP}."/inc/"){
	$aligner_path=$ENV{PAR_TEMP}."/inc/bowtie_progs/";
}
$| = 1;

my ($script_name,$script_version,$script_date,$script_years) = ('cld','1.4.6','2017-02-17','2013-2015');


###################################################################################################################################################################################################
# initialize all GUI and command line variables
###################################################################################################################################################################################################
my(
    %parameter_presets,
    %something,
    $annotation_options_frame,
    $annotation_options_page,
    $chk_CDS_only,
    $chk_CpG_exclusive,
    $annonymous_funct,
    $chk_exclude_overlapping_genes,
    $chk_exon_exclusive,
    $chk_gene_exclusive,
    $chk_ignore_intergenic,
    $chk_ignore_missing_id,
    $chk_out_gff,
    $chk_purpose_exclusive,
    $chk_retrieve_recomb_matrix,
    $chk_sec_off_target,
    $chosen_kit,
    $databasepath_entry,
    $general_options_frame,
    $general_options_page,
    $io_options_frame,
    $io_options_page,
    $key,
    $kit_frame,
    $kit_list,
    $kit_page,
    $lab_crisprai,
    $lab_databasepath,
    $lab_kit_list,
    $lab_knockout,
    $lab_paired,
    $lab_ref_organism,
    $lab_specific_exon,
    $lab_specific_transcript,
    $lab_tagging,
    $list_lab,
    $mw,
    $nb,
    $off_target_options_page,
    $offtarget_options_frame,
    $opt_bowtie_mode,
    $opt_bowtie_version,
    $opt_data_type,
    $gene_list_entry,
    $opt_kind,
    $opt_offtargetdb,
    $opt_PAM,
    $opt_preceding,
    $opt_purpose,
    $param,
    $ref_organism_entry,
    $scl_crispra_downstream,
    $scl_crispra_upstream,
    $scl_crispri_downstream,
    $scl_crispri_upstream,
    $scl_downstream_window,
    $scl_edit_distance_allowed,
    $scl_left_homology,
    $scl_max_A,
    $scl_max_C,
    $scl_max_G,
    $scl_max_length,
    $scl_max_per_exon,
    $scl_max_T,
    $scl_maxspacerlength,
    $scl_min_A,
    $scl_min_C,
    $scl_min_G,
    $scl_min_length,
    $scl_min_T,
    $scl_minspacerlength,
    $scl_number_of_CDS,
    $scl_off_targets_allowed,
    $scl_right_homology,
    $scl_unspecific_leading_bases,
    $scl_upstream_window,
    $specific_exon_entry,
    $specific_options_frame,
    $specific_options_page,
    $specific_options_placeholder_frame,
    $specific_transcripts_entry,
    $Start_but,
    $temp_file,
    $theme,
    $tk_dressing,
    $trigger_frame,
    $trigger_page,
    $higher_order_functions,
    $lower_order_functions,
    $make_database,
    $target_ident,
    $end_to_end,
    $gene_list,
    $database_download_but,
    $database_creation_but,
    $opt_rsync,
    $opt_precalc,
    $lab_gene_entry,
    @string,
    $isok,
    $lab_output_dir,
    $output_dir_entry ,
    $funct_results
);


###################################################################################################################################################################################################
# read in command line options
###################################################################################################################################################################################################

GetOptions(
	    'task=s'		=> \$something{"task"},
	    'output-dir=s'	=> \$something{"working_path"},	
	    'parameter-file=s'	=> \$something{"param_file_name"},
	    'gene-list=s'	=> \$something{"gene_list_file_name"},	        
	    'cov=i'		=> 		\$something{"coverage"},
	    'lib-size=i'	=> \$something{"total_lib_size"},
	    'library_name=s'	=> \$something{"library_name"},
	    '5-prime=s'		=> \$something{"5_adapt"},
	    '3-prime=s'		=> \$something{"3_adapt"},
	    'cor-5-prime=s'	=> \$something{"correct_5_prime_G"},
	    'input-folder=s'	=> \$something{"input_folder"},
	    'organism=s'	=> \$something{"organism_db"},
	    'rsync-link=s'	=> \$something{"rsync_link"},
	    'version'		=> \$something{"version"},
	    'help'		=> \$something{"help"},
		'GUI'		=> \$something{"GUI"},
		'spread-over-transcripts=s'=> \$something{"cover_many_transcripts"},
		'scoring-module=s'=>\$something{"scoring_module"}
	);

###################################################################################################################################################################################################
# define help and version strings
###################################################################################################################################################################################################

$something{"version_string"} = "$script_name, version $script_version, $script_date\nAuthor $script_years Florian Heigwer\n";
$something{"help_string"} = qq{Usage: cld --task=end_to_end [options=value] ...
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
		    --library_name=<string>				- Prefix for the final library as <string> default(test_lib).
		    --5-prime=<string>				- Define the adapter to be put in 5' before the target site.
												default(CTGAGCTCATAGAAGACCTCACC)
		    --3-prime=<string>				- Define the adapter to be put in 3' behind the target site.
								   				default(GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG)
		    --cor-5-prime=<string>			- Specify if the first 5' baspair should be corrected to a G
								   				 may be 1 or 0 default :1.	    
		    --input-folder=<path/to/dir>		- Specify the input folder for library assembly.
								    			this folder must be prepared by --task= target_ident
			--spread-over-transcripts=<string>	- should the designs be equally spread oer the different transcripts of the gene
													-an be : 1 or 0 (default:0)

		 end_to_end 							to perform and end-to-end analysis from target identification to library formatting
		    --output-dir=<path/to/dir>			- a working directory as unix path to directory.
		    --parameter-file=<path/to/dir>		- a parameter file in cld format as path to file.
		    --gene-list=<path/to/dir>			- a gene list file with ENSEMBL IDs new-line seprated as path to file. 
		    --cov=<int>						- Specify the minimum gene coverage as <int> default(15)
		    --lib-size=<int>				- Specify the maximum library size as <int> default(2000)
		    --library_name=<string>				- Prefix for the final library as <string> default(test_lib).
		    --5-prime=<string>				- Define the adapter to be put in 5' before the target site.
								    			default(CTGAGCTCATAGAAGACCTCACC)
		    --3-prime=<string>				- Define the adapter to be put in 3' behind the target site.
								    			default(GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG)
		    --cor-5-prime=<string>			- Specify if the first 5' baspair should be corrected to a G.
			--spread-over-transcripts=<string>	- should the designs be equally spread oer the different transcripts of the gene
													-can be : 1 or 0 (default:1)
			--scoring-module=<path/to/dir>		- the path and filename of a file defining a perl scoring function
			
	    --version							- Show version.
	    --help								- Show this message.
};

print (($something{"version"} ? $something{"version_string"} : ''), ($something{"help"} ? $something{"help_string"} : ''));
###################################################################################################################################################################################################
# read in command line options
###################################################################################################################################################################################################

#if(!can_run('bowtie') ){
#	print "Cannot find bowtie. Please download and install from\n http://bowtie-bio.sourceforge.net/index.shtml\n";
#	print $something{"version_string"}."\n";
#	exit;
#	}elsif(!can_run('bowtie2')){
#	print "Cannot find bowtie2. Please download and install from\n http://bowtie-bio.sourceforge.net/bowtie2/index.shtml\n";
#	print $something{"version_string"}."\n";
#	exit;
#}
if ($something{"GUI"}) {
			###################################################################################################################################################################################################
			# define the default parameters for the program
			###################################################################################################################################################################################################
			
				#hardcode the choosable Kit presets
				 %parameter_presets=();
				%{$parameter_presets{"strict"}}					=(
																	"preceding"=>"G",
																	"PAM"=>"NGG",
																	"min_length"=>20,
																	"max_length"=>20,
																	"min_G"=>1,
																	"max_G"=>80,
																	"min_A"=>1,
																	"max_A"=>80,
																	"min_C"=>1,
																	"max_C"=>80,
																	"min_T"=>1,
																	"max_T"=>80,
																	"bowtie_version"=>"bowtie",
																	"offtargetdb"=>"genomeDNA",
																	"targets-allowed"=>5,
																	"unspecific_leading_bases"=>6,
																	"edit_distance_allowed"=>2,
																	"bowtie_mode"=>"very-sensitive",
																	"ignore_intergenic"=>0,
																	"sec_off_target"=>0,
																	"purpose"=>"knockout",
																	"gene_exclusive"=>1,
																	"exon_exclusive"=>1,
																	"CDS_only"=>1,
																	"specific_exon"=>"any",
																	"specific_transcript"=>"any",
																	"exclude_overlapping_genes"=>1,
																	"CpG_exclusive"=>1,
																	"ignore_missing_id"=>1,
																	"kind"=>"single",
																	"match_info"=>0,
																	"max_per_exon"=>40,
																	"purpose_exclusive"=>1,
																	"downstream_window"=>50,
																	"upstream_window"=>50,
																	"number_of_CDS"=>1,
																	"minspacerlength"=>15,
																	"maxspacerlength"=>17,
																	"crispra_upstream"=>400,
																	"crispra_downstream"=>50,
																	"crispri_upstream"=>400,
																	"crispri_downstream"=>50
																	);
				%{$parameter_presets{"medium"}}					=(
																	"preceding"=>"N",
																	"PAM"=>"NGG",
																	"min_length"=>20,
																	"max_length"=>20,
																	"min_G"=>1,
																	"max_G"=>80,
																	"min_A"=>1,
																	"max_A"=>80,
																	"min_C"=>1,
																	"max_C"=>80,
																	"min_T"=>1,
																	"max_T"=>80,
																	"bowtie_version"=>"bowtie",
																	"offtargetdb"=>"genomeDNA",
																	"targets-allowed"=>10,
																	"unspecific_leading_bases"=>1,
																	"edit_distance_allowed"=>1,
																	"bowtie_mode"=>"sensitive",
																	"ignore_intergenic"=>0,
																	"purpose"=>"knockout",
																	"gene_exclusive"=>1,
																	"exon_exclusive"=>1,
																	"CDS_only"=>0,
																	"exclude_overlapping_genes"=>0,
																	"CpG_exclusive"=>1,
																	"ignore_missing_id"=>1,
																	"kind"=>"single",
																	"match_info"=>0,
																	"max_per_exon"=>40,
																	"purpose_exclusive"=>0,
																	"downstream_window"=>50,
																	"upstream_window"=>50,
																	"number_of_CDS"=>3,
																	"minspacerlength"=>13,
																	"maxspacerlength"=>19,
																	"crispra_upstream"=>500,
																	"crispra_downstream"=>70,
																	"crispri_upstream"=>500,
																	"crispri_downstream"=>70
																	);
				%{$parameter_presets{"relaxed"}}				=(
																	"preceding"=>"N",
																	"PAM"=>"NRG",
																	"min_length"=>20,
																	"max_length"=>20,
																	"min_G"=>1,
																	"max_G"=>90,
																	"min_A"=>1,
																	"max_A"=>90,
																	"min_C"=>1,
																	"max_C"=>90,
																	"min_T"=>1,
																	"max_T"=>90,
																	"bowtie_version"=>"bowtie2",
																	"offtargetdb"=>"gDNA",
																	"targets-allowed"=>20,
																	"unspecific_leading_bases"=>0,
																	"edit_distance_allowed"=>0,
																	"bowtie_mode"=>"fast",
																	"ignore_intergenic"=>1,
																	"purpose"=>"knockout",
																	"gene_exclusive"=>0,
																	"exon_exclusive"=>0,
																	"CDS_only"=>0,
																	"exclude_overlapping_genes"=>0,
																	"CpG_exclusive"=>0,
																	"ignore_missing_id"=>0,
																	"kind"=>"single",
																	"match_info"=>0,
																	"max_per_exon"=>400,
																	"purpose_exclusive"=>0,
																	"downstream_window"=>500,
																	"upstream_window"=>500,
																	"number_of_CDS"=>10,
																	"minspacerlength"=>11,
																	"maxspacerlength"=>21,
																	"crispra_upstream"=>700,
																	"crispra_downstream"=>200,
																	"crispri_upstream"=>700,
																	"crispri_downstream"=>200);
				%{$parameter_presets{"CRISPRi"}}				=(
																	"preceding"=>"G",
																	"PAM"=>"NRG",
																	"min_length"=>20,
																	"max_length"=>20,
																	"min_G"=>1,
																	"max_G"=>90,
																	"min_A"=>1,
																	"max_A"=>90,
																	"min_C"=>1,
																	"max_C"=>90,
																	"min_T"=>1,
																	"max_T"=>90,
																	"bowtie_version"=>"bowtie",
																	"offtargetdb"=>"gDNA",
																	"targets-allowed"=>5,
																	"unspecific_leading_bases"=>3,
																	"edit_distance_allowed"=>2,
																	"bowtie_mode"=>"very-sensitive",
																	"ignore_intergenic"=>0,
																	"sec_off_target"=>0,
																	"purpose"=>"CRISPRi",
																	"gene_exclusive"=>0,
																	"exon_exclusive"=>0,
																	"CDS_only"=>0,
																	"specific_exon"=>"any",
																	"specific_transcript"=>"any",
																	"exclude_overlapping_genes"=>1,
																	"CpG_exclusive"=>0,
																	"ignore_missing_id"=>1,
																	"kind"=>"single",
																	"match_info"=>0,
																	"max_per_exon"=>40,
																	"purpose_exclusive"=>1,
																	"downstream_window"=>50,
																	"upstream_window"=>50,
																	"number_of_CDS"=>1,
																	"minspacerlength"=>15,
																	"maxspacerlength"=>17,
																	"crispra_upstream"=>400,
																	"crispra_downstream"=>50,
																	"crispri_upstream"=>400,
																	"crispri_downstream"=>50
																	);
				%{$parameter_presets{"CRISPRa"}}				=(
																	"preceding"=>"G",
																	"PAM"=>"NRG",
																	"min_length"=>20,
																	"max_length"=>20,
																	"min_G"=>1,
																	"max_G"=>90,
																	"min_A"=>1,
																	"max_A"=>90,
																	"min_C"=>1,
																	"max_C"=>90,
																	"min_T"=>1,
																	"max_T"=>90,
																	"bowtie_version"=>"bowtie",
																	"offtargetdb"=>"gDNA",
																	"targets-allowed"=>5,
																	"unspecific_leading_bases"=>3,
																	"edit_distance_allowed"=>2,
																	"bowtie_mode"=>"very-sensitive",
																	"ignore_intergenic"=>0,
																	"sec_off_target"=>0,
																	"purpose"=>"CRISPRa",
																	"gene_exclusive"=>0,
																	"exon_exclusive"=>0,
																	"CDS_only"=>0,
																	"specific_exon"=>"any",
																	"specific_transcript"=>"any",
																	"exclude_overlapping_genes"=>1,
																	"CpG_exclusive"=>0,
																	"ignore_missing_id"=>1,
																	"kind"=>"single",
																	"match_info"=>0,
																	"max_per_exon"=>40,
																	"purpose_exclusive"=>1,
																	"downstream_window"=>50,
																	"upstream_window"=>50,
																	"number_of_CDS"=>1,
																	"minspacerlength"=>15,
																	"maxspacerlength"=>17,
																	"crispra_upstream"=>400,
																	"crispra_downstream"=>50,
																	"crispri_upstream"=>400,
																	"crispri_downstream"=>50
																	);
				%{$parameter_presets{"non-coding"}}				=(
																	"preceding"=>"G",
																	"PAM"=>"NRG",
																	"min_length"=>20,
																	"max_length"=>20,
																	"min_G"=>1,
																	"max_G"=>90,
																	"min_A"=>1,
																	"max_A"=>90,
																	"min_C"=>1,
																	"max_C"=>90,
																	"min_T"=>1,
																	"max_T"=>90,
																	"bowtie_version"=>"bowtie",
																	"offtargetdb"=>"gDNA",
																	"targets-allowed"=>5,
																	"unspecific_leading_bases"=>3,
																	"edit_distance_allowed"=>2,
																	"bowtie_mode"=>"very-sensitive",
																	"ignore_intergenic"=>0,
																	"sec_off_target"=>0,
																	"purpose"=>"non-coding",
																	"gene_exclusive"=>1,
																	"exon_exclusive"=>0,
																	"CDS_only"=>0,
																	"specific_exon"=>"any",
																	"specific_transcript"=>"any",
																	"exclude_overlapping_genes"=>1,
																	"CpG_exclusive"=>0,
																	"ignore_missing_id"=>1,
																	"kind"=>"single",
																	"match_info"=>0,
																	"max_per_exon"=>400,
																	"purpose_exclusive"=>1,
																	"downstream_window"=>50,
																	"upstream_window"=>50,
																	"number_of_CDS"=>1,
																	"minspacerlength"=>15,
																	"maxspacerlength"=>17,
																	"crispra_upstream"=>400,
																	"crispra_downstream"=>50,
																	"crispri_upstream"=>400,
																	"crispri_downstream"=>50
																	);
				
				#hardcode the parameter defaults for the following scripts										
				#general sgRNA properties											;#
				#$something{"preceding"}="any"										;#
				#$something{"PAM"}="NGG"												;#
				$something{"min_length"}=20											;#
				$something{"max_length"}=20											;#
				$something{"min_G"}=1												;#
				$something{"max_G"}=90												;#
				$something{"min_A"}=1												;#
				$something{"max_A"}=90												;#
				$something{"min_C"}=1												;#
				$something{"max_C"}=90												;#
				$something{"min_T"}=1												;#
				$something{"max_T"}=90												;#
				
				#off-target specific options										;#
				#$something{"bowtie_version"}="bowtie"								;#
				#$something{"offtargetdb"}="gDNA"									;#
				$something{"targets-allowed"}=5									;#
				$something{"unspecific_leading_bases"}=5							;#
				$something{"edit_distance_allowed"}=2								;#
				#$something{"bowtie_mode"}="very-sensitive"							;#
				$something{"ignore_intergenic"}=0									;#
				$something{"sec_off_target"}=0										;#
				
				#locus annotation specific options									;#
				#$something{"purpose"}="knockout"									;#
				$something{"gene_exclusive"}=1										;#
				$something{"exon_exclusive"}=1										;#
				$something{"CDS_only"}=1											;#
				$something{"specific_exon"}="any"									;#
				$something{"specific_transcript"}="any"								;#
				$something{"exclude_overlapping_genes"}=0							;#
				$something{"CpG_exclusive"}=0										;#
				
				#input/output options												;#
				$something{"databasepath"}="/Users/heigwer/Desktop"						;#
				$something{"ref_organism"}="homo_sapiens"							;#
			#	$something{"data_type"}="ensemble_acc"								;#
				$something{"ignore_missing_id"}=1									;#
			#	$something{"kind"}="single"											;#
				$something{"match_info"}=0											;#
				$something{"max_per_exon"}=4000										;#
				$something{"out_gff"}=1												;#
				
				#specific options													;#
				$something{"purpose_exclusive"}=0									;#determines if purpose specific options takes effect
				#options for tagging												;#
				$something{"retrieve_recomb_matrix"}=0								;#
				$something{"right_homology"}=500									;#
				$something{"left_homology"}=500										;#
				#option for tagging and KO											;#
				$something{"downstream_window"}=50									;#
				$something{"upstream_window"}=50									;#
				$something{"number_of_CDS"}=1										;#
				#options for double design after Ran et al. 2014					;#
				$something{"minspacerlength"}=15									;#
				$something{"maxspacerlength"}=17									;#
				#options for CRISPRa/i												;#
				$something{"crispra_upstream"}=400									;#
				$something{"crispra_downstream"}=50									;#
				$something{"crispri_upstream"}=400									;#
				$something{"crispri_downstream"}=50									;#
				$something{"coverage"}=15                                                ;#
				$something{"total_lib_size"}=2000                                         ;#
				$something{"library_name"}="test_lib"                                   ;#
				$something{"5_adapt"}="CTGAGCTCATAGAAGACCTCACC"                     ;#
				$something{"3_adapt"}="TTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG";#
				$something{"correct_5_prime_G"}=0                                         ;#
				$something{"cover_many_transcripts"}=0                             ;#
				$something{"sort_by_rank"}=0;
				$something{"working_path"}=".";
				$something{"custom_score"}=$something{"default_custom_score"}='
			##########################################################################################
			# Adapted to perl after Doench et al. 2014, NAt. Biotechnology
			##########################################################################################
			sub {
				my $score;  
				if (length($_[0])==30) {
				my %sing_nuc_hash = ("G2"=>-0.275377128,"A3"=>-0.323887456,"C3"=>0.172128871,"C4"=>-0.100666209,"C5"=>-0.20180294, 
								"G5"=>0.245956633,"A6"=>0.036440041,"C6"=>0.098376835,"C7"=>-0.741181291,
								"G7"=>-0.393264397,"A12"=>-0.466099015,"A15"=>0.085376945,"C15"=>-0.013813972,
								"A16"=>0.272620512,"C16"=>-0.119022648,"T16"=>-0.285944222,"A17"=>0.097454592,
								"G17"=>-0.17554617,"C18"=>-0.345795451,"G18"=>-0.678096426,"A19"=>0.22508903,
								"C19"=>-0.507794051,"G20"=>-0.417373597,"T20"=>-0.054306959,"G21"=>0.379899366,
								"T21"=>-0.090712644,"C22"=>0.057823319,"T22"=>-0.530567296,"T23"=>-0.877007428,
								"C24"=>-0.876235846,"G24"=>0.278916259,"T24"=>-0.403102218,"A25"=>-0.077300704,
								"C25"=>0.287935617,"T25"=>-0.221637217,"G28"=>-0.689016682,"T28"=>0.117877577,
								"C29"=>-0.160445304,"G30"=>0.386342585);
				my %dinuc_hash = ("GT2"=>-0.625778696,"GC5"=>0.300043317,"AA6"=>-0.834836245,"TA6"=>0.760627772,"GG7"=>-0.490816749,
								  "GG12"=>-1.516907439,"TA12"=>0.7092612,"TC12"=>0.496298609,"TT12"=>-0.586873894,"GG13"=>-0.334563735,
								  "GA14"=>0.76384993,"GC14"=>-0.53702517,"TG17"=>-0.798146133,"GG19"=>-0.66680873,"TC19"=>0.353183252,
								  "CC20"=>0.748072092,"TG20"=>-0.367266772,"AC21"=>0.568209132,"CG21"=>0.329072074,"GA21"=>-0.836456755,
								  "GG21"=>-0.782207584,"TC22"=>-1.029692957,"CG23"=>0.856197823,"CT23"=>-0.463207679,"AA24"=>-0.579492389,
								  "AG24"=>0.649075537,"AG25"=>-0.077300704,"CG25"=>0.287935617,"TG25"=>-0.221637217,"GT27"=>0.117877577,
								  "GG29"=>-0.697740024);
				my $gc = ( substr($_[0],4,20) =~ tr/GC/GC/);
				if ($gc < 10){
					$score=0.597636154+(abs($gc-10)*-0.202625894)
				}else{
					$score=0.597636154+(($gc-10)*-0.166587752)
				}        
				foreach my $i (0..29){        
				   my $key = substr($_[0],$i,1).($i+1);
				   if ($sing_nuc_hash{$key}) {
					$score+=$sing_nuc_hash{$key};
				   }
				   if($i<29){
					$key =substr($_[0],$i,2).($i+1);
					if ($dinuc_hash{$key}){
							$score+=$dinuc_hash{$key};
					}
				   }
				}
				return(1/(1+exp(-$score)))
				  #code
				}else{
					return(0);
				}
			};';
									
					
			###################################################################################################################################################################################################
			###################################################################################################################################################################################################
			
			open( $temp_file, ">", "temp") or die $!;
			print $temp_file ';=================================================
			; Theme  : snow
			; Author : Djibril Ousmanou
			; Date   : 01/01/2011 00:00:00
			;=================================================
			
			[Balloon]
			-background: #FFFFFF
			-foreground: #0047B9
			
			[BrowseEntry]
			-background: #FFFFFF
			-disabledbackground: #FFFFFF
			-disabledforeground: #5F5E5E
			-foreground: #0047B9
			
			[Button]
			-activebackground: #A8A8A7
			-activeforeground: #0047B9
			-background: #A8A8A7
			-disabledforeground: #46562C
			-foreground: #0047B9
			
			[Canvas]
			-background: #FFFFFF
			
			[Checkbutton]
			-activebackground: #FFFFFF
			-activeforeground: #0047B9
			-background: #FFFFFF
			-disabledforeground: #0047B9
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-selectcolor: #0047B9
			
			[ColoredButton]
			-autofit: 1
			-background: #FFFFFF
			-highlightbackground: #FFFFFF
			
			[DirTree]
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			-selectbackground: #A8A8A7
			-selectforeground: #0047B9
			
			[Entry]
			-background: #FFFFFF
			-disabledbackground: #FFFFFF
			-disabledforeground: #5F5E5E
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-insertbackground: #0047B9
			-readonlybackground: #FFFFFF
			-selectbackground: #A8A8A7
			
			[EntryCheck]
			-background: #FFFFFF
			-disabledbackground: #FFFFFF
			-disabledforeground: #5F5E5E
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-insertbackground: #0047B9
			-readonlybackground: #FFFFFF
			-selectbackground: #A8A8A7
			
			[Frame]
			-background: #FFFFFF
			-highlightbackground: #FFFFFF
			
			[HList]
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			-selectbackground: #A8A8A7
			-selectforeground: #0047B9
			
			[LabEntry]
			-background: #FFFFFF
			-disabledbackground: #FFFFFF
			-disabledforeground: #5F5E5E
			
			[LabFrame]
			-background: #FFFFFF
			-foreground: #0047B9
			
			[Label]
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			
			[Labelframe]
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			
			[Listbox]
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			-selectbackground: #A8A8A7
			-selectforeground: #0047B9
			
			[MainWindow]
			-background: #FFFFFF
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			
			[Menu]
			-activebackground: #FFFFFF
			-activeforeground: #0047B9
			-background: #FFFFFF
			-foreground: #0047B9
			-selectcolor: #FFFFFF
			
			[Menubutton]
			-activebackground: #A8A8A7
			-activeforeground: #0047B9
			-background: #A8A8A7
			-disabledforeground: #46562C
			-foreground: #0047B9
			
			[NoteBook]
			-background: #FFFFFF
			-backpagecolor: #FFFFFF
			-disabledforeground: #5F5E5E
			-focuscolor: #FFFFFF
			-foreground: #0047B9
			-inactivebackground: #FFFFFF
			
			[Optionmenu]
			-activebackground: #A8A8A7
			-activeforeground: #0047B9
			-background: #A8A8A7
			-borderwidth: 0
			-disabledforeground: #46562C
			-foreground: #0047B9
			
			[ProgressBar]
			-colors: 0
			-colors: #00FF7B
			-troughcolor: #FFFFFF
			
			[ProgressBarPlus]
			-colors: 0
			-colors: #00FF7B
			-troughcolor: #FFFFFF
			
			[ROText]
			-background: #FFFFFF
			-foreground: #0047B9
			-insertbackground: #0047B9
			
			[Radiobutton]
			-activebackground: #FFFFFF
			-activeforeground: #0047B9
			-background: #FFFFFF
			-disabledforeground: #0047B9
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-selectcolor: #FFFFFF
			
			[Scale]
			-activebackground: #A8A8A7
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			-troughcolor: #FFFFFF
			
			[Spinbox]
			-background: #FFFFFF
			-buttonbackground: #A8A8A7
			-foreground: #0047B9
			-readonlybackground: #FFFFFF
			
			[TList]
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightcolor: #FFFFFF
			-selectbackground: #A8A8A7
			-selectforeground: #0047B9
			
			[Table]
			-background: #FFFFFF
			-highlightbackground: #FFFFFF
			
			[Text]
			-background: #FFFFFF
			-foreground: #0047B9
			-insertbackground: #0047B9
			
			[TextUndo]
			-background: #FFFFFF
			-foreground: #0047B9
			-insertbackground: #0047B9
			
			[Toplevel]
			-background: #FFFFFF
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			
			[Tree]
			-background: #FFFFFF
			-foreground: #0047B9
			-highlightbackground: #FFFFFF
			-highlightcolor: #FFFFFF
			-selectbackground: #A8A8A7
			-selectforeground: #0047B9
			';
									close($temp_file);
				
			  # Set it to e frame widget
			  
			###################################################################################################################################################################################################
			# create a graphical user interface (GUI) to render inputs and output in a human radble form
			###################################################################################################################################################################################################
			  
			 $mw = MainWindow->new( -title => "cld parameter", );
			$mw->minsize( 500, 500 );
			 $higher_order_functions=$mw->NoteBook();
				$make_database = $higher_order_functions ->add('page1', -label => 'Genome Data')->NoteBook();
				$gene_list = $higher_order_functions ->add('page2', -label => 'Gene List')->NoteBook();
				$end_to_end = $higher_order_functions ->add('page3', -label => 'Design Parameter')->NoteBook();
				
				
				$higher_order_functions -> grid(-row=>1,-column=>1,-sticky=>"nw");
				$gene_list		-> grid(-row=>1,-column=>1,-columnspan=>2,-sticky=>"nw");
				$end_to_end		-> grid(-row=>1,-column=>1,-sticky=>"nw");
				$make_database	-> grid(-row=>1,-column=>1,-sticky=>"nw");
				
				$lab_databasepath = $make_database -> Button(-text=>"Select folder for \ndeposition genome data",
																	  -command=>\&choosepathtoindex,
																	-width=>20,-anchor =>"w");	#make the button to trigger the file choosing function
				  $databasepath_entry = $make_database -> PathEntry(
																		 -textvariable=>\$something{"databasepath"},
																		 -width=>40,
																		 -background=>"white"
																		 );					#make the Text/Path entry widget for assisted entry
				  sub choosepathtoindex{							#function triggered by the button to choose an openable existing file
						 $something{"databasepath"}=$mw->chooseDirectory(-title=>"Please choose an folder containing organism data",-initialdir=>"~"); #open the file opening dialog
							$something{"ref_organism"}=$ref_organism_entry->get('1.0','end-1c');
                            if ($something{"databasepath"}=~/[\s]+/) {
                                $something{"databasepath"}="";
                                $databasepath_entry -> configure( -background => "#d73027" );
                                $make_database->messageBox(
                                        -icon => 'error',
                                        -type => 'ok',
                                        -title => 'Error',
                                        -message => 'Folder names may not contain characters other than A-Za-z0-9_-:',
                                    );
                                
                            }
                            
							if(-d $something{"databasepath"}."/".$something{"ref_organism"}){
                                $databasepath_entry -> configure( -background => "#FFFFFF" );
                                $ref_organism_entry -> configure( -background => "#FFFFFF" );
                                }else{
                                $databasepath_entry -> configure( -background => "#d73027" );
                                $ref_organism_entry -> configure( -background => "#d73027" );
							}
						 
				  }
				 
				  $lab_ref_organism  = $make_database -> Label(-text=>"please type the name of the \nreference organism as given in the database\n(e.g. drosophila_melanogaster)",-anchor =>"w");		#create a label object
				  $ref_organism_entry = $make_database -> Text( -height=>1,
																-width=>40);			#create a text entry object
				  
				   $ref_organism_entry->bind('<<Modified>>'=>sub {
							if($ref_organism_entry->editModified) {
								$something{"ref_organism"}=$something{"organism_db"}=$ref_organism_entry->get('1.0','end-1c');
									if(-d $something{"databasepath"}."/".$something{"ref_organism"}){
										$databasepath_entry -> configure( -background => "#FFFFFF" );
										$ref_organism_entry -> configure( -background => "#FFFFFF" );
										}else{
										$databasepath_entry -> configure( -background => "#d73027" );
										$ref_organism_entry -> configure( -background => "#d73027" );
									} 
								$ref_organism_entry->editModified(0)
							}
					});
				   #$ref_organism_entry->insert('1.0',$something{"ref_organism"});
					$opt_rsync=$make_database->Optionmenu(
					   -options => [
										["ensembl-pub"=>"rsync://ftp.ensembl.org/ensembl/pub/release-81/"],
										["ensemblgenomes-protists"=>"rsync://ftp.ensemblgenomes.org/all/pub/protists/current/"],
										["ensemblgenomes-plants"=>"rsync://ftp.ensemblgenomes.org/all/pub/plants/current/"],
										["ensemblgenomes-fungi"=>"rsync://ftp.ensemblgenomes.org/all/pub/fungi/current/"],
										["ensemblgenomes-metazoa"=>"rsync://ftp.ensemblgenomes.org/all/pub/metazoa/current/"]
											],
					   -variable => \$something{"rsync_link"}
					);
					$opt_precalc=$make_database->Optionmenu(
					   -options => [
										["drosophila_melanogaster.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/drosophila_melanogaster.tar.gz"],
										["aedes_aegypti.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/aedes_aegypti.tar.gz"],
										["anopheles_darlingi.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/anopheles_darlingi.tar.gz"],
										["anopheles_gambiae.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/anopheles_gambiae.tar.gz"],
										["arabidopsis_thaliana.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/arabidopsis_thaliana.tar.gz"],
										["aspergillus_nidulans.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/aspergillus_nidulans.tar.gz"],
										["aspergillus_niger.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/aspergillus_niger.tar.gz"],
										["brachypodium_distachyon.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/brachypodium_distachyon.tar.gz"],
										["chlamydomonas_reinhardtii.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/chlamydomonas_reinhardtii.tar.gz"],
										["ciona_intestinalis.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/ciona_intestinalis.tar.gz"],
										["cricetulus_griseus.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/cricetulus_griseus.tar.gz"],
										["dictyostelium_discoideum.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/dictyostelium_discoideum.tar.gz"],
										["drosophila_pseudoobscura.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/drosophila_pseudoobscura.tar.gz"],
										["emiliania_huxleyi.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/emiliania_huxleyi.tar.gz"],
										["homo_sapiens.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/homo_sapiens.tar.gz"],
										["hordeum_vulgare.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/hordeum_vulgare.tar.gz"],
										["komagataella_pastoris.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/komagataella_pastoris.tar.gz"],
										["leishmania_major.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/leishmania_major.tar.gz"],
										["macaca_fascicularis.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/macaca_fascicularis.tar.gz"],
										["macaca_mulatta.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/macaca_mulatta.tar.gz"],
										["magnaporthe_oryzae.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/magnaporthe_oryzae.tar.gz"],
										["mus_musculus.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/mus_musculus.tar.gz"],
										["nematostella_vectensis.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/nematostella_vectensis.tar.gz"],
										["neurospora_crassa.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/neurospora_crassa.tar.gz"],
										["oryctolagus_cuniculus.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/oryctolagus_cuniculus.tar.gz"],
										["oryza_indica.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/oryza_indica.tar.gz"],
										["oryza_nivara.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/oryza_nivara.tar.gz"],
										["oryza_sativa.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/oryza_sativa.tar.gz"],
										["physcomitrella_patens.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/physcomitrella_patens.tar.gz"],
										["plasmodium_falciparum.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/plasmodium_falciparum.tar.gz"],
										["plasmodium_vivax.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/plasmodium_vivax.tar.gz"],
										["populus_trichocarpa.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/populus_trichocarpa.tar.gz"],
										["puccinia_graminis.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/puccinia_graminis.tar.gz"],
										["schizosaccharomyces_pombe.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/schizosaccharomyces_pombe.tar.gz"],
										["taeniopygia_guttata.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/taeniopygia_guttata.tar.gz"],
										["toxoplasma_gondii.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/toxoplasma_gondii.tar.gz"],
										["trichoderma_reesei.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/trichoderma_reesei.tar.gz"],
										["triticum_aestivum.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/triticum_aestivum.tar.gz"],
										["ustilago_maydis.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/ustilago_maydis.tar.gz"],
										["vitis_vinifera.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/vitis_vinifera.tar.gz"],
										["zea_mays.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/zea_mays.tar.gz"],
										["zymoseptoria_tritici.tar.gz"=>"http://www.dkfz.de/signaling/crispr-downloads/DATABASES/zymoseptoria_tritici.tar.gz"]
											],
					   -textvariable => \$something{"precalc"}
					);
					$database_creation_but =$make_database-> Button(	-text=>"Create Genome Data de novo", 	#the button's text
									   -command => sub {
                                                        print "Making a database...\n";
                                                        if (!-d $something{"databasepath"}) {
                                                                $make_database->messageBox(
                                                                    -icon => 'error',
                                                                    -type => 'ok',
                                                                    -title => 'Error',
                                                                    -message => 'The Base folder you selected '.$something{"databasepath"}. "\n does not exist!",
                                                                );
                                                        }else{
                                                            chdir $something{"databasepath"} ;
                                                            make_database(
                                                                defined($something{"ref_organism"}) ? $something{"ref_organism"} : "drosophila_melanogaster",
                                                                defined($something{"rsync_link"}) ? $something{"rsync_link"} : "rsync://ftp.ensembl.org/ensembl/pub/release-77/"
                                                            );
                                                        }
                                            },		#the function triggered by the button
									   -background=>"lightgreen");		#the button's background
					$database_download_but =$make_database-> Button(	-text=>"Download and unzip pre-formed library", 	#the button's text
									   -command => sub {
                                                            $something{"precalc"}=~m/\/([^\/]+)$/;
                                                            my $file_name=$1;
                                                            if (!-d $something{"databasepath"}) {
                                                                $make_database->messageBox(
                                                                    -icon => 'error',
                                                                    -type => 'ok',
                                                                    -title => 'Error',
                                                                    -message => 'The Base folder you selected '.$something{"databasepath"}. "\n does not exist!",
                                                                );
                                                            }else{
                                                                system("curl ".$something{"precalc"}." -o ".$something{"databasepath"}."/$file_name");
                                                                if (!-e $something{"databasepath"}."/$file_name") {
                                                                     $make_database->messageBox(
                                                                         -icon => 'error',
                                                                         -type => 'ok',
                                                                         -title => 'Error',
                                                                         -message => 'The file you selected '.$something{"databasepath"}."/$file_name". "\n failed to download. \n Please check your internet connection and/or proxy settings.",
                                                                     );
                                                                 }else{
                                                                     system("tar -xzvf ".$something{"databasepath"}."/$file_name -C ".$something{"databasepath"});
                                                                    $make_database->messageBox(
                                                                        -icon => 'info',
                                                                        -type => 'ok',
                                                                        -title => 'Infor',
                                                                        -message => 'Your database has been created in\n'.$something{"databasepath"},
                                                                    );
                                                                 }
                                                                
                                                            }
														},		#the function triggered by the button
									   -background=>"lightgreen");		#the button's background
					
					$lab_databasepath		-> grid(-row=>1,-column=>1);
					$databasepath_entry		-> grid(-row=>1,-column=>2,-sticky=>"w");
					$lab_ref_organism		-> grid(-row=>2,-column=>1,-pady =>10);
					$ref_organism_entry		-> grid(-row=>2,-column=>2,-sticky=>"w",-pady =>10);
                   
                    my $lab_option1  = $make_database -> Label(-text=>"Option 1: Choose your organism from the pre-calculated databases",-anchor =>"n");		#create a label object
                    $lab_option1            -> grid(-row=>3,-column=>1,-columnspan=>2,-pady =>10);
                    $opt_precalc			-> grid(-row=>4,-column=>1);
                    $database_download_but		-> grid(-row=>4,-column=>2);
                    
                    my $lab_option2  = $make_database -> Label(-text=>"Option 2: Download your organisms data directly from ensembl and format a database \n (might take some time, requires wget and rsync)",-anchor =>"n");		#create a label object
                    $lab_option2            -> grid(-row=>5,-column=>1,-columnspan=>2,-pady =>10);
                    $database_creation_but		-> grid(-row=>6,-column=>2);
					$opt_rsync				-> grid(-row=>6,-column=>1);
				
					
					
					$gene_list_entry = $gene_list -> Scrolled('Text',
															  -scrollbars => "e",
															  -width=>40
																);
						$gene_list_entry->bind('<<Modified>>'=>sub {
							if($gene_list_entry->editModified) {
																$something{"genes"}=$gene_list_entry->get('1.0','end-1c');
																if ($something{"genes"}=~m/\n/) {
																	@string=split '\n' , $something{"genes"} ;
																	$isok=0;
																	foreach my $sub (@string) {
																		if ($something{"data_type"} eq "ensemble_acc" || $something{"data_type"} eq "gene_symbol") {
																			if ($sub!~/^[\w\.\+\-]+$/) {
																				$isok=1;
																			}  
																		}else{
																			if ($sub!~/^[\w\.\+\-]+\t[\w\.\+\-]+\t\d+\t\d+$/) {
																				$isok=1;
																			}
																		}      
																	}
																	if($isok==0){
																		$gene_list_entry -> configure( -background => "#FFFFFF" );
																	}else{
																		$gene_list_entry -> configure( -background => "#d73027" );
																	}
																}
																
																$gene_list_entry->editModified(0)
																}
						});
					
					$opt_data_type=$gene_list->Optionmenu(
					   -options => [	["gene_symbol" => "gene_symbol"],
                                        ["ensemble_acc" => "ensemble_acc"],
										["coordinates" => "coordinates"]
										],
					   -textvariable => \$something{"data_type"}
					);
					$lab_gene_entry  = $gene_list -> Label(-anchor=>'w',-text=>"Please paste the gene list or\nlist of genomic coordinates here \n\nand select the type of list you chose:\nthe window will turn red if you input was not correct\nallowed are:\nENSG0001234\nENSG001231924\n\nor\ngene1\t5\t12391113\t12391137\ngene2\t5\t12392900\t12392924\ngene3\t5\t12393310\t12393334\ngene4\t5\t12394006\t12394030\n");		#create a label object	
					my $lab_data_type  = $gene_list -> Label(-anchor=>'w',-text=>"Please choose one of the entry options here\n");		#create a label object	
					$lab_data_type          -> grid(-row=>1,-column=>1,-sticky=>"w");
					$opt_data_type			-> grid(-row=>1,-column=>2,-sticky=>"w");
					$lab_gene_entry			-> grid(-row=>2,-column=>1,-sticky=>"w");
					$gene_list_entry			-> grid(-row=>2,-column=>2,-rowspan=>2,-sticky=>"w");
			
					$kit_page=$end_to_end->add('page1', -label => 'General');
					$kit_frame=$kit_page->Frame(-borderwidth => 2, -relief => 'groove')	;
					$lab_kit_list  = $kit_frame -> Label(-text=>"choose a predefined parameter set \nor customise each parameter on the next pages");		#create a label object	
					$kit_list = $kit_frame -> Scrolled(
									   "Listbox", 						#create a scollbar decorated listbox
									   -scrollbars 		=>	"e",				#only one scollbar on the right (east) side should be made
									   -selectmode 		=>	"extended",			#only on can be selected but also unselectd and stuff
									   -activestyle   		=>	"underline",			#the selected entity appears underlined
									   -selectbackground	=>	"lightgreen",			#the selected entity appears in lightgreen
									   -width				=>	20,				#the window has a width of 55 units
									   -background			=>	"white",
									   -height				=> 6
									   );			#the windows has a white background
					#bind a function for the left click on the listbox (what happens when a entity gets selected)
					$kit_list -> bind(
						   '<ButtonPress-1>'=> sub {							#bind a function for the left click on the listbox (what happens when a entity gets selected)
								$chosen_kit=$_[0]->get($_[0]->curselection);				#the selection get read out
                                
							   if (exists $parameter_presets{$chosen_kit}){					#chekc is the selection exists in the preset
                                    my $tmptext="";
                                   foreach  $param (sort keys %{$parameter_presets{$chosen_kit}}){
                                    if($something{$param}!= ${$parameter_presets{$chosen_kit}}{$param} || $something{$param} ne ${$parameter_presets{$chosen_kit}}{$param}){
                                        $tmptext.=$param." = ".${$parameter_presets{$chosen_kit}}{$param}."\n";
									   if ($param eq "PAM") {
										   $opt_PAM -> configure( -textvariable => \$something{$param} );
									   }elsif($param eq "preceding"){
										   $opt_preceding -> configure( -textvariable => \$something{$param} );
									   }elsif($param eq "offtargetdb"){
										   $opt_offtargetdb -> configure( -textvariable => \$something{$param} );
									   }elsif($param eq "bowtie_mode"){
										   $opt_bowtie_mode -> configure( -textvariable => \$something{$param} );
									   }elsif($param eq "purpose"){
										   $opt_purpose -> configure( -textvariable => \$something{$param} );
									   }elsif($param eq "bowtie_version"){
										   $opt_bowtie_version -> configure( -textvariable => \$something{$param} );
									   }
									   $something{$param}= ${$parameter_presets{$chosen_kit}}{$param};
                                    }
								   }
                                   $kit_frame->messageBox(
                                        -icon => 'info',
                                        -type => 'ok',
                                        -title => 'Info',
                                        -message => "Default paramters have been set to the following:\n$tmptext",
                                    );
							   }
						   },
					);
					#write entities in the list beginning with the current "end " of the list so the emptyness
					$kit_list -> insert( 'end' , "strict" , "medium" , "relaxed" , "CRISPRi" , "CRISPRa" , "non-coding");
					my $kind_lab = $kit_frame -> Label(	-text=>"Select if reagents should be designed for \nconventional targeting or paired targeting after Ran et al.",		#the label's text
									   -background=>"white"					#the label's background
					);
					$opt_kind=$kit_frame->Optionmenu(
					   -options => [	["single" => "single"],
									   ["paired" => "paired"]],
					   -textvariable => \$something{"kind"}
					);
					
					$kit_frame 				-> grid(-row=>1,-column=>1,-rowspan=>10,-columnspan=>2,-sticky=>"nw");
					$lab_kit_list			-> grid(-row=>2,-column=>1,-sticky=>"w");
					$kit_list				-> grid(-row=>2,-column=>2,-sticky=>"nw");
                    $kind_lab				-> grid(-row=>1,-column=>1,-sticky=>"nw");
					$opt_kind				-> grid(-row=>1,-column=>2,-sticky=>"nw");
					
					#general options
					$general_options_page=$end_to_end->add('page2', -label => 'Sequence');
					$general_options_frame=$general_options_page->Frame(-borderwidth => 2, -relief => 'groove')									;#
					
					$opt_PAM=$general_options_frame->Optionmenu(
					   -options => [	
									   ["NGG"=>"NGG"],
									   ["NAG"=>"NAG"],
									   ["TTN"=>"TTN"],
									   ["CC"=>"CC"],
									   ["NG"=>"NG"],
									   ["AWG"=>"AWG"],
									   ["TTC"=>"TTC"],
                                       ["TTTV"=>"TTTV"],
									   ["NRG"=>"NRG"],
									   ["CCN"=>"CCN"],
									   ["NGGNG"=>"NGGNG"],
									   ["NAAAAC"=>"NAAAAC"],
									   ["NNGRRT"=>"NNGRRT"],
									   ["NNGRRV"=>"NNGRRV"],
									   ["NNAGAA"=>"NNAGAA"],
									   ["NNNNACA"=>"NNNNACA"],
									   ["NNNNGATT"=>"NNNNGATT"],
									   ["GNNNCNNA"=>"GNNNCNNA"]],
					   -textvariable => \$something{"PAM"}
					)	;
                    my $opt_PAM_location=$general_options_frame->Optionmenu(
					   -options => [	
									   ["3'"=>"3_prime"],
									   ["5'"=>"5_prime"]],
					   -textvariable => \$something{"PAM_location"}
					)	;
					$opt_preceding=$general_options_frame->Optionmenu(
					   -options => [	["G"=>"G"],
									   ["A"=>"A"],
									   ["C"=>"C"],
									   ["T"=>"T"],
									   ["U"=>"U"],
									   ["M"=>"M"],
									   ["R"=>"R"],
									   ["W"=>"W"],
									   ["S"=>"S"],
									   ["Y"=>"Y"],
									   ["K"=>"K"],
									   ["V"=>"V"],
									   ["H"=>"H"],
									   ["D"=>"D"],
									   ["B"=>"B"],
									   ["N"=>"N"]],
					   -textvariable => \$something{"preceding"}
					)	;
					
					$scl_min_length 	= $general_options_frame -> Scale(-label=>"min protospacer [nt]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>13,	-to=>22,	-variable=>\$something{"min_length"}, -tickinterval=>5,	-resolution=>1,-command=> sub {if($something{"min_length"}>$something{"max_length"}){$something{"max_length"}=$something{"min_length"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_max_length 	= $general_options_frame -> Scale(-label=>"max protospacer [nt]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>13,	-to=>22,	-variable=>\$something{"max_length"}, -tickinterval=>5,	-resolution=>1,-command=> sub {if($something{"min_length"}>$something{"max_length"}){$something{"min_length"}=$something{"max_length"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_min_G 	= $general_options_frame -> Scale(-label=>"min guanine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"min_G"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_G"}>$something{"max_G"}){$something{"max_G"}=$something{"min_G"}}}	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_max_G 	= $general_options_frame -> Scale(-label=>"max guanine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"max_G"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_G"}>$something{"max_G"}){$something{"min_G"}=$something{"max_G"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_min_A 	= $general_options_frame -> Scale(-label=>"min adenine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"min_A"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_A"}>$something{"max_A"}){$something{"max_A"}=$something{"min_A"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_max_A 	= $general_options_frame -> Scale(-label=>"max adenine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"max_A"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_A"}>$something{"max_A"}){$something{"min_A"}=$something{"max_A"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_min_C 	= $general_options_frame -> Scale(-label=>"min cytosine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"min_C"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_C"}>$something{"max_C"}){$something{"max_C"}=$something{"min_C"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_max_C 	= $general_options_frame -> Scale(-label=>"max cytosine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"max_C"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_C"}>$something{"max_C"}){$something{"min_C"}=$something{"max_C"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_min_T 	= $general_options_frame -> Scale(-label=>"min thymine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"min_T"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_T"}>$something{"max_T"}){$something{"max_T"}=$something{"min_T"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_max_T 	= $general_options_frame -> Scale(-label=>"max thymine [%]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"max_T"}, -tickinterval=>33,	-resolution=>5 ,-command=> sub {if($something{"min_T"}>$something{"max_T"}){$something{"min_T"}=$something{"max_T"}}});	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					my $PAM_lab = $general_options_frame -> Label(	-text=>"PAM",		#the label's text
									   -background=>"white"					#the label's background
					);
                    my $preceeding_lab = $general_options_frame -> Label(	-text=>"3' base",		#the label's text
									   -background=>"white"					#the label's background
					);
                     my $PAM_loc_lab = $general_options_frame -> Label(	-text=>"PAM location",		#the label's text
									   -background=>"white"					#the label's background
					);
					$general_options_frame 	-> grid(-row=>1,-column=>1,-rowspan=>6,-columnspan=>2,-sticky=>"nw");
					$PAM_lab				-> grid(-row=>1,-column=>3,-sticky=>"nw");
					$opt_PAM				-> grid(-row=>1,-column=>4,-sticky=>"nw");
                    $PAM_loc_lab			-> grid(-row=>1,-column=>5,-sticky=>"nw");
                    $opt_PAM_location       -> grid(-row=>1,-column=>6,-sticky=>"nw");
                    $preceeding_lab			-> grid(-row=>1,-column=>1,-sticky=>"nw");
					$opt_preceding			-> grid(-row=>1,-column=>2,-sticky=>"nw");
					
					$scl_min_length			-> grid(-row=>2,-column=>1,-sticky=>"nw",-columnspan=>2);
					$scl_max_length			-> grid(-row=>2,-column=>3,-sticky=>"nw",-columnspan=>2);
					$scl_min_G				-> grid(-row=>3,-column=>1,-sticky=>"nw",-columnspan=>2);
					$scl_max_G				-> grid(-row=>3,-column=>3,-sticky=>"nw",-columnspan=>2);
					$scl_min_A				-> grid(-row=>4,-column=>1,-sticky=>"nw",-columnspan=>2);
					$scl_max_A				-> grid(-row=>4,-column=>3,-sticky=>"nw",-columnspan=>2);
					$scl_min_C				-> grid(-row=>5,-column=>1,-sticky=>"nw",-columnspan=>2);
					$scl_max_C				-> grid(-row=>5,-column=>3,-sticky=>"nw",-columnspan=>2);
					$scl_min_T				-> grid(-row=>6,-column=>1,-sticky=>"nw",-columnspan=>2);
					$scl_max_T				-> grid(-row=>6,-column=>3,-sticky=>"nw",-columnspan=>2);
                    #on-target specific options
					
					#off-target specific options
					$off_target_options_page=$end_to_end->add('page3', -label => 'Specificity score');
					$offtarget_options_frame=$off_target_options_page->Frame(-borderwidth => 2, -relief => 'groove');
					$opt_bowtie_mode=$offtarget_options_frame->Optionmenu(
					   -options => [	["very-sensitive" => "very-sensitive"],
									   ["sensitive" => "sensitive"],
									   ["fast" => "fast"],
									   ["very-fast" => "very-fast"]],
					   -textvariable => \$something{"bowtie_mode"}
					)	;#
					$opt_bowtie_version=$offtarget_options_frame->Optionmenu(
					   -options => [	["bowtie" => "bowtie"],
									   ["bowtie2" => "bowtie2"],
                                       ["blast" => "blast"]
                                       ],
					   -textvariable => \$something{"bowtie_version"},
					   -command => sub {
								   if($something{"bowtie_version"} eq "bowtie2"){
										   $opt_bowtie_mode					-> grid(-row=>1,-column=>5,-sticky=>"nw");
									   }else{
										   $opt_bowtie_mode					-> gridRemove();
									   }
								   }
					)																	;
					$opt_offtargetdb=$offtarget_options_frame->Optionmenu(
					   -options => [	["all genomic DNA" => "genomeDNA"],
                                        ["annotated genes" => "gDNA"],
                                        ["cDNA" => "cDNA"]
                                        ],
					   -textvariable => \$something{"offtargetdb"}
					)																	;
																				   ;
					$scl_off_targets_allowed 	= $offtarget_options_frame -> Scale(-label=>"targets allowed [#]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>1,	-to=>100,	-variable=>\$something{"targets-allowed"}, -tickinterval=>33,	-resolution=>1	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_unspecific_leading_bases 	= $offtarget_options_frame -> Scale(-label=>"unspecific PAM distal bases [nt]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>0,	-to=>10,	-variable=>\$something{"unspecific_leading_bases"}, -tickinterval=>3,	-resolution=>1	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_edit_distance_allowed 	= $offtarget_options_frame -> Scale(-label=>"mismatches allowed [nt]:",-orient=>'h',	-length=>150,	-digit=>1,	-from=>0,	-to=>5,	-variable=>\$something{"edit_distance_allowed"}, -tickinterval=>2,	-resolution=>1	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					
					$chk_ignore_intergenic = $offtarget_options_frame 		-> Checkbutton(
																				   -text=>"ignore intergenic off targets",	#the button's text
																				   -variable=>\$something{"ignore_intergenic"}
																   );			#the button's background
					$chk_sec_off_target = $offtarget_options_frame 		-> Checkbutton(
																				   -text=>"check off targets in a secondary fasta file",	#the button's text
																				   -variable=>\$something{"sec_off_target"},
                                                                                   -command=>sub{
                                                                                        if (!-e $something{"databasepath"}."/secondary_off_targets.fasta") {
                                                                                            $offtarget_options_frame->messageBox(
                                                                                                    -icon => 'error',
                                                                                                    -type => 'ok',
                                                                                                    -title => 'Error',
                                                                                                    -message => "secondary_off_targets.fasta\nneeds to be provided in ".$something{"databasepath"},
                                                                                                );
                                                                                        }
                                                                                        
                                                                                   }
																   );			#the button's background						;#
					 my $bowtie_version_lab = $offtarget_options_frame -> Label(	-text=>"Off-target search program",		#the label's text
									   -background=>"white"					#the label's background
					);
                      my $offtargetdb_lab = $offtarget_options_frame -> Label(	-text=>"Off-target reference",		#the label's text
									   -background=>"white"					#the label's background
					);
                       my $secondary_lab = $offtarget_options_frame -> Label(	-text=>"(a file named secondary_off_targets.fasta\nneeds to be placed in your genome data path)",		#the label's text
									   -background=>"white"					#the label's background
					);
					$offtarget_options_frame 			-> grid(-row=>1,-column=>1,-columnspan=>5,-sticky=>"nw");
					
                    $bowtie_version_lab 				-> grid(-row=>1,-column=>1,-sticky=>"nw");
                    $opt_bowtie_version 				-> grid(-row=>1,-column=>2,-sticky=>"nw");
					
                    $offtargetdb_lab					-> grid(-row=>1,-column=>3,-sticky=>"nw");
                    $opt_offtargetdb					-> grid(-row=>1,-column=>4,-sticky=>"nw");
					
                    $chk_ignore_intergenic				-> grid(-row=>2,-column=>1,-sticky=>"nw",-columnspan=>2,-pady =>10);
					$chk_sec_off_target					-> grid(-row=>2,-column=>3,-sticky=>"nw",-columnspan=>2,-pady =>10);
                    $secondary_lab					    -> grid(-row=>3,-column=>3,-sticky=>"nw",-columnspan=>2);
					
                    $scl_off_targets_allowed			-> grid(-row=>4,-column=>1,-sticky=>"nw");
					$scl_unspecific_leading_bases		-> grid(-row=>4,-column=>2,-sticky=>"nw");
					$scl_edit_distance_allowed			-> grid(-row=>4,-column=>3,-sticky=>"nw");
					
					#locus annotation specific options									;#
					$annotation_options_page=$end_to_end->add('page4', -label => 'Annotation score');
					$annotation_options_frame=$annotation_options_page->Frame(-borderwidth => 2, -relief => 'groove')							;#
					$opt_purpose=$annotation_options_frame->Optionmenu(
					   -options => [	["knockout" => "knockout"],
									   ["CRISPRi" => "CRISPRi"],
									   ["CRISPRa" => "CRISPRa"],
									   ["non-coding" => "non-coding"],
									   ["c-tagging" => "c-tagging"],
									   ["n-tagging" => "n-tagging"]],
					   -textvariable => \$something{"purpose"}
					);
					$chk_gene_exclusive = $annotation_options_frame 		-> Checkbutton(
																				   -text=>"target only inside genes",	#the button's text
																				   -variable=>\$something{"gene_exclusive"}
																   );
					$chk_exon_exclusive = $annotation_options_frame 		-> Checkbutton(
																				   -text=>"target only inside exons",	#the button's text
																				   -variable=>\$something{"exon_exclusive"}
																   );
					$chk_CDS_only = $annotation_options_frame 		-> Checkbutton(
																				   -text=>"target only inside coding regions",	#the button's text
																				   -variable=>\$something{"CDS_only"}
																   );
					$chk_exclude_overlapping_genes = $annotation_options_frame 		-> Checkbutton(
																				   -text=>"exclude designs in overlapping genes",	#the button's text
																				   -variable=>\$something{"exclude_overlapping_genes"}
																   );
					$chk_CpG_exclusive = $annotation_options_frame 		-> Checkbutton(
																				   -text=>"exclude designs in CpG islands",	#the button's text
																				   -variable=>\$something{"CpG_exclusive"}
																   );
					
					#text entry for the exon specificity
					$lab_specific_exon  = $annotation_options_frame -> Label(-text=>"enter a number if a specific exon should be targeted:");		#create a label object
					$specific_exon_entry = $annotation_options_frame -> Entry(-textvariable=>\$something{"specific_exon"},-width=>20);			#create a text entry object
					
					$chk_purpose_exclusive = $annotation_options_frame 		-> Checkbutton(
																				   -text=>"select if designs should be purpose specific",	#the button's text
																				   -variable=>\$something{"purpose_exclusive"},
																				   -command => sub {
																							if($something{"purpose_exclusive"} == 1){
																									$opt_purpose		 				-> grid(-row=>2,-column=>2,-sticky=>"nw");
																								}else{
																									$opt_purpose					-> gridRemove();
																								}
																							}
																   );
					$annotation_options_frame 			-> grid(-row=>1,-column=>1,-sticky=>"nw");
					
					$chk_purpose_exclusive		 		-> grid(-row=>1,-column=>1,-sticky=>"nw");
					$chk_gene_exclusive					-> grid(-row=>3,-column=>1,-sticky=>"nw");
					$chk_exon_exclusive					-> grid(-row=>4,-column=>1,-sticky=>"nw");
					$chk_CDS_only						-> grid(-row=>5,-column=>1,-sticky=>"nw");
					$chk_exclude_overlapping_genes		-> grid(-row=>6,-column=>1,-sticky=>"nw");
					$chk_CpG_exclusive					-> grid(-row=>7,-column=>1,-sticky=>"nw");
					$lab_specific_exon					-> grid(-row=>8,-column=>1,-sticky=>"nw");
					$specific_exon_entry				-> grid(-row=>8,-column=>2,-sticky=>"nw");
					my $scoring_page=$end_to_end->add('page5', -label => 'On-target score');
						my $scoring_frame=$scoring_page->Frame(-borderwidth => 2, -relief => 'groove');
							my $lab_custom_score= $scoring_frame -> Label(-text=>"Enter a specific scoring funtion\n Only funtions are allowed which:\nhave 30 bp including target site as input\nhave a numeric score between 0 and 1 as output");		#create a label object	
							my $custom_score_entry = $scoring_frame -> Scrolled('SuperText',-width=>50,-height=>30,-wrap=>"none");
							$custom_score_entry->insert('1.0',$something{"custom_score"});
							my $test_but =$scoring_frame-> Button(	-text=>"Test and set function", 	#the button's text
										   -command => sub {
											
												$something{"custom_score"}=$custom_score_entry->get('1.0','end-1c');
												$annonymous_funct=eval($something{"custom_score"});
												if ($annonymous_funct) {
													my $funct_seq=rndStr(30, 'A','C','G','T');
													$funct_results=$funct_seq." : ".$annonymous_funct->($funct_seq);
												}else{
													$funct_results="Function did not return result";
												}
											
											});
							my $default_but =$scoring_frame-> Button(	-text=>"Reset function to default", 	#the button's text
										   -command => sub {
											$custom_score_entry->delete('1.0','end');
											$custom_score_entry->insert('1.0',$something{"default_custom_score"});
											$funct_results="";
											});
							my $lab_results= $scoring_frame -> Entry(-textvariable=>\$funct_results,-width=>50);		#create a label object	
							my $lab_score= $scoring_frame -> Label(-text=>"Choose a Score:");		#create a label object		
							my $opt_scores=$scoring_frame->Optionmenu(
							   -options => [	["Xu" => "xu_score"],
												["Doench (2014)" => "doench_old"],
												["Customized" => "custom"]
											   ],
							   -variable => \$something{"scores"},
							   -command => sub {
									if($something{"scores"} eq "custom"){
											$lab_custom_score                   -> grid(-row=>3,-column=>1,-columnspan=>2,-sticky=>"nw");
											$test_but                           -> grid(-row=>4,-column=>1,-columnspan=>1,-sticky=>"nw");
											$default_but                        -> grid(-row=>4,-column=>2,-columnspan=>1,-sticky=>"nw");
											$lab_results                        -> grid(-row=>5,-column=>1,-columnspan=>2,-sticky=>"nw");
											$custom_score_entry		 			-> grid(-row=>6,-column=>1,-columnspan=>2,-sticky=>"nw");                                             
											
										}else{
											$lab_custom_score                   -> gridRemove();
											$test_but                           -> gridRemove();
											$default_but                        -> gridRemove();
											$lab_results                        -> gridRemove();
											$custom_score_entry		 			-> gridRemove();
										}
									}
							);
							my $chk_sort_by_scores = $scoring_frame 		-> Checkbutton(-text=>"select if designs should be ranked by an efficacy score",
																					   -variable=>\$something{"sort_by_rank"},
																					   -command => sub {
																									if($something{"sort_by_rank"} == 1){
																											$opt_scores		 -> grid(-row=>2,-column=>2,-sticky=>"nw");
																											$lab_score      -> grid(-row=>2,-column=>1,-sticky=>"nw");
																										}else{
																											$opt_scores             -> gridRemove();
																											$custom_score_entry    -> gridRemove();
																											$lab_custom_score		-> gridRemove();
																											$lab_score				-> gridRemove();
                                                                                                            $test_but                           -> gridRemove();
                                                                                                            $default_but                        -> gridRemove();
                                                                                                            $lab_results                        -> gridRemove();
																										}
																									});
							
							$scoring_frame-> grid(-row=>1,-column=>1,-sticky=>"nw");
							$chk_sort_by_scores-> grid(-row=>1,-column=>1,-columnspan=>2,-sticky=>"nw");
					#input/output options												;#
					$io_options_page=$end_to_end->add('page6', -label => 'Input/Output');
					$io_options_frame=$io_options_page->Frame(-borderwidth => 2, -relief => 'groove');#
					$lab_output_dir = $io_options_frame -> Button(-text=>"Select output directory",
																	  -command=>\&choosepathtoindex2,
																	-width=>20,-anchor =>"w");	#make the button to trigger the file choosing function
						$output_dir_entry = $io_options_frame -> PathEntry(
																			   -textvariable=>\$something{"working_path"},
																			   -width=>40,
																			   -background=>"white"
																			   );					#make the Text/Path entry widget for assisted entry
						sub choosepathtoindex2{							#function triggered by the button to choose an openable existing file
							   $something{"working_path"}=$io_options_frame->chooseDirectory(-title=>"Please choose the output directory",-initialdir=>"~"); #open the file opening dialog
						}
					$chk_ignore_missing_id = $io_options_frame 		-> Checkbutton(
																				   -text=>"ignore if ids could not be found",	#the button's text
																				   -variable=>\$something{"ignore_missing_id"}
																   );
					$chk_out_gff = $io_options_frame 		-> Checkbutton(
																				   -text=>"print gff output",	#the button's text
																				   -variable=>\$something{"out_gff"}
																   );
					$scl_max_per_exon 			= $io_options_frame -> Scale(-label=>"max number of designs per exon:",-orient=>'h',	-length=>600,	-digit=>1,	-from=>1,	-to=>4000,	-variable=>\$something{"max_per_exon"}, -tickinterval=>1000,	-resolution=>10	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					
					$io_options_frame 			-> grid(-row=>1,-column=>1,-rowspan=>6,-columnspan=>2,-sticky=>"nw");
					$lab_output_dir             -> grid(-row=>2,-column=>1,-sticky=>"nw");	#make the button to trigger the file choosing function
					$output_dir_entry           -> grid(-row=>2,-column=>2,-sticky=>"nw");
					$chk_ignore_missing_id		-> grid(-row=>3,-column=>1,-sticky=>"nw");
					$chk_out_gff		 		-> grid(-row=>4,-column=>1,-sticky=>"nw");
					#$scl_max_per_exon		 	-> grid(-row=>5,-column=>1,-columnspan=>2,-sticky=>"nw");
					#specific options
					$specific_options_page=$end_to_end->add('page7', -label => 'Specific Options');
					$specific_options_frame=$specific_options_page->Frame(-borderwidth => 2, -relief => 'groove')									;#
					
					$lab_tagging  = $specific_options_frame -> Label(-text=>"Select specific options for tagging proteins:");		#create a label object	
					$chk_retrieve_recomb_matrix = $specific_options_frame 		-> Checkbutton(
																				   -text=>"select if the sequences flanking the target site should be retrieved",	#the button's text
																				   -variable=>\$something{"retrieve_recomb_matrix"}
																   );
					$scl_right_homology 			= $specific_options_frame -> Scale(-label=>"3' homology length:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>400,	-variable=>\$something{"right_homology"}, -tickinterval=>100,	-resolution=>10	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_left_homology 			= $specific_options_frame -> Scale(-label=>"5' homology length:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>400,	-variable=>\$something{"left_homology"}, -tickinterval=>100,	-resolution=>10	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					
					$lab_knockout  	= $specific_options_frame -> Label(-text=>"Select specific options for knockout purpose:");		#create a label object
					
					$scl_downstream_window 		= $specific_options_frame -> Scale(-label=>"bp downstream of start codon:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>50,	-variable=>\$something{"downstream_window"}, -tickinterval=>20,	-resolution=>1	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_upstream_window 		= $specific_options_frame -> Scale(-label=>"bp upstream of start codon:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>50,	-variable=>\$something{"upstream_window"}, -tickinterval=>20,	-resolution=>1	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_number_of_CDS 			= $specific_options_frame -> Scale(-label=>"# exons downstream of the start codon:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>50,	-variable=>\$something{"number_of_CDS"}, -tickinterval=>20,	-resolution=>1	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					
					$lab_paired  	= $specific_options_frame -> Label(-text=>"Select specific options for paired design purpose:");		#create a label object
					$scl_minspacerlength			= $specific_options_frame -> Scale(-label=>"min spacer-length:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>50,	-variable=>\$something{"minspacerlength"}, -tickinterval=>20,	-resolution=>1 ,-command=> sub {if($something{"minspacerlength"}>$something{"maxspacerlength"}){$something{"maxspacerlength"}=$something{"minspacerlength"}}}	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_maxspacerlength 		= $specific_options_frame -> Scale(-label=>"max spacer-length:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>50,	-variable=>\$something{"maxspacerlength"}, -tickinterval=>20,	-resolution=>1 ,-command=> sub {if($something{"minspacerlength"}>$something{"maxspacerlength"}){$something{"minspacerlength"}=$something{"maxspacerlength"}}}	);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					
					my $lab_crispra  	= $specific_options_frame -> Label(-text=>"Select specific options for CRISPRa design purpose:");		#create a label object
					$scl_crispra_upstream		= $specific_options_frame -> Scale(-label=>"bp upstream of TSS:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>500,	-variable=>\$something{"crispra_upstream"}, -tickinterval=>100,	-resolution=>10);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					$scl_crispra_downstream 		= $specific_options_frame -> Scale(-label=>"bp downstream of TSS:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>500,	-variable=>\$something{"crispra_downstream"}, -tickinterval=>100,	-resolution=>10);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					my $lab_crispri  	= $specific_options_frame -> Label(-text=>"Select specific options for CRISPRi design purpose:");		#create a label object
					my $scl_crispri_upstream		= $specific_options_frame -> Scale(-label=>"bp upstream of TSS:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>500,	-variable=>\$something{"crispri_upstream"}, -tickinterval=>100,	-resolution=>10);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					my $scl_crispri_downstream 		= $specific_options_frame -> Scale(-label=>"bp downstream of TSS:",-orient=>'h',	-length=>200,	-digit=>1,	-from=>1,	-to=>500,	-variable=>\$something{"crispri_downstream"}, -tickinterval=>100,	-resolution=>10);	#make slider with a minimum, a maximum,an orientation, a modifed variable, a certain bckground color, the intervals ticks should be drawn, the number by that the interval should be devided 1..100 will be devided by ten so that every tenth value can be selected
					
                    $specific_options_placeholder_frame=$specific_options_frame->Frame(-height=>10)		;
					$specific_options_frame		-> grid(-row=>1,-column=>1,-rowspan=>14,-columnspan=>2,-sticky=>"nw");
					$specific_options_placeholder_frame		-> grid(-row=>1,-column=>1,-sticky=>"nw");
					$chk_purpose_exclusive		-> grid(-row=>2,-column=>1,-sticky=>"nw");
					$lab_tagging				-> grid(-row=>3,-column=>1,-columnspan=>2,-sticky=>"nw",-pady =>20);
					$chk_retrieve_recomb_matrix	-> grid(-row=>4,-column=>1,-columnspan=>2,-sticky=>"nw");
					$scl_right_homology			-> grid(-row=>5,-column=>1,-sticky=>"nw");
					$scl_left_homology			-> grid(-row=>5,-column=>2,-sticky=>"nw");
                    $scl_downstream_window		-> grid(-row=>6,-column=>1,-sticky=>"nw");
					$scl_upstream_window		-> grid(-row=>6,-column=>2,-sticky=>"nw");                    
					$lab_knockout				-> grid(-row=>7,-column=>1,-columnspan=>2,-sticky=>"nw",-pady =>20);					
					$scl_number_of_CDS			-> grid(-row=>8,-column=>1,-sticky=>"nw");
					$lab_paired					-> grid(-row=>9,-column=>1,-columnspan=>2,-sticky=>"nw",-pady =>20);
					$scl_minspacerlength		-> grid(-row=>10,-column=>1,-sticky=>"nw");
					$scl_maxspacerlength		-> grid(-row=>10,-column=>2,-sticky=>"nw");
					$lab_crispra				-> grid(-row=>11,-column=>1,-columnspan=>2,-sticky=>"nw",-pady =>20);
					$scl_crispra_upstream		-> grid(-row=>12,-column=>1,-sticky=>"nw");
					$scl_crispra_downstream		-> grid(-row=>12,-column=>2,-sticky=>"nw");
                    $lab_crispri				-> grid(-row=>13,-column=>1,-columnspan=>2,-sticky=>"nw",-pady =>20);
					$scl_crispri_upstream		-> grid(-row=>14,-column=>1,-sticky=>"nw");
					$scl_crispri_downstream		-> grid(-row=>14,-column=>2,-sticky=>"nw");
					
					my $library_page=$end_to_end->add('page8', -label => 'Library');
						my $library_frame=$library_page->Frame(-borderwidth => 2, -relief => 'groove');
							
                            my  $lab_cov  	= $library_frame -> Label(-text=>"Enter number of sgRNAs per gene:");		#create a label object
							my $cov_backup=20;
							my $cov_entry = $library_frame -> Text(-width=>20,-height=>1);
								$cov_entry->bind('<<Modified>>'=>sub {
									if($cov_entry->editModified) {
                                            $something{"coverage"}=$cov_entry->get('1.0','end-1c');
                                            if ($something{"coverage"}=~m/\D/) {
                                                $something{"coverage"}=$cov_backup;
                                                $cov_entry->delete('0.0','end');
                                                $cov_entry->insert('1.0',$something{"coverage"});
                                            }
                                            $cov_entry->editModified(0)
                                    }
								});
                                
                            my  $lab_lib_size  	= $library_frame -> Label(-text=>"Enter library size in total number of sgRNAs:");		#create a label object
							my $backup=2000;
							my $lib_size_entry = $library_frame -> Text(-width=>20,-height=>1);
								$lib_size_entry->bind('<<Modified>>'=>sub {
									if($lib_size_entry->editModified) {
																		$something{"total_lib_size"}=$lib_size_entry->get('1.0','end-1c');
																		if ($something{"total_lib_size"}=~m/\D/) {
																			$something{"total_lib_size"}=$backup;
																			$lib_size_entry->delete('0.0','end');
																			$lib_size_entry->insert('1.0',$something{"total_lib_size"});
																		}
																		$lib_size_entry->editModified(0)
																		}
								});
							my  $lab_name_entry 	= $library_frame -> Label(-text=>"Enter a name for this library:");		#create a label object		
							my $lib_name_entry = $library_frame -> Entry(-textvariable=>\$something{"library_name"},-width=>20);       
							my $x5_prime_entry = $library_frame -> Scrolled('Entry',-textvariable=>\$something{"5_adapt"},-width=>60);
							my $x3_prime_entry = $library_frame -> Scrolled('Entry',-textvariable=>\$something{"3_adapt"},-width=>60);
							my  $lab_x5_prime 	= $library_frame -> Label(-text=>"Enter cloning adapter 5' of target sequence:");		#create a label object
							my  $lab_x3_prime  	= $library_frame -> Label(-text=>"Enter cloning adapter 3' of target sequence:");		#create a label object
							my $chk_cor_5_prime = $library_frame 		-> Checkbutton(-text=>"Define if the 5' end of each target site should be corrected to a G",	-variable=>\$something{"correct_5_prime_G"});
							my $chk_spread_over_transcripts = $library_frame 		-> Checkbutton(-text=>"should the designs be equally spread over different transcripts of the gene",	-variable=>\$something{"cover_many_transcripts"});
							
							 $library_frame                 -> grid(-row=>1,-column=>1,-columnspan=>2,-sticky=>"nw");
							 $lab_name_entry                -> grid(-row=>1,-column=>1,-sticky=>"nw");
							 $lib_name_entry                -> grid(-row=>1,-column=>2,-sticky=>"nw");
							 $lib_size_entry                -> grid(-row=>2,-column=>2,-sticky=>"w");         
							 $lab_lib_size                  -> grid(-row=>2,-column=>1,-sticky=>"nw");
							
                             $cov_entry                     -> grid(-row=>3,-column=>2,-sticky=>"w");         
							 $lab_cov                       -> grid(-row=>3,-column=>1,-sticky=>"nw");
							 $lab_x5_prime                  -> grid(-row=>4,-column=>1,-sticky=>"nw");
							 $x5_prime_entry                -> grid(-row=>4,-column=>2,-sticky=>"nw");
							 $lab_x3_prime                  -> grid(-row=>5,-column=>1,-sticky=>"nw");
							 $x3_prime_entry                -> grid(-row=>5,-column=>2,-sticky=>"nw");
							 $chk_cor_5_prime               -> grid(-row=>7,-column=>1,-columnspan=>2,-sticky=>"nw");
							 $chk_spread_over_transcripts   -> grid(-row=>8,-column=>1,-columnspan=>2,-sticky=>"nw");
					
					
			###########################################################################################################################################################################
			 ###########################################################################################################################################################################
					
			###########################################################################################################################################################################
					$trigger_page=$end_to_end->add('page9', -label => 'Start Analysis');
					$trigger_frame=$trigger_page->Frame(-borderwidth => 2, -relief => 'groove');#
					#create button objects that will trigger the search program
					$Start_but =$trigger_frame-> Button(	-text=>"Perform\nanalysis", 	#the button's text
									   -command => sub {
                                                    if (!-d $something{"databasepath"}."/".$something{"ref_organism"}) {
                                                        $offtarget_options_frame->messageBox(
                                                                -icon => 'error',
                                                                -type => 'ok',
                                                                -title => 'Error',
                                                                -message => "Please select a reference organism with a database in ".$something{"databasepath"} ,
                                                            );
                                                    }else{
                                                        filter_library(make_a_crispr_library(),
																   defined($something{"library_name"}) ? $something{"library_name"} : "test_lib",
																	defined($something{"coverage"}) ? $something{"coverage"} : 15,
																	defined($something{"total_lib_size"}) ? $something{"total_lib_size"} : 200000, 
																	defined($something{"5_adapt"}) ? $something{"5_adapt"} : "CTGAGCTCATAGAAGACCTCACC", 
																	defined($something{"3_adapt"}) ? $something{"3_adapt"} : "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG",
																	defined($something{"correct_5_prime_G"}) ? $something{"correct_5_prime_G"} : 1, 
																	$something{"working_path"},
																	$something{"gene_list_file_name"},
																	defined($something{"cover_many_transcripts"}) ? $something{"cover_many_transcripts"} : 0,
                                                                    %something
																	);
                                                        }
                                                    });
                                                    
													
					$trigger_frame 			-> grid(-row=>1,-column=>1,-rowspan=>6,-columnspan=>2,-sticky=>"nw");
					$Start_but		-> grid(-row=>1,-column=>1,-sticky=>"nw");
					my $results_text = $trigger_frame -> Scrolled('Text',
															  -scrollbars => "se",
															  -width=>100,
															  -wrap => "none"
																);
				tie *STDOUT, 'Tk::Text', $results_text;
				$results_text		-> grid(-row=>2,-column=>1);
					####################################################################################################################################################################################################################################
				 $tk_dressing = Tk::Dressing->new();
				 $theme = 'theme';#$tk_dressing->get_current_theme;
				$tk_dressing->load_theme_file($theme, "temp");
				$tk_dressing->design_widget(
					-widget => $mw,
					-theme  => $theme,
				   );
			MainLoop();
}else{
		for (my $i=0; $i<scalar(@ARGV); $i++){
			if (substr($ARGV[$i],0,1) eq '-' and $i < scalar(@ARGV)-1){
			if ($ARGV[$i] eq '-task'		){$something{"task"}		= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-output-dir'		){$something{"working_path"}= int($ARGV[++$i]); }	
			if ($ARGV[$i] eq '-parameter-file'	){$something{"param_file_name"}	= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-gene-list'		){$something{"gene_list_file_name"}	= int($ARGV[++$i]); }	        
			if ($ARGV[$i] eq '-cov'			){$something{"coverage"}		= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-lib-size'		){$something{"total_lib_size"}	= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-library_name'		){$something{"library_name"}	= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-5-prime'		){$something{"5_adapt"}		= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-3-prime'		){$something{"3_adapt"}		= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-cor-5-prime'		){$something{"correct_5_prime_G"}	= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-organism'		){$something{"organism_db"}		= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-rsync-link'		){$something{"rsync_link"}	= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-input-folder'	){$something{"input_folder"}		= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-spread-over-transcripts'	){$something{"cover_many_transcripts"}		= int($ARGV[++$i]); }
			if ($ARGV[$i] eq '-scoring-module'	){$something{"scoring_module"}		= int($ARGV[++$i]); }
			}
		}
		
		if(!defined($something{"task"})){
		   print "\nThere must be any specific task specified.\n\n".$something{"version_string"};
		   exit;
		}
		
		
		
		if($something{"task"} eq "make_database"){   
			make_database(
				  defined($something{"organism_db"}) ? $something{"organism_db"} : "drosophila_melanogaster",
				  defined($something{"rsync_link"}) ? $something{"rsync_link"} : "rsync://ftp.ensembl.org/ensembl/pub/release-77/"
				  );
		}elsif($something{"task"} eq "target_ident"){
			if(!defined($something{"gene_list_file_name"})){
				print "\nThere must be an gene list specified.\n\n".$something{"version_string"};
			   exit;
			}
			if(!defined($something{"gene_list_file_name"}) or !(-f $something{"gene_list_file_name"})){ die "The gene list file ".$something{"gene_list_file_name"}." could not be opened. Either the user has no rights the read it or the file does not exist." }
			if(!defined($something{"param_file_name"}) or !(-f $something{"param_file_name"})){ die "The parameter list file ".$something{"param_file_name"}." could not be opened. Either the user has no rights the read it or the file does not exist." }
			if(!defined($something{"working_path"}) or !(-d $something{"working_path"})){ $something{"working_path"}="." }
			make_a_crispr_library(
						 $something{"param_file_name"},
						 $something{"gene_list_file_name"},
						 $something{"working_path"}
					   ),
		}elsif($something{"task"} eq "end_to_end"){
			if(!defined($something{"gene_list_file_name"})){
				print "\nThere must be an gene list specified.\n\n".$something{"version_string"};
			   exit;
			}
			if(!defined($something{"gene_list_file_name"}) or !(-f $something{"gene_list_file_name"})){ die "The gene list file ".$something{"gene_list_file_name"}." could not be opened. Either the user has no rights the read it or the file does not exist." }
			if(!defined($something{"param_file_name"}) or !(-f $something{"param_file_name"})){ die "The parameter list file ".$something{"param_file_name"}." could not be opened. Either the user has no rights the read it or the file does not exist." }
			if(!defined($something{"working_path"}) or !(-d $something{"working_path"})){ $something{"working_path"}="." }
			
			if(defined($something{"coverage"})		and	$something{"coverage"}		=~m/([^\d]+)/g){ die "The gene coverage must be specified as an integer. Your input is not integer as it contains: $1." }
			if(defined($something{"total_lib_size"})	and	$something{"total_lib_size"}	=~m/([^\d]+)/g){ die "The total library size must be specified as an integer. Your input is not integer as it contains: $1." }
			if(defined($something{"5_adapt"})		and	$something{"5_adapt"}		=~m/([^ACGT]+)/g){ die "The 5' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
			if(defined($something{"3_adapt"})		and	$something{"3_adapt"}		=~m/([^ACGT]+)/g){ die "The 3' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
			if(defined($something{"correct_5_prime_G"})	and	$something{"correct_5_prime_G"} ne 1 and $something{"correct_5_prime_G"} ne 0 ){ die "Correcting the 5' basepair to a G must be 1 or 0" }    
			
			filter_library(		
					make_a_crispr_library(
								$something{"param_file_name"},#"params.txt",
								$something{"gene_list_file_name"},# "still_left.tab"
								$something{"working_path"} #~
							),
					defined($something{"library_name"}) ? $something{"library_name"} : "test_lib",
					defined($something{"coverage"}) ? $something{"coverage"} : 15,
					defined($something{"total_lib_size"}) ? $something{"total_lib_size"} : 2000, 
					defined($something{"5_adapt"}) ? $something{"5_adapt"} : "CTGAGCTCATAGAAGACCTCACC", 
					defined($something{"3_adapt"}) ? $something{"3_adapt"} : "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG",
					defined($something{"correct_5_prime_G"}) ? $something{"correct_5_prime_G"} : 1, 
					$something{"working_path"},
					$something{"gene_list_file_name"},
					defined($something{"cover_many_transcripts"}) ? $something{"cover_many_transcripts"} : 0,
                    %something
				);
		}elsif($something{"task"} eq "library_assembly"){
			if(!defined($something{"gene_list_file_name"})){
				print "\nThere must be an gene list specified.\n\n".$something{"version_string"};
			   exit;
			}
			if(!defined($something{"gene_list_file_name"}) or !(-f $something{"gene_list_file_name"})){ die "The gene list file ".$something{"gene_list_file_name"}." could not be opened. Either the user has no rights the read it or the file does not exist." }
			if(!defined($something{"param_file_name"}) or !(-f $something{"param_file_name"})){ die "The parameter list file ".$something{"param_file_name"}." could not be opened. Either the user has no rights the read it or the file does not exist." }
			if(!defined($something{"working_path"}) or !(-d $something{"working_path"})){ $something{"working_path"}="." }
			if(!defined($something{"input_folder"}) or !(-d $something{"input_folder"})){die "The input folder ".$something{"input_folder"}." could not be opened. Either the user has no rights the read it or the file does not exist." }
			if(defined($something{"coverage"})		and	$something{"coverage"}		=~m/([^\d]+)/g){ die "The gene coverage must be specified as an integer. Your input is not integer as it contains: $1." }
			if(defined($something{"total_lib_size"})	and	$something{"total_lib_size"}	=~m/([^\d]+)/g){ die "The total library size must be specified as an integer. Your input is not integer as it contains: $1." }
			if(defined($something{"5_adapt"})		and	$something{"5_adapt"}		=~m/([^ACGT]+)/g){ die "The 5' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
			if(defined($something{"3_adapt"})		and	$something{"3_adapt"}		=~m/([^ACGT]+)/g){ die "The 3' adapter sequence must contain only ACGT charcters. Your input contains: $1." }
			if(defined($something{"correct_5_prime_G"})	and	$something{"correct_5_prime_G"} ne 1 and $something{"correct_5_prime_G"} ne 0 ){ die "Correcting the 5' basepair to a G must be 1 or 0" }    
			
			filter_library(
					$something{"input_folder"},#"Thu_Dec_4_17:28:18_20141417710498",
					defined($something{"library_name"}) ? $something{"library_name"} : "test_lib",
					defined($something{"coverage"}) ? $something{"coverage"} : 15,
					defined($something{"total_lib_size"}) ? $something{"total_lib_size"} : 200000, 
					defined($something{"5_adapt"}) ? $something{"5_adapt"} : "CTGAGCTCATAGAAGACCTCACC", 
					defined($something{"3_adapt"}) ? $something{"3_adapt"} : "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTGGGTCTTCGTTCG",
					defined($something{"correct_5_prime_G"}) ? $something{"correct_5_prime_G"} : 1, 
					$something{"working_path"},
					$something{"gene_list_file_name"},
					defined($something{"cover_many_transcripts"}) ? $something{"cover_many_transcripts"} : 0,
                    %something
				   );
		}else{
			die "Some parameters are missing!\nOne of the four options make_database , target_ident , end_to_end or library assembly must be set.\n";
		}	
}



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
      my $strand="+";
        $score{"CRISPRi"}=0;
        $score{"CRISPRa"}=0;
      my $expression="[";
      if (($_[1]->{"number_of_CDS"}>0) && ($_[1]->{"purpose_exclusive"}==1)) {
         foreach my $number(1..$_[1]->{"number_of_CDS"}){
            $expression.=$number;
            }
         $expression.="]";
      }else{
             $expression="\\d*?";
      }
     
	  if ( exists $_[0]->{$_[4]} ) { # check wethere the tree exists
            #search for annotations in the intervall from start (2) to end (3) and count them in score
            my $annotations = $_[0]->{$_[4]}->fetch( int($_[2]), int($_[3]) );            
			my $tmp;
            my $gene_tmp;
            my $tmpstr;
            foreach  my $anno ( @{$annotations} ) {
                  if ( $anno =~ m/gene_(\S+)::([\+\-])_(\d+)_(\d+)/ ) {
                        $new_score[1]++;
                        $tmp=$1;
                        $tmpstr=$2;
                        ${ $score{"gene"} }{$tmp}++;
                        if ($tmp=~m/$gene_name/) {
                              ${ $score{"this_gene"} }{$tmp}++;
                              $gene_tmp=$tmp;
                              $strand=$tmpstr;
                        }           
                  } elsif ( $anno =~ m/exon::(\S+)::(\d+)::(\S+)\_(\d+)_(\d+)/) {
                        ${ $score{"exon"} }{$2}++;
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
                  } elsif ( $anno =~ m/CDS::(\S+)::($expression)::(\S+)\_(\d+)_(\d+)/ ) {
                        $new_score[1]=$new_score[1]+5/$2;
                        $score{"CDS"}++;
                        if($_[1]->{"specific_transcript"} ne "any"){
                              if ($1 eq $_[1]->{"specific_transcript"}) {
                                    $score{"CDS_1"}++;
                              }
                        } else{
                              $score{"CDS_1"}++;
                        }
                  } elsif ( $anno =~ m/CDS/ ) {
                        $score{"CDS"}++;
                        $new_score[1]++;
                  }
            }
            
           
            if ($strand eq "+") {
                $annotations = $_[0]->{$_[4]}->fetch( int($_[2]-$_[1]->{"crispri_downstream"}),
                                                     int($_[3]+$_[1]->{"crispri_upstream"})
                                                     );     
            }else{
                $annotations = $_[0]->{$_[4]}->fetch( int($_[2]-$_[1]->{"crispri_upstream"}),
                                                     int($_[3]+$_[1]->{"crispri_downstream"})
                                                     );     
            }
            foreach  my $anno ( @{$annotations} ) {
                 if ( $anno =~ m/TSS::\S+::(\S+)_(\d+)_(\d+)/ ) {
                    $tmp=$1;
                    #if ($gene_tmp=~m/$tmp/) {
                        $score{"CRISPRi"}=1;
                    #}           
              }
            }
           if ($strand eq "+") {
                $annotations = $_[0]->{$_[4]}->fetch( int($_[2]-$_[1]->{"crispra_downstream"}),
                                                     int($_[3]+$_[1]->{"crispra_upstream"})
                                                     );     
            }else{
                $annotations = $_[0]->{$_[4]}->fetch( int($_[2]-$_[1]->{"crispra_upstream"}),
                                                     int($_[3]+$_[1]->{"crispra_downstream"})
                                                     );     
            }
            foreach  my $anno ( @{$annotations} ) {
                 if ( $anno =~ m/TSS::\S+::(\S+)_(\d+)_(\d+)/ ) {
                    $tmp=$1;
                    #if ($gene_tmp=~m/$tmp/) {
                        $score{"CRISPRa"}=1;
                    #}           
              }
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
      if ( $_[0]->{"gene_exclusive"} ==1 && !exists $_[1]->{"this_gene"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they did not hit any gene"}++;
            return 1;
      }
	  if ( $_[0]->{"exclude_overlapping_genes"} ==1  && scalar(keys(%{$_[1]->{"gene"}}))>1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they hit multiple genes"}++;
            return 1;
      }
      if (  $_[0]->{"exon_exclusive"} ==1 && !exists $_[1]->{"exon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they did not hit any exon"}++;
            return 1;
      }
      if (  $_[0]->{"CpG_exclusive"} ==1 && exists $_[1]->{"CpG"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were located in an CpG island"}++;
            return 1;
      }
      if (  $_[0]->{"purpose_exclusive"} ==1 && ($_[0]->{"purpose"} eq "knockout") && !(exists $_[1]->{"CDS_1"}) && ($_[2] != 1) ) {
            $_[3]->{"Number of designs excluded because they were not directly behind the ATG of the specified transcript"}++;
            return 1;
      }
      if (  $_[0]->{"purpose_exclusive"} ==1 && ($_[0]->{"purpose"} eq "n-tagging") && !exists $_[1]->{"start_codon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located at the start codon"}++;
            return 1;
      }
      if (  $_[0]->{"purpose_exclusive"} ==1 && ($_[0]->{"purpose"} eq "c-tagging") && !exists $_[1]->{"stop_codon"} && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located at the stop codon"}++;
            return 1;
      }
		if ( $_[0]->{"purpose_exclusive"} ==1 && ($_[0]->{"purpose"} eq "CRISPRa") && $_[1]->{"CRISPRa"}!=1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not amenable for CRISPRa"}++;
            return 1;
      }
      if ( $_[0]->{"purpose_exclusive"} ==1 && ($_[0]->{"purpose"} eq "CRISPRi") && $_[1]->{"CRISPRi"}!=1 && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not amenable for CRISPRi"}++;
            return 1;
      }
      if (  $_[0]->{"CDS_only"} ==1 && !(exists $_[1]->{"CDS"}) && $_[2] != 1 ) {
            $_[3]->{"Number of designs excluded because they were not located in a coding sequence"}++;
            return 1;
      }
      if (  $_[0]->{"exon_exclusive"} ==1 && !( $_[0]->{"specific_exon"} eq "any") && $_[2] != 1  ) {
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
      if($_[3] eq "3_prime"){
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
      }else{
        if ($_[2] eq "fw") {            
            foreach (1..$_[1]){
                    push @matchstring , "n";
              }
        }
        else {
              foreach (1..$_[1]){
                  unshift @matchstring , "n";
            }
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
	my $output_dir=$_[7];
	my $gene_list_file=$_[8];
	my $many_transcripts=$_[9];
    my %parameters=$_[10];
    my @info=();
	my @ids;
    my @line=();
    my %designs=();
    my %genes=();
	my %genes_from_list=();
	
	if ($output_dir!~/\/$/) {
        $output_dir.="/";
    }
    my $gene="";
	
	if (defined $something{"GUI"}) {
		if ($something{"data_type"} eq "ensemble_acc" || $something{"data_type"} eq "gene_symbol") {
            @ids=split("\n",$something{"genes"});
			foreach my $id (@ids){
				$genes_from_list{$id}++;
			}
        }else{
			@ids=split("\n",$something{"genes"});
			foreach my $id (@ids){
				my @coords=split("\t",$id);
				$genes_from_list{$coords[0]}++;
			}
		}   
	}else{
		open (my $gene_list, "<", $gene_list_file) or die $!;
		while(<$gene_list>){
			chomp $_;
			$genes_from_list{$_}++;
		}
	close $gene_list;
	}
	
	my %genes_avail;
    open (my $libgff, "<", $temp_dir . "/all_results_together.gff") or die $!;
       while(<$libgff>){
            if(!($_=~/\#/)){
				my $cline = $_;
                @line=split("\t",$cline);
                @info=split(";\ ",$line[8]);
				foreach my $key (keys(%genes_from_list)){
					if ($cline=~/id=$key/) {
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
	my $gene_name;
    my @new_file=(); 
    my %ids=();
    my %count=();
    my %missing=();
	my %all_ids;
    my $general_coverage=0;
	my %id_for_lib;
	my %id_with_info;
	open (my $libtab, "<", $temp_dir . "/all_results_together.tab") or die $!;
	while(<$libtab>){		
        @line=split("\t",$_);
		if($line[0]=~m/(\S+)_\d+_\d+$/){
			$gene_name=$1;
			my @transcripts=split("_",$line[8]);
			foreach my $trans (@transcripts){
				${${$ids{$gene_name}}{$trans}}{$line[0]}++;
			}
			${${$id_with_info{$gene_name}}{$line[0]}}{"anno_score"}=$line[14];
			${${$id_with_info{$gene_name}}{$line[0]}}{"spec_score"}=$line[13];
			${${$id_with_info{$gene_name}}{$line[0]}}{"custom_score"}=$line[15];
			${${$id_with_info{$gene_name}}{$line[0]}}{"seq"}=$line[6];
		}
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
	if ($many_transcripts ==1) {
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
					sort { if($something{"sort_by_rank"}==1){${$id_with_info{$key}}{$b}->{"custom_score"} <=> ${$id_with_info{$key}}{$a}->{"custom_score"} }else{return 1} }
					keys %{$id_with_info{$key}}
					){
					if( 		$count{$key} < $coverage  
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
		if (($general_coverage+scalar keys %{$id_for_lib{$key}})>$limit) {
			delete $id_for_lib{$key};
			$isgone{$key}="was cut because of overall library size\n";
		}else{
			$general_coverage+=scalar keys %{$id_for_lib{$key}}
		}		
	}

    $lib_name=$lib_name.".".$script_version.".".$coverage.".".scalar(keys(%id_for_lib)).".".$general_coverage;
    open (my $libtab_out, ">", $output_dir.$lib_name.".tab") or die $!;
    open (my $libfa_out, ">", $output_dir.$lib_name.".fasta") or die $!;
    open ($libtab, "<", $temp_dir . "/all_results_together.tab") or die $!;
    my %fasta=();
        while (<$libtab>) {
			@line=split("\t",$_);
			if($line[0]=~m/(\S+)_\d+_\d+$/){
				$gene_name=$1;
				if (${$id_for_lib{$gene_name}}{$line[0]}) {
					$all_ids{$line[0]}++;
					if (!exists($fasta{$line[0]})) {
						print $libtab_out $_;
                        if($something{"PAM_location"} eq "3_prime"){ 
                            if ($correct_five_prime ==1) {
                                if ($line[6]=~m/\w(\w+)\s[NACGT]+_\w(\w+)\s[NACGT]+/) {
                                    print $libfa_out ">".$line[0]."_left\n".$five_prime_extension."G".$1."$three_prime_extension\n";
                                    print $libfa_out ">".$line[0]."_right\n".$five_prime_extension."G".$2."$three_prime_extension\n";
                                }else{
                                    $line[6]=~m/\w(\w+)\s\w+/;
                                    print $libfa_out ">".$line[0]."\n".$five_prime_extension."G".$1."$three_prime_extension\n";
                                }
                            }else{
                                if ($line[6]=~m/(\w+)\s[NACGT]+_(\w+)\s[NACGT]+/) {
                                    print $libfa_out ">".$line[0]."_left\n".$five_prime_extension.$1."$three_prime_extension\n";
                                    print $libfa_out ">".$line[0]."_right\n".$five_prime_extension.$2."$three_prime_extension\n";
                                }else{
                                    $line[6]=~m/(\w+)\s\w+/;
                                    print $libfa_out ">".$line[0]."\n"."$five_prime_extension".$1."$three_prime_extension\n";
                                }
                            }
                        }else{
                                if ($line[6]=~m/[NACGT]+\s(\w+)_[NACGT]+\s(\w+)/) {
                                    print $libfa_out ">".$line[0]."_left\n".$five_prime_extension.$1."$three_prime_extension\n";
                                    print $libfa_out ">".$line[0]."_right\n".$five_prime_extension.$2."$three_prime_extension\n";
                                }else{
                                    $line[6]=~m/\w+\s(\w+)/;
                                    print $libfa_out ">".$line[0]."\n"."$five_prime_extension".$1."$three_prime_extension\n";
                                }                            
                        }
						
						$fasta{$line[0]}++;
					}
				}
			}else{
				print $libtab_out $_;
			}
        }
    close ($libtab);
    close ($libfa_out);
    close ($libtab_out);    
    open ($libgff, "<", $temp_dir . "/all_results_together.gff") or die $!;
    open (my $libgff_out, ">", $output_dir.$lib_name.".gff") or die $!;
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
    open ($libtab_out, ">", $output_dir.$lib_name.".large.tab") or die $!;
    open ($libtab, "<", $temp_dir . "/all_results_together.tab") or die $!;
        while (<$libtab>) {
                print $libtab_out $_;
        }
    close $libtab_out;
    close $libtab;
    open ($libgff_out, ">", $output_dir.$lib_name.".large.gff") or die $!;
    open ($libgff, "<", $temp_dir . "/all_results_together.gff") or die $!;
        while (<$libgff>) {
                print $libgff_out $_;
        }
    close $libgff_out;
    close $libgff;
    open (my $coverage_file, ">", $output_dir.$lib_name.".coverage.tab") or die $!;   
        foreach my $key (keys %id_for_lib) {
            print $coverage_file $key."\t".int(scalar(keys %{$id_for_lib{$key}} ))."\n";
        }
    close $coverage_file;
    open (my $mis, ">", $output_dir.$lib_name.".missing.tab") or die $!;   
        foreach my $key (keys %isgone) {
			if (looks_like_number($isgone{$key})) {
                print $mis 	$key." is missing from the library. It was covered by ".
							$isgone{$key}." designs. Maybe it was covered to low or not found in the cld database.\n";
				print $key." is missing from the library. It was covered by ".
							$isgone{$key}." designs. Maybe it was covered to low or not found in the cld database.\n";
            }else{
				print $mis 	$key," ",$isgone{$key};
				print $key," ",$isgone{$key};
			}
		}
		print $mis (scalar keys %genes_from_list) - (scalar keys %id_for_lib ) + (scalar keys %isgone)." genes are missing because of to harsh design criteria or because they were not found by CLD in the data base.\n";
		print " ",((scalar keys %genes_from_list) - (scalar keys %id_for_lib ) + (scalar keys %isgone))," genes are missing because of to harsh design criteria or because they were not found by CLD in the data base.\n";
    close $mis;
	open (my $parameters, ">", $output_dir.$lib_name.".parameters.tab") or die $!;
		foreach my $key (sort keys %something){
			print $parameters $key,"\t",$something{$key},"\n";
		}
	close($parameters);
    if (defined $something{"GUI"}) {
        $trigger_frame->messageBox(
               -icon => 'info',
               -type => 'ok',
               -title => 'Info',
               -message => "Your library $lib_name has been build in $output_dir.",
           );
        #rmdir $temp_dir;
   }
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
      my $match         = $line[21];
      my $querry        = "";
      my $adjust        = 0;
      my $adjustleft    = 0;
      my $start         = 0;
      my $end           = 0;
      my $insert        = 0;
      
      if ($line[24] eq "rc") { #turn and translate the match and querry string if direction is reverse complementary
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
      my @target_name   = split("::", $line[18]);
      my $obj           = $db->get_Seq_by_id($target_name[0]);
      if (!defined $obj) { #have to search in whole chromosom
            $db         = Bio::DB::Fasta->new( $_[1] . '.dna.toplevel.fa', -makeid => \&make_my_id );
            $obj        = $db->get_Seq_by_id($line[18]);
      }
      $start   = $line[19]-$adjust-$_[3]-1;
      if ($start < 0 ) {
            $adjustleft = abs($start);
            $start = 0;
      }
      
      $end     = ($line[20]-$adjust+$_[3]);
      $target  = substr $obj->seq, $start, $end - $start - 1;
      #count insertions and adjust the end property to display
      $insert  = $line[21] =~ tr/I/x/;
      
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
      my @match         = split("-", $line[20]);
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
      if ($line[23] eq "rc") {
            $querry[0]   = reverse_comp($querry[0]);
            $match0 = $match[0];
            $match1 = $match[1];
            #count insertions and deletions to adjust the end/start property
            $deletion0   = $match0 =~ tr/D/o/;
            $deletion1   = $match1 =~ tr/D/o/;
            #calc start and end of the target sequence
            $start       = $line[18] - $_[3] - 1 - $deletion1; # start of the match - the wanted overflow on the left side
            $end         = $line[19] + $line[24] + length($match[0]) + $_[3] - ($leadingbases*2) - ($deletion0*2); # end of the match + spacer size + length of the matchstring + the wanted overflow on the right side - the count of unspec leadingbases*2 (1st and 2nd)
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
            $start       = $line[18] - $line[24] - length($match[0]) - $_[3] + $leadingbases - 1 - $deletion1; # swaped of the if-clause and reverted (-/+), but only the leadingbases of the 2nd string count here
            $end         = $line[19] + $_[3] - $leadingbases - ($deletion0*2); # swaped of the if-clause and reverted (-/+), but only the leadingbases of the 1st string count here
      }
      
      if ($start < 0 ) {
            $adjustleft = abs($start);
            $start = 0;
      }
      
      #count insertions and deletions to adjust the end/start property
      $insert0   = $match0 =~ tr/I/o/;
      $insert1   = $match1 =~ tr/I/o/;
      for (my $i = 0; $i < ($line[24] - $leadingbases*2 - $deletion0 + $deletion1); $i++){
            $spacer .= " ";
      }
      
      #search the target-string in database - substr workaround used, because the Bioperl function does not work proper with large numbers
      my $db             = Bio::DB::Fasta->new( $_[1] . '.all.dna.fa', -makeid => \&make_my_id );
      my @target_name    = split("::", $line[17]);
      my $obj            = $db->get_Seq_by_id($target_name[0]);
      if (!defined $obj) { #have to search in whole chromosom
            $db          = Bio::DB::Fasta->new( $_[1] . '.dna.toplevel.fa', -makeid => \&make_my_id );
            $obj         = $db->get_Seq_by_id($line[17]);
      }
      
      $target  = substr $obj->seq, $start, $end - $start - 1;
      
      #build the table for the popup
      $report .= '<table style="font-size:75%">';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry[0], $_[3], 1, 0, $match[0], 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match[0], $_[3], 1, 1, $match[0], 0, 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Target:</td><td>|'.$start.'*|</td>'.create_popup_string($target, $_[3], 0, 0, $match[0].$spacer.$match[1], 0, 1, $adjustleft).'<td>|'.($end - 1 - $insert1).'*|</td></tr>';
      $report .= '<tr><td>Matchstring:</td><td></td>'.create_popup_string($match[1], $_[3], 1, 1, $match[1], (length($match[0]) + int($line[24]) + $insert0 - ($leadingbases*2) - $deletion0 + $deletion1), 0, $adjustleft).'</tr>';
      $report .= '<tr><td>Query:</td><td></td>'.create_popup_string($querry[1], $_[3], 1, 0, $match[1], (length($match[0]) + int($line[24]) + $insert0 - ($leadingbases*2) - $deletion0 + $deletion1), 0, $adjustleft).'</tr>';
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
      if ( open (my $failfile, "<", $_[4] . "/failfile.tab") && !( $_[2]->{"ignore_missing_id"} ==1)) {
            my $error="";
            while (<$failfile>) {
                  $error.=$_;
            }
            close($failfile);
            die "No Database entry found for \n\"". $error."\" in the \" ".$_[2]->{"ref_organism"}."\" genome.\n Please enter a valid ensembl ID or gene symbol (case sensitive) or change the input option above to FASTA sequence.\n";
      }
      
}
=head2 make_temp_fasta_file_from_coords
#########################################################################################
#name:      make_temp_fasta_file
#function:  creates a temporary fasta file for the bowtie index and builds a trees
#input:     (given id-Array, tree as referemce, something-Hashreference,
#           enzyme db, temp_dir, 1/0 if file or not)
#output:    N/A
#########################################################################################
=cut

sub make_temp_fasta_file_from_coords {
      my %something=%{$_[2]};
      if ( !( $_[2]->{"specific_transcript"} eq "any") && scalar(@{$_[0]}) >1) {
            die "Transcript specificity is only defined for single gene analyses.\n";
      }
	  my @coords;
	  my $header;
	   my $chrom;
      open (my $tempfile, ">", $_[4] . "/tempfile.fasta");
            foreach my $id (@{$_[0]}) { 
                  @coords=split("\t",$id);
                  my $seq_obj = $_[3]->seq($coords[1],$coords[2],$coords[3]); # get a PrimarySeq obj
				  print $coords[1],"\t",$coords[2],"\t",$coords[3],"\n";
                  if ($seq_obj) {
                        $header = "chrom:".$coords[1].":".$coords[2]."..".$coords[3];
                        $chrom = $coords[1];
                        if ( !exists $_[1]->{$chrom} ) {
                              $_[1]->{$chrom} = build_tree( $something{"databasepath"} . $_[2]->{"ref_organism"} . "/" . $chrom . "_indexed" );
                        }
                        print $tempfile ">", $coords[0], " ", $header, "\n", $seq_obj, "\n";
                  } else {
                        open (my $failfile, ">>", $_[4] . "/failfile.tab");
                              print $failfile substr($id,0)."\n";
                        close($failfile);
                  }
            }
      close $tempfile;
      if ( open (my $failfile, "<", $_[4] . "/failfile.tab") && !( $_[2]->{"ignore_missing_id"} ==1)) {
            my $error="";
            while (<$failfile>) {
                  $error.=$_;
            }
            close($failfile);
            die "No Database entry found for \n\"". $error."\" in the \" ".$_[2]->{"ref_organism"}."\" genome.\n Please enter a valid ensembl ID or gene symbol (case sensitive) or change the input option above to FASTA sequence.\n";
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
	  if (!defined $something{"GUI"}) {		
		$something{"input_file"}=$_[1];
		#create a time stamped output folder		
		$temp_dir = scalar localtime();
		$temp_dir =~ s/\s+/\_/ig;
		$temp_dir =~ s/\W+/\_/ig;
		$temp_dir =~ s/[^\w]/_/ig;
		if ($_[2]=~m/\/$/) {
		  $temp_dir=$_[2].$temp_dir;
		}else{
		  $temp_dir=$_[2]."/".$temp_dir;
		}
		print $temp_dir	,"\n";
		mkdir($temp_dir) or die $!;
		#read in the parameter file, must contain each parameter with name and value '=' separated, without any quotes
		open(my $parameterfile, "<", $_[0]) or die $_[0]."could not be opened. No such file or directory.";
			  foreach my $element (<$parameterfile>) {
					chomp $element;
					my @line=split("=",$element);
					$something{$line[0]} = $line[1];
                    #print $line[0],"\t",$something{$line[0]},"\n";
			  }
		close $parameterfile;
      }else{
		$temp_dir = scalar localtime();
		$temp_dir =~ s/\s+/\_/ig;
		$temp_dir =~ s/\W+/\_/ig;
		$temp_dir =~ s/[^\w]/_/ig;
		if ($something{"working_path"}=~m/\/$/) {
		  $temp_dir=$something{"working_path"}.$temp_dir;
		}else{
		  $temp_dir=$something{"working_path"}."/".$temp_dir;
		}
		print $temp_dir	,"\n";
		mkdir($temp_dir) or die $!;
	  }
		my %trees               = ();
		my $seqio_obj           = "";
		my $parallel_number     = 2;
		my $databasepath;
      if ($something{"databasepath"}!~m/\/$/) {
		  $something{"databasepath"}.= "/" ;
		}
	  $databasepath = $something{"databasepath"} . $something{"ref_organism"} . "/" . $something{"ref_organism"};
      #################################################################################################################################################################################
      # upload a file and save it in the temp directory for bowtie
      #################################################################################################################################################################################
      if ($something{"sec_off_target"} ==1) {
			if (-e $something{"databasepath"} .'secondary_off_targets.fasta') {
				system($aligner_path.'bowtie2-build '.$something{"databasepath"}.'secondary_off_targets.fasta '.$temp_dir.'/temp_sec ;');
			}else{
				die $something{"databasepath"} .'secondary_off_targets.fasta'."could not be opened. No such file or directory.";
			}
      }
      #################################################################################################################################################################################
      # For ENSEMBLE: define the path to the bowtie index and do checks if an database entry is found for the ensemble accesion number - if all checks passed, create the $seqio_obj
      #################################################################################################################################################################################
      if ( $something{"data_type"} eq "ensemble_acc" || $something{"data_type"} eq "gene_symbol") {
            my $db = Bio::DB::Fasta->new( $databasepath . '.all.dna.fa', -makeid => \&make_my_id );
            my @ids = ();
			if (defined $something{"GUI"}) {
                @ids=split("\n",$something{"genes"})
            }else{
				open (my $infile, "<", $something{"input_file"});
                  while (<$infile>) {
                        my $line = $_;
                        chomp $line;
                        push @ids, $line;
                  }
				close $infile;
			}
            my %seen;
            @id = sort @id;
            foreach my $string (@ids) {            
                next unless $seen{$string}++;
                 print $string." is a duplicated ID.\n Please ensure that all Sequence/Gene/Coordiante IDs are unique.\n" ;
                 if (defined $something{"GUI"}) {$mw->update;};
                die $string." is a duplicated ID.\n Please ensure that all Sequence/Gene/Coordiante IDs are unique.\n" ;
            }            
            make_temp_fasta_file(\@ids, \%trees, \%something, $db, $temp_dir, 1);
            $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ); #read the temporary fasta file            
      } elsif($something{"data_type"} eq "fasta") { 
            #################################################################################################################################################################################
            # For FASTA: define the path to the bowtie index and do checks if the fasta sequence/file is in the right format - if all checks passed, create the $seqio_obj
            #################################################################################################################################################################################
            my $count=0;
            my $temp="";
            my %seen;
            open(my $infile, "<",$something{"input_file"});
                        while (my $line = <$infile>){
                              if ($line=~m/^(>.+)/) {
                                    $count++;
                                    next unless $seen{$1}++;
                                     print $string." is a duplicated ID.\n Please ensure that all Sequence/Gene/Coordiante IDs are unique.\n" ;
                                     if (defined $something{"GUI"}) {$mw->update;};
                                    die $string." is a duplicated ID.\n Please ensure that all Sequence/Gene/Coordiante IDs are unique.\n" ;                  
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
      }elsif($something{"data_type"} eq "coordinates") { 
            my $db = Bio::DB::Fasta->new( $databasepath . '.dna.toplevel.fa', -makeid => \&make_my_id );
            my @ids = ();
			if (defined $something{"GUI"}) {
                @ids=split("\n",$something{"genes"});
				print join("_",@ids),"\n";
            }else{
				open (my $infile, "<", $something{"input_file"});
                  while (<$infile>) {
                        my $line = $_;
                        chomp $line;
                        push @ids, $line;
                  }
				close $infile;
			}
            my %seen;
            @id = sort @id;
            foreach my $string (@id) {            
                next unless $seen{$string}++;
                print $string." is a duplicated ID.\n Please ensure that all Sequence/Gene/Coordiante IDs are unique.\n" ;
                if (defined $something{"GUI"}) {$mw->update;};
                die $string." is a duplicated ID.\n Please ensure that all Sequence/Gene/Coordiante IDs are unique.\n" ;
            }
            make_temp_fasta_file_from_coords(\@ids, \%trees, \%something, $db, $temp_dir, 1);
            $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ); #read the temporary fasta file
		}
      #################################################################################################################################################################################
      # Start the creation of the report 
      #################################################################################################################################################################################
      #define empty object to be filled in the process
      my $fname               = "";
      my $dont_asses_context  = 0;
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
            #$statistics{$fname}{"Number of successful designs"}   = 0;
            if ( $chrom eq "" ) { $chrom = $fname; } #if $chrom is still empty fill it with the sequence' id
            my $whole_seq = $seq_obj->seq; #deduce the complete nucleotide sequence as a alphanumerical string from the SeqIo object
            if ($statistics{$fname}{"seq_length"}!=0) {
                  if (int(log($statistics{$fname}{"seq_length"})) <= 8) {
                        if (int(log($statistics{$fname}{"seq_length"}))>=2) {
                              $parallel_number=int(log($statistics{$fname}{"seq_length"}));
                        }else{
                             $parallel_number=8;
                        }
                  }else{
                        $parallel_number=$max_parallel;
                  }
            }else{
                  $parallel_number=2;
            }
            if ($parallel_number<2) {
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
            if (defined $something{"GUI"}) {$mw->update;};
      } # end the loop after Hash and statistic were created, so the Bowtie Magic will be executed only once
      %statistics    = ();
      %CRISPR_hash   = ();
      foreach my $fname ( @fname_array ) {
            my $json = read_file( $temp_dir . "/" .$fname. '.json', { binmode => ':raw' } );
            %CRISPR_hash = ( %CRISPR_hash, %{ decode_json $json } );
            unlink $temp_dir . "/" . $fname. ".json";
            $json = read_file( $temp_dir . "/" . $fname. 'stats.json', { binmode => ':raw' } );
            %statistics=( %statistics, %{ decode_json $json } );
            unlink $temp_dir . "/" . $fname. "stats.json";
      }
	  if ($something{"PAM"} eq "any") {
			if ($something{"textpam"}=~m/([^ACGTUKMSWRYBVDHN]+)/g) {
				  print "The PAM you entered \: \" ".$something{"textpam"}." \"must contain only IUPAC code ACGTUKMSWRYBVDHN";
				}else{
				  $something{"PAM"}=$something{"textpam"} ;
			}        
		
	  }
      ###########################################################################################################################################################################
      #Bowtie for single sequence
      ###########################################################################################################################################################################
                        if ($something{"kind"} eq "single") {
                               open (my $crisprs, ">", $temp_dir . "/temp_CRISPRS.fasta");
                                    foreach my $seq ( sort( keys(%CRISPR_hash) ) ) {
                                          foreach my $key ( sort( keys(%{$CRISPR_hash{$seq}}) ) ) {
                                                my $counter=0;
                                                foreach my $pam (from_pam_to_fasta_combis($something{"PAM"})){
                                                           if ($something{"bowtie_version"} eq "blast") {
                                                                print $crisprs "\>" . $key ."_$counter". "\n";
                                                           }else{
                                                                print $crisprs "\>" . $key . "\n";
                                                           }
                                                           if(${ ${ $CRISPR_hash{$seq} } {$key} }{"strand"} eq "minus"){
                                                                if ($something{"PAM_location"} eq "5_prime") {
                                                                    print $crisprs    $pam.substr(  reverse_comp(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}),
                                                                                                    length($pam),
                                                                                                    (length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-length($pam)-$something{"unspecific_leading_bases"})
                                                                                                    )
                                                                                    . "\n"; #$whole_CRISPR_seq
                                                                }else{
                                                                    print $crisprs    substr(     reverse_comp(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"}),
                                                                                                $something{"unspecific_leading_bases"},
                                                                                                (length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-length($pam)-$something{"unspecific_leading_bases"}))
                                                                                    .$pam
                                                                                    . "\n"; #$whole_CRISPR_seq
                                                                }
                                                                
                                                                  
                                                            }else{
                                                             if ($something{"PAM_location"} eq "5_prime") {
                                                                 print $crisprs $pam.substr(        ${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"},
                                                                                                    length($pam),
                                                                                                    (length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-length($pam)-$something{"unspecific_leading_bases"})
                                                                                        )
                                                                                                      . "\n"; #$whole_CRISPR_seq
                                                             }else{
                                                                print $crisprs substr(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"},
                                                                                       $something{"unspecific_leading_bases"},
                                                                                       (length(${ ${ $CRISPR_hash{$seq} } {$key} }{"nucseq"})-length($pam)-$something{"unspecific_leading_bases"}))
                                                                                                      .$pam. "\n"; #$whole_CRISPR_seq
                                                             }
                                                            }
                                                $counter++;
                                                }
                                          }
                                    }
                              close $crisprs;
                              my $k = 30;
                              if ($something{"targets-allowed"}>30) {
                                $k=$something{"targets-allowed"}+1;
                              }                              
                              #####################################################################################################################################################################
                              #teemp_sec
                              if ($something{"bowtie_version"} eq "bowtie2") {
                                    if ($something{"offtargetdb"} eq "gDNA") {
                                          system( $aligner_path.'bowtie2 -p '.$parallel_number.' -f -k '.$k.' --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath .".dna". ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }elsif($something{"offtargetdb"} eq "cDNA"){
                                          system( $aligner_path.'bowtie2 -p '.$parallel_number.' -f -k '.$k.' --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath.".cdna" . ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }else{
                                          system( $aligner_path.'bowtie2 -p '.$parallel_number.' -f -k '.$k.' --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $databasepath.".genome" . ' -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    }
                              }elsif($something{"bowtie_version"} eq "bowtie"){
                                    if ($something{"offtargetdb"} eq "gDNA") {
                                          system( $aligner_path.'bowtie ' . $databasepath .".dna". ' ' . $temp_dir . "/" . 'temp_CRISPRS.fasta -f -v 3 -y -k '.$k.' -S --sam-nohead --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir . '/temp_out.bwt' );
                                    }elsif($something{"offtargetdb"} eq "cDNA"){
                                          system($aligner_path. 'bowtie ' . $databasepath .".cdna". ' ' . $temp_dir . "/" . 'temp_CRISPRS.fasta -f -v 3 -y -k '.$k.' -S --sam-nohead --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir . '/temp_out.bwt' );
                                    }else{
                                          if (-e $databasepath.".genome.1.ebwtl") {
                                                system( $aligner_path.'bowtie ' . $databasepath.".genome" . ' ' . $temp_dir . "/" .'temp_CRISPRS.fasta -f -v 3 -y -k '.$k.' -S --sam-nohead --large-index  --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir .'/temp_out.bwt' );
                                          }else{
                                                system( $aligner_path.'bowtie ' . $databasepath.".genome" . ' ' . $temp_dir . "/" .'temp_CRISPRS.fasta -f -v 3 -y -k '.$k.' -S --sam-nohead --sam-nosq -p '.$parallel_number.'  > ' . $temp_dir .'/temp_out.bwt' );
                                          }
                                    }
                              }else{
                                my $length=$something{"min_length"}-$something{"unspecific_leading_bases"}+length($something{"PAM"});
                                 if ($something{"offtargetdb"} eq "gDNA") {
                                        if (!-e  $databasepath.".blast.gDNA.nsq") {
                                            system($aligner_path."makeblastdb -in ".$databasepath.".all.dna.fa -dbtype nucl -title ".$something{"ref_organism"}.".blast.gDNA -out ".$databasepath.".blast.gDNA -parse_seqids");
                                        }
                                        system( $aligner_path.'blastn -db '.$databasepath.'.blast.gDNA -query '. $temp_dir . "/" . 'temp_CRISPRS.fasta -task blastn-short -outfmt 15 -parse_deflines -num_threads '.$parallel_number.' | sed "s/lcl|//g" | grep "AS:i:'.$length.'"  > ' . $temp_dir .'/temp_out.bwt' );
                                    }elsif($something{"offtargetdb"} eq "cDNA"){
                                        if (!-e  $databasepath.".blast.cDNA.nsq") {
                                                system($aligner_path."makeblastdb -in ".$databasepath.".cdna.all.fa -dbtype nucl -title ".$something{"ref_organism"}.".blast.cDNA -out ".$databasepath.".blast.cDNA -parse_seqids");
                                        }
                                         system( $aligner_path.'blastn -db '.$databasepath.'.blast.cDNA -query '. $temp_dir . "/" . 'temp_CRISPRS.fasta -task blastn-short -outfmt 15 -parse_deflines -num_threads '.$parallel_number.' | sed "s/lcl|//g" | grep "AS:i:'.$length.'"  > ' . $temp_dir .'/temp_out.bwt' );
                                    }else{
                                        if (!-e  $databasepath.".blast.genome.nsq") {
                                                system($aligner_path."makeblastdb -in ".$databasepath.".dna.toplevel.fa -dbtype nucl -title ".$something{"ref_organism"}.".blast.genome -out ".$databasepath.".blast.genome -parse_seqids");
                                         }
                                        system( $aligner_path.'blastn -db '.$databasepath.'.blast.genome -query '. $temp_dir . "/" . 'temp_CRISPRS.fasta -task blastn-short -outfmt 15 -parse_deflines -num_threads '.$parallel_number.' | sed "s/lcl|//g" | grep "AS:i:'.$length.'"  > ' . $temp_dir .'/temp_out.bwt' );
                                          
                                    }
                                
                              }
                              my $id="";
                                my $seq = "";
                              open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                                    while (my $line = <$bowtie>) {
                                          chomp $line;
                                          my @line = split( "\t", $line );
                                          
                                          if ($line[2] eq "*") {
                                                if($line[0] =~m/^(\S+?)_(\d+)_(\d+)$/){
                                                   $id=$1."_".$2."_".$3;
                                                    $seq=$1;
                                               }
                                                ${ ${ $CRISPR_hash{$seq} } {$id} }{"hits"}.=";;".$line[2]."§§0§§0§§NA§§0§§NA";
                                          }else{
                                            if($line[0] =~m/^(\S+?)_(\d+)_(\d+)$/){
                                                $id=$1."_".$2."_".$3;
                                                 $seq=$1;
                                            }
                                                $line=~m/NM:i:(\d+)/;
                                                my $edit_distance=$1;
                                                if (($edit_distance <= $something{"edit_distance_allowed"})) {
                                                      my @line = split( "\t", $line );
                                                      #decide if it is forward (fw) or backward (bw) query-sequence
                                                      my $direction = "rc";
                                                      if ($line[1] == 0 || $line[1] == 256) {
                                                            $direction = "fw";
                                                      }
                                                      if ( $line[0] =~ m/(\S+?_[^_]+_[^_]+)$/ig ) {
                                                            my @matchstringo=make_mismatch_string (\$line,$something{"unspecific_leading_bases"}, $direction, $something{"PAM_location"});
                                                            my $cond=0;
                                                            if ($something{"PAM_location"} eq "3_prime") {
                                                                if ( (  $direction eq "fw"
                                                                    && !(substr(join("",@matchstringo),scalar(@matchstringo)-length($something{"PAM"}))=~m/X/)
                                                                    )
                                                                || ($direction eq "rc"
                                                                    && !(substr(join("",@matchstringo),0,length($something{"PAM"}))=~m/X/)
                                                                    )
                                                                ) {
                                                                    $cond=1;
                                                                }
                                                            }else{
                                                                if ( (  $direction eq "fw"
                                                                    && !(substr(join("",@matchstringo),0,length($something{"PAM"}))=~m/X/)
                                                                    )
                                                                || ($direction eq "rc"
                                                                    && !(substr(join("",@matchstringo),scalar(@matchstringo)-length($something{"PAM"}))=~m/X/)
                                                                    )
                                                                ) {
                                                                    $cond=1;
                                                                }
                                                            }
                                                            
                                                            if ( $cond==1 ) {                                                                
                                                                  my $startcoordinate=0;
                                                                  
                                                                  if ($direction eq "fw"){
                                                                    $startcoordinate=$something{"unspecific_leading_bases"};
                                                                  }                                                                 

                                                                if ($something{"offtargetdb"} eq "genomeDNA") {
                                                                        my $namestuff="";
                                                                        if ( !exists $trees{$line[2]} ) { #TODO if und else fast identisch, aber relativ kurzer Part und daher Funktion performancetechnisch nachteilig
                                                                              $trees{$line[2]} = build_tree( $something{"databasepath"}  . $something{"ref_organism"} . "/" . $line[2] . "_indexed" );
                                                                        }
                                                                        my $annotations = $trees{$line[2]}->fetch( int($line[3]), int(($line[3])) );
                                                                        foreach  my $anno ( @{$annotations} ) {
                                                                              if ( $anno =~ m/gene_(\S+?)_([0-9]+)_([0-9]+)/ ) {
                                                                                    $namestuff=$1;
                                                                              }
                                                                        }
                                                                    if ($namestuff eq "" && $something{"ignore_intergenic"} ==1) {                                                                              
                                                                    }elsif($namestuff ne ""){
                                                                              ${ ${ $CRISPR_hash{$seq} } {$id} }{"hits"}.=
                                                                                ";;".
                                                                                    $namestuff.
                                                                                    "§§".($line[3]-$startcoordinate).
                                                                                    "§§".($line[3]+@matchstringo-$startcoordinate).
                                                                                    "§§".join("",@matchstringo).
                                                                                    "§§".$edit_distance.
                                                                                    "§§".$direction.
                                                                                    "§§".$line[2];
                                                                    }elsif($namestuff eq "" && $something{"ignore_intergenic"} eq 0){
                                                                              ${ ${ $CRISPR_hash{$seq} } {$id} }{"hits"}.=
                                                                                ";;".$line[2].
                                                                                    "§§".($line[3]-$startcoordinate).
                                                                                    "§§".($line[3]+@matchstringo-$startcoordinate).
                                                                                    "§§".join("",@matchstringo).
                                                                                    "§§".$edit_distance.
                                                                                    "§§".$direction.
                                                                                    "§§".$line[2];
                                                                    }
                                                                }else{
                                                                        ${ ${ $CRISPR_hash{$seq} } {$id} }{"hits"}.=
                                                                                ";;".$line[2].
                                                                                    "§§".($line[3]).
                                                                                    "§§".($line[3]+@matchstringo).
                                                                                    "§§".join("",@matchstringo).
                                                                                    "§§".$edit_distance.
                                                                                    "§§".$direction.
                                                                                    "§§".$line[2];
                                                                }
                                                            }
                                                      }    
                                                }
                                          }
                                    }
                              close $bowtie;
                              
                              unlink $temp_dir . "/temp_out.bwt";
                              if ( $something{"sec_off_target"} ==1 ) { #ckeck if this is wanted
                                    ###############################################################################################################################################################
                                    
                                    #do send a bowtie2 job fot the two temporary written fasta files as if they were paired seqencing reads and save the result in a ~out.bwt file
                                    system( 'bowtie2 -p '.$parallel_number.'  -f -k '.$k.' --'.$something{"bowtie_mode"}.' --end-to-end --no-hd -x ' . $temp_dir .'/temp_sec -U ' . $temp_dir . '/temp_CRISPRS.fasta > ' . $temp_dir . '/temp_out.bwt' );
                                    open (my $bowtie, "<", $temp_dir . "/temp_out.bwt");
                                          my $edit_distance=0;
                                          while (my $line = <$bowtie>) {
                                                chomp $line;
                                                my @line = split( "\t", $line );
                                                if($line=~m/NM:i:(\d+)/){
                                                      $edit_distance=$1;
                                                }
                                                if (($edit_distance <= $something{"edit_distance_allowed"}) && ($line[2] ne "*")) {
                                                      if ( $line[0] =~ m/(\S+?)_([^_]+)_([^_]+)$/ ) {
                                                            my $key = $1."_".$2."_".$3;
                                                            my $seq = $1;
                                                            push @{ ${ ${ $CRISPR_hash{$seq} } {$key} }{"sec_hits"} }, $line[2];
                                                      }
                                                }
                                          }
                                    close $bowtie;
                                 unlink $temp_dir . "/temp_out.bwt";
                              }
                              unlink $temp_dir . "/temp_CRISPRS.fasta";
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
                                    system( 'bowtie2 -p '.$parallel_number.' -f -k '.$k.' --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath .".dna". ' -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }elsif($something{"offtargetdb"} eq "cDNA"){
                                    system( 'bowtie2 -p '.$parallel_number.' -f -k '.$k.' --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath.".cdna". ' -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
                              }else{
                                    system( 'bowtie2 -p '.$parallel_number.' -f -k '.$k.' --'.$something{"bowtie_mode"}.' --end-to-end --no-mixed --no-discordant --no-hd -x ' . $databasepath.".genome".'  -I '.$something{"minspacerlength"}.' -X 100 -1 ' . $temp_dir . '/temp_LEFTCRISPRS.fasta'. ' -2 ' . $temp_dir . '/temp_RIGHTCRISPRS.fasta'.' > ' . $temp_dir . '/temp_out.bwt' );
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
                                                if ( $line[0] =~ m/(\S+?_[^_]+_[^_]+)/ig ) {
                                                      $line[0] =~m/(\S+)_(\S+)_/;
                                                      my $seq = $1;
                                                      my @matchstringo=make_mismatch_string (\$line,$something{"unspecific_leading_bases"}, $direction,$something{"PAM_location"});
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
																			  if ($namestuff eq "" && $something{"ignore_intergenic"}==1) {                                                                              
																				}elsif($namestuff ne ""){
																					  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=
                                                                                      ";;".
                                                                                        $namestuff.
                                                                                        "§§".($line[3]-$startcoordinate).
                                                                                        "§§".($line[3]+@matchstringo-$startcoordinate).
                                                                                        "§§".join("",@matchstringo).
                                                                                        "§§".$edit_distance.
                                                                                        "§§".$direction.
                                                                                        "§§".$spacer.
                                                                                        "§§".$line[2];
																				}elsif($namestuff eq "" && $something{"ignore_intergenic"}==0){
																					  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=
                                                                                      ";;".$line[2].
                                                                                        "§§".($line[3]-$startcoordinate).
                                                                                        "§§".($line[3]+@matchstringo-$startcoordinate).
                                                                                        "§§".join("",@matchstringo).
                                                                                        "§§".$edit_distance.
                                                                                        "§§".$direction.
                                                                                        "§§".$spacer.
                                                                                        "§§".$line[2];
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
                                                                                    my @lasthitarray=split("§§",$hitarray[-1]);
                                                                                    $lasthitarray[3].="-".join("",@matchstringo);
                                                                                    $hitarray[-1]=join("§§",@lasthitarray);
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
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=";;".$line[2]."§§".($line[3]-$startcoordinate)."§§".($line[3]+@matchstringo-$startcoordinate)."§§".join("",@matchstringo)."§§".$edit_distance."§§".$direction."§§".$spacer;
                                                            }elsif($line[1]==147){
                                                                  my @hitarray=split(";;",${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"});
                                                                  my @lasthitarray=split("§§",$hitarray[-1]);
                                                                  $lasthitarray[3].="-".join("",@matchstringo);
                                                                  $hitarray[-1]=join("§§",@lasthitarray);
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}=join(";;",@hitarray);
                                                            }else{
                                                                  ${ ${ $CRISPR_hash{$seq} } {$line[0]} }{"hits"}.=
                                                                    ";;".
                                                                        $line[2].
                                                                        "§§".($line[3]-$startcoordinate).
                                                                        "§§".($line[3]+@matchstringo-$startcoordinate).
                                                                        "§§".join("",@matchstringo).
                                                                        "§§".$edit_distance.
                                                                        "§§".$direction.
                                                                        "§§".$spacer.
                                                                        "§§".$line[2];
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
                  #open(my $debug_log, ">", "debug_log.txt") or die $!;
                  my $number_of_hits = 1;
                  SEQLOOP: foreach my $fname ( keys(%CRISPR_hash) ) {
                        CRISPRHASHLOOP: foreach my $key ( keys(%{$CRISPR_hash{$fname}}) ) {
                              if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} ) {
                                    $number_of_hits = scalar(split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"}));
                                    ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"} = $number_of_hits-1;
                                    if (( $number_of_hits > $something{"targets-allowed"} + 2 || $number_of_hits < 1 ) ) {
                                          $statistics{$fname}{"Number of designs excluded because they hit multiple targets or none"}++;
                                          delete $CRISPR_hash{$fname}{$key};
                                          next CRISPRHASHLOOP;
                                    }else{
                                          $statistics{$fname}{"Number of designs that hit a specific target"}++;
                                    }
                              }else{
                                    delete $CRISPR_hash{$fname}{$key};
                                    next CRISPRHASHLOOP;
                              }
                              
                              
                              if ($something{"sec_off_target"} ==1 ){
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
                              
                              my $whole_crisp_seq = "A";
                              
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
                  #close($debug_log);
                  #################################################################################################################################################################################
                  #reopen the main loop which is looping through the sequence found in the SeqIO object, may be one or many more
                  #################################################################################################################################################################################
                  $seqio_obj = Bio::SeqIO->new( -file => $temp_dir . "/tempfile.fasta", -format => "fasta" ); #read the temporary fasta file
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
                              foreach my $hit (@targets){
                                    if ($hit ne "") {
                                          #print $hit."\n";
                                          my @splithit=split("§§",$hit);
                                          if (${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}-((20-(100/${ ${ $CRISPR_hash{$fname} } {$key} }{"length"}*$splithit[4])))>0) {
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}-((20-(100/${ ${ $CRISPR_hash{$fname} } {$key} }{"length"}*$splithit[4])));
                                          }else{
                                                ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=0;
                                          }
                                    }
                              }
                              if(exists(${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"})){
                                    @{${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"}}[0]=${ ${ $CRISPR_hash{$fname} } {$key} }{"score"};
                                    ${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}=${ ${ ${ $CRISPR_hash{$fname} } {$key} } {"context"} }{"new_score"};
                              }
                              foreach my $i (0..4){
                                    @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[$i]=round_digits(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[$i],4);
                              }
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"spec_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[0];
                              ${ ${ $CRISPR_hash{$fname} } {$key} }{"anno_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[1];
							  ${ ${ $CRISPR_hash{$fname} } {$key} }{"custom_score"}= @{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}}[2];
                        }
                       
                        #####################################################################################################################################################################
                        # Create the Table for the findings
                        #####################################################################################################################################################################
                              
                              open (my $outfiletab, ">", $temp_dir . "/" . $fname . "_" . "table.tab");
                                    if ($something{"kind"} eq "single") {
                                          print $outfiletab "Name\tLength\tChromosome\tStart\tEnd\tStrand\tNucleotide sequence\tGene Name\tTranscripts\tTranscript:: Exon\tNumber of Cpg Islands hit\tSequence around the cutside\t%A %C %T %G\tS-Score\tA-Score\tCustom-Score\tDoench-Score\tXu-Score\tpercent of total transcripts hit\tMatch-Target\tMatch-Chromosome\tMatch-Start\tMatch-End\tMatchstring\tEditdistance\tNumber of Hits\tDirection\tStart_rti\tEnd_rti\n";
                                    }
                                    else {
                                          print $outfiletab "Name\tLength\tStart\tEnd\tStrand\tNucleotide sequence\tGene Name\tTranscripts\tTranscript:: Exon\tNumber of Cpg Islands hit\tSequence around the cutside\t%A %C %T %G\tS-Score\tA-Score\tCustom-Score\tpercent of total transcripts hit\tTarget\tMatch-start\tMatch-end\tMatchstring\tEditdistance\tNumber of Hits\tDirection\tSpacer\tStart_rti\tEnd_rti\n";
                                    }
									if ($something{"purpose"} eq "non-coding") {
											PRINTLOOP: foreach my $key (
																	sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} }																	
																 keys(%{$CRISPR_hash{$fname}}) ) {
											  $statistics{$fname}{"Number of successful designs"}++;
											  my @targets=split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} );
											  #write the tab-delimited file
											  HITLOOP: foreach my $hit (@targets){
													if ($hit ne "") {
														  #print the candidates name
														  print $outfiletab "$key\t";
														  my @splithit = split("§§",$hit);
														  #print its length on this whole sequence these are genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} . "\t";
														  }
														  #print its start on this whole sequence these are  genomic coordinates
														  my @locus=split("::",$statistics{$fname}{"seq_location"});
                                                          if ( exists $statistics{$fname}{"seq_location"}) {
																print $outfiletab $locus[0] . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}+$locus[1] . "\t";
														  }
														  #print its end on this whole sequence these are  genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}+$locus[1] . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"} ) {
																if ($something{"kind"} eq "single") {
																	  if($something{"PAM_location"} eq "3_prime"){
                                                                            if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                                print $outfiletab reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length($something{"PAM"}))) . "\t";
                                                                            }else{
                                                                                  print $outfiletab substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"})) . "\t";
                                                                            }
                                                                        }else{
                                                                            if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                                print $outfiletab reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"}))) . "\t";
                                                                            }else{
                                                                                print $outfiletab substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length($something{"PAM"})) . "\t";
                                                                            }
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
														  print $outfiletab $splithit[0]."\t".$splithit[-1]."\t".$splithit[1]."\t".$splithit[2]."\t".$splithit[3]."\t".$splithit[4]."\t";
														  print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"}, "\t";
														  print $outfiletab $splithit[5]."\t";
														  if (!($something{"kind"} eq "single")) {
																print $outfiletab $splithit[6]."\t";
														  }
                                                          if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}."\t";
														  }
														  #print its end on this whole sequence these are  genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}."\t";
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
									}else{									
										PRINTLOOP: foreach my $key (
																	sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} }
																	sort { $CRISPR_hash{$fname}{$b}->{"anno_score"} <=> $CRISPR_hash{$fname}{$a}->{"anno_score"} }
																	sort { if($something{"sort_by_rank"}==1){$CRISPR_hash{$fname}{$b}->{"custom_score"} <=> $CRISPR_hash{$fname}{$a}->{"custom_score"} }else{return 1} }
																 keys(%{$CRISPR_hash{$fname}}) ) {
											  $statistics{$fname}{"Number of successful designs"}++;
											  my @targets=split(";;",${ ${ $CRISPR_hash{$fname} } {$key} }{"hits"} );
											  #write the tab-delimited file
											  HITLOOP: foreach my $hit (@targets){
													if ($hit ne "") {
														  #print the candidates name
														  print $outfiletab "$key\t";
														  my @splithit = split("§§",$hit);
														  #print its length on this whole sequence these are genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"length"} . "\t";
														  }
														  #print its start on this whole sequence these are  genomic coordinates
														  my @locus=split("::",$statistics{$fname}{"seq_location"});
                                                          if ( exists $statistics{$fname}{"seq_location"}) {
																print $outfiletab $locus[0] . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}+$locus[1] . "\t";
														  }
														  #print its end on this whole sequence these are  genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}+$locus[1] . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} ) {
																print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} . "\t";
														  }
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"} ) {
																if ($something{"kind"} eq "single") {
																	  if($something{"PAM_location"} eq "3_prime"){
                                                                            if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                                print $outfiletab reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length($something{"PAM"}))) . "\t";
                                                                            }else{
                                                                                  print $outfiletab substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"})) . "\t";
                                                                            }
                                                                        }else{
                                                                            if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                                print $outfiletab reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"}))) . "\t";
                                                                            }else{
                                                                                print $outfiletab substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length($something{"PAM"}))." ".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length($something{"PAM"})) . "\t";
                                                                            }
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
														  print $outfiletab $splithit[0]."\t".$splithit[-1]."\t".$splithit[1]."\t".$splithit[2]."\t".$splithit[3]."\t".$splithit[4]."\t";
														  print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"number_of_hits"}, "\t";
														  print $outfiletab $splithit[5]."\t";
														  if (!($something{"kind"} eq "single")) {
																print $outfiletab $splithit[6]."\t";
														  }
                                                          if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}."\t";
														  }
														  #print its end on this whole sequence these are  genomic coordinates
														  if ( exists ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"} ) {
																#print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1] . "\t";
                                                                print $outfiletab ${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}."\t";
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
                              
                              if ($something{"out_gff"} ==1) { 
                                    open(my $gfffile, ">",$temp_dir . "/" . $fname . ".gff" ) or die $!;
                                     print $gfffile "##gff-version 3\n";
                                          PRINTLOOP: foreach my $key (
																	  sort { $CRISPR_hash{$fname}{$b}->{"spec_score"} <=> $CRISPR_hash{$fname}{$a}->{"spec_score"} }
																	  sort { $CRISPR_hash{$fname}{$b}->{"anno_score"} <=> $CRISPR_hash{$fname}{$a}->{"anno_score"} }
																	  sort { if($something{"sort_by_rank"}==1){$CRISPR_hash{$fname}{$b}->{"custom_score"} <=> $CRISPR_hash{$fname}{$a}->{"custom_score"} }else{return 1} }
																	keys(%{$CRISPR_hash{$fname}}) ) {
                                                my @locus=split("::",$statistics{$fname}{"seq_location"});
                                               # print $gfffile $locus[0]."\tcld\tCRISPRtarget\t".(${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}-500+$locus[1])."\t";
                                            #print $gfffile (${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}-500+$locus[1])."\t".sum(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}})."\t";
                                             print $gfffile $locus[0]."\tcld\tCRISPRtarget\t".(${ ${ $CRISPR_hash{$fname} } {$key} }{"start"}+$locus[1])."\t";
                                            print $gfffile (${ ${ $CRISPR_hash{$fname} } {$key} }{"end"}+$locus[1])."\t".sum(@{${ ${ $CRISPR_hash{$fname} } {$key} }{"score"}})."\t";
                                                if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                      print $gfffile "-"."\t."."\t";
                                                }else{
                                                      print $gfffile "+"."\t."."\t";
                                                }
                                                print $gfffile "id=".$key."; ";
                                                print $gfffile "spec_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"spec_score"}."; ";
                                                print $gfffile "anno_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"anno_score"}."; ";
												 print $gfffile "eff_score=".${ ${ $CRISPR_hash{$fname} } {$key} }{"custom_score"}."; ";
												#die "correct the PAM in the output gff\nALSO CORRECT IN E_CRISP V5\nAdd the option to sort by new DOENCH score\nadd option to also sort by DOENCH off-target score\nStore all other options in something\nManage GUI option\n";
                                                if ($something{"kind"} eq "single") {
                                                               if($something{"PAM_location"} eq "5_prime"){
                                                                       
                                                                        if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                                    print $gfffile "seq=".$something{"PAM"}."_".reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"})))."; ";
                                                                        }else{
                                                                                    print $gfffile "seq=".$something{"PAM"}."_".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length($something{"PAM"})). "; ";
                                                                        }
                                                                    }else{
                                                                        if (${ ${ $CRISPR_hash{$fname} } {$key} }{"strand"} eq "minus") {
                                                                            print $gfffile "seq=".reverse_comp(substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},length($something{"PAM"})))."_".$something{"PAM"}. "; ";
                                                                        }else{
                                                                            print $gfffile "seq=".substr(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"},0,length(${ ${ $CRISPR_hash{$fname} } {$key} }{"nucseq"})-length($something{"PAM"}))."_".$something{"PAM"}. "; ";
                                                                        }
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
                            
                              #unlink $temp_dir . "/tempfile.fasta";
                              my $zip    = Archive::Zip->new();
                              my $member = "";
                              if ( $something{"out_gff"} ==1) { $member = $zip->addFile( $temp_dir . "/" . $fname . ".gff", $fname . "_CRISPR.gff" ); }
                              $member = $zip->addFile( $temp_dir . "/" . $fname . "_" . "table.tab", $fname . "_CRISPR.tab" );
                              $zip->writeToFileNamed( $temp_dir . '/' . $fname . '.zip' );
                              
                    #####################################################################################################################################################################
                    # Print the report site with the results (header was created earlier in the loop)
                    #####################################################################################################################################################################
                    $CRISPR_cnt{$fname}++;
                  print "$fname is completed 100%\n";
                  if (defined $something{"GUI"}) {$mw->update;};
                  } #end Sequence loop
				  my %all_stats;
                  open( my $missing_log, ">", "missing_log.txt") or die $!;
                    foreach my $key (keys %statistics){
                        if( $statistics{$key} =~m/[a-zA-Z]+/){
                            if ($statistics{$key}{"Number of successful designs"}==0) {
                                print $missing_log 'Query name: '.$statistics{$key}{"seq_name"}.'   Query length: '.$statistics{$key}{"seq_length"}.'   Query location: '.$statistics{$key}{"seq_location"}."\n";
                            
                            }
                            foreach my $subkey (sort keys(%{$statistics{$key}})){
                                if ($subkey =~ m/Number/) {
                                    if ($statistics{$key}{"Number of successful designs"}==0) {
                                        print $missing_log $key."\t".$subkey.' = '.$statistics{$key}{$subkey}."\n";
                                    }
                                    $all_stats{$subkey}+=$statistics{$key}{$subkey};
                                }
                            }
                        }
                    }
                    close($missing_log);
				   foreach my $key (sort keys %all_stats){
						print $key.' = '.$all_stats{$key}."\n";
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
                                   }elsif($filename=~m/\.zip/){
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
	  my $gene_id					= $seq_obj->display_id;
      my $whole_seq                 = $seq_obj->seq;
      my $count                     = 0;
      my %finished_CRISPR_hash      = ();
      my $pm                        = Parallel::ForkManager->new($parallel_number);
      my $cutnumber                 = int( length($whole_seq) / int( $parallel_number - 1 ) );
      my $cut                       = 0;
      my @cuts                      = ();
      my $correction=500;
        if ($something{"data_type"} eq "coordinates") {
            $correction=0;
        }
        
      my %tempstatistics            = ();
      my $input5=$something{"preceding"};
      my $input3=$something{"PAM"};
      if ($something{"PAM"} eq "any") {
            if ($something{"textpam"}=~m/([^ACGTUKMSWRYBVDHN]+)/g) {
                  print "The PAM you entered \: \" ".$something{"textpam"}." \" must contain only IUPAC code ACGTUKMSWRYBVDHN\n But it contais \"$1\"" ;
                }else{
                  $input3=$something{"PAM"}=$something{"textpam"} ;
            }        
        
      }
      
      my $minlength=$something{"min_length"}-1;
      my $maxlength=$something{"max_length"}-1;
      
      my $input5_rev=reverse $input5;
      my $input3_rev=reverse $input3;
      
        my $prime_5="";
        my $prime_3="";
        my $prime_5_comp="";
        my $prime_3_comp="";
         
      if ($something{"PAM_location"} eq "3_prime") {
         $prime_5=translate_IUPAC($input5);
         $prime_3=translate_IUPAC($input3);
         $prime_5_comp=comp(translate_IUPAC($input5_rev));
         $prime_3_comp=comp(translate_IUPAC($input3_rev));
      }else{
         $prime_3=translate_IUPAC($input5);
         $prime_5=translate_IUPAC($input3);
         $prime_3_comp=comp(translate_IUPAC($input5_rev));
         $prime_5_comp=comp(translate_IUPAC($input3_rev));
      }
      
      
      
      while ( $cut <= length($whole_seq) ) {
            push @cuts, $cut;
            $cut = $cut + $cutnumber;
      }
      #################################################################################################################################################################################
      # cut the sequence into equal peaces, so that each forked child can work on one part (paralell!)
      #################################################################################################################################################################################
      
      foreach $cut (@cuts) {
            my $seq = substr( $whole_seq, $cut, $cutnumber );  
            $pm->start and next;            
            my %CRISPR_hash = ();
            
            ###########################################################################################################################################################################
            # Single Sequence
            ###########################################################################################################################################################################
            
            if ($something{"kind"} eq "single") {
                  
                  #####################################################################################################################################################################
                  # Foward Sequence Calculations
                  #####################################################################################################################################################################
                  
                  my @lengths;
                  LENGTHLOOP:foreach my $length ($minlength..$maxlength){
                        my $re_fwd=$prime_5.'.{'.$length.'}'.$prime_3;
                        my $re_rev=$prime_3_comp.'.{'.$length.'}'.$prime_5_comp;
                        POSLOOP:while ($seq =~ m/$re_rev|$re_fwd/g) {
                                    pos $seq = $-[0] + 1 ;
                                    my @temp=($-[0],length($&));                    
                                    my $crisprseq = substr( $seq, $temp[0], $temp[1] );
                                    my @flank_array = find_base_count( $crisprseq );
                                    $tempstatistics{"Number of total possible designs"}++;
                                    if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                          $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                          $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                          $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                          !($crisprseq=~m/TTTTT/) 
                                    ) {
                                          my $name = $seq_obj->display_id;
                                          $name .= "_" . $count . "_" . $cut;
                                          my @new_score=(0,0,0,0,0);
										  my $doench2014_seq="";
											if ($crisprseq=~m/$re_rev/) {
                                                ${ $CRISPR_hash{$name} }{"strand"} = "minus";
												 if($crisprseq=~m/^CC/){
													if (length($crisprseq)==23) {
																$doench2014_seq=reverse_comp(substr( $seq, $temp[0]-3, 30 ));
																my $Xu_seq=reverse_comp(substr( $seq, $temp[0]-7, 30 ));
																if (length($doench2014_seq)==30 && length($Xu_seq)==30) {
																	$new_score[3]+=calc_doench_score($doench2014_seq);
																	$new_score[4]+=calc_XU_score($Xu_seq); 
																}
													}
												 }
                                          }else{
                                                ${ $CRISPR_hash{$name} }{"strand"} = "plus";
                                                if($crisprseq=~m/GG$/){
                                                      if (length($crisprseq)==23) {
															$doench2014_seq=substr( $seq, $temp[0]-4, 30 );
															my $Xu_seq=substr( $seq, $temp[0], 30 );
															if (length($doench2014_seq)==30 && length($Xu_seq)==30) {
																$new_score[3]+=calc_doench_score($doench2014_seq);
																$new_score[4]+=calc_XU_score($Xu_seq);
															}
                                                      }
                                                }    
                                          }                                       
                                          
											if($something{"scores"} eq "xu_score"){
												$new_score[2]=$new_score[4];
											}elsif($something{"scores"} eq "doench_old"){
												$new_score[2]=$new_score[3];
											}else{
												$annonymous_funct=eval($something{"custom_score"});
												$new_score[2]=$annonymous_funct->($doench2014_seq);
											}
                                          ${ $CRISPR_hash{$name} }{"start"} = ($temp[0]) + $cut-$correction;
                                          ${ $CRISPR_hash{$name} }{"end"} = ( $temp[0] + $temp[1] ) + $cut-$correction;
                                          ${ $CRISPR_hash{$name} }{"length"} = $temp[1];
                                          my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset;
                                          my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset;
                                          my %scoring ;
                                          if(${ $CRISPR_hash{$name} }{"strand"} == "plus"){
                                            %scoring = calculate_CRISPR_score(\%trees, \%something, ($end-5), ($end-5), $chrom, 1, \@new_score,$gene_id);
                                          }else{
                                            %scoring = calculate_CRISPR_score(\%trees, \%something, ($start-5), ($start-5), $chrom, 1, \@new_score,$gene_id);
                                          }
                                          
                                          
                                          #############################################################################################################################################
                                          #Statistics
                                          #############################################################################################################################################
                                          
                                          if ($something{"retrieve_recomb_matrix"} ==1 ) {
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                          }
                                          
                                          %{ ${ $CRISPR_hash{$name} }{"context"} } = %scoring;
                                          ${ $CRISPR_hash{$name} }{"nucseq"} = $crisprseq;
                                         
                                          $count++;
                                          
                                          if (make_CRISPR_statistics(\%something, \%scoring, $dont_asses_context, \%tempstatistics) == 1){
                                                delete $CRISPR_hash{$name};
                                          }
                                          #############################################################################################################################################
                                          
                                    } else {
                                          $tempstatistics{"Number of designs excluded because their nucleotide composition was too invariable or contained TTTTT"}++;
                                    }
                              }
                        }                  
            } else{
                  
                  #####################################################################################################################################################################
                  # Double Sequence - only forward calculations needed
                  #####################################################################################################################################################################
                  my %Gpos = make_pos_index( \$seq, "G" );
					my %Cpos = make_pos_index( \$seq, "C" );
					my %Apos = make_pos_index( \$seq, "A" );
					my %Tpos = make_pos_index( \$seq, "T" );
					my %combined;
					% {$combined{"G"}}=%Gpos;
					% {$combined{"A"}}=%Apos;
					% {$combined{"C"}}=%Cpos;
					% {$combined{"T"}}=%Tpos;
				  my %dont_care_ind;
				  my %dont_care_ind_right;
				  my %PAMindex;
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
                  my %PAMindex_right;
                  if ($something{"PAM"} eq "NAG") {
                        %PAMindex=%Tpos;
                        %PAMindex_right=%Apos;
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
                                          $tempstatistics{"Number of total possible designs"}++;
                                          if (  $something{"min_A"} < $flank_array[0] && $something{"max_A"} > $flank_array[0] &&
                                                $something{"min_C"} < $flank_array[1] && $something{"max_C"} > $flank_array[1] &&
                                                $something{"min_T"} < $flank_array[2] && $something{"max_T"} > $flank_array[2] &&
                                                $something{"min_G"} < $flank_array[3] && $something{"max_G"} > $flank_array[3] &&
                                                !($completeseq=~/TTTTT/)
                                          ) {
                                          my $name = ($seq_obj->display_id)."_" . $count . "_" . $cut;
											 my @new_score=(0,0,0,0);
												if (defined $something{"scoring_module"}) {
													$something{"custom_score"}="";
													open my $custom_score_file ,"<", $something{"scoring_module"};
														foreach my $line (<$custom_score_file>){
															$something{"custom_score"}.=$line;
														}
													close $custom_score_file;
													$annonymous_funct=eval($something{"custom_score"});
													$new_score[3]=$annonymous_funct->(reverse_comp(substr( $seq, $Cposind-3, 30)))+$annonymous_funct->(substr( $seq, ( $Cposind + $length + 1 + $spacerlength-4),30));
												}
												$new_score[2]=calc_doench_score(reverse_comp(substr( $seq, $Cposind-3, 30)))+calc_doench_score(substr( $seq, ( $Cposind + $length + 1 + $spacerlength-4),30));
                                                if($something{"scores"} eq "custom"){
													$new_score[2]=$new_score[3];
												}
												@{${ $CRISPR_hash{$name} }{"lengthcombo"}}=($length,$spacerlength);
                                                ${ $CRISPR_hash{$name} }{"start"} = ($Cposind) + $cut;
                                                ${ $CRISPR_hash{$name} }{"end"} = ( $Cposind + $length+$spacerlength+$length+2 + 2 ) + $cut;
                                                ${ $CRISPR_hash{$name} }{"length"} =  $length+$spacerlength+$length+2 + 2;
                                                my $start = ${ $CRISPR_hash{$name} }{"start"} + $location_offset - 500;
                                                my $end = ${ $CRISPR_hash{$name} }{"end"} + $location_offset - 500;
                                                my %scoring;
                                                 %scoring = calculate_CRISPR_score(\%trees, \%something, ($end-5), ($end-5), $chrom, 0, \@new_score, $gene_id);
                                                
                                                #######################################################################################################################################
                                                #Statistics
                                                #######################################################################################################################################
                                                
                                                if (make_CRISPR_statistics(\%something, \%scoring, $dont_asses_context, \%tempstatistics) == 1){
                                                      delete $CRISPR_hash{$name};
                                                      next LENGTHLOOP;
                                                }
                                                
                                                if ($something{"retrieve_recomb_matrix"} ==1 ) {
                                                      ${ ${ $CRISPR_hash{$name} }{"homology"} }{"left"} = substr( $whole_seq, ( ${ $CRISPR_hash{$name} }{"start"} - $something{"left_homology"} ), ($something{"left_homology"}) );
                                                      ${ ${ $CRISPR_hash{$name} }{"homology"} }{"right"} = substr( $whole_seq, ${ $CRISPR_hash{$name} }{"end"}, $something{"right_homology"} );
                                                }
                                                
                                                %{ ${ $CRISPR_hash{$name} }{"context"} } = %scoring;
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
                  write_file( $temp_dir . "/" .$seq_obj->display_id .'_'. $cut . '.json', { binmode => ':raw' }, $json );
                  $json = JSON::XS::encode_json(\%tempstatistics);
                  write_file( $temp_dir . "/" . $seq_obj->display_id .'_'. $cut . 'stats.json', { binmode => ':raw' }, $json );
            }
            
            $pm->finish();
      }
      
      ##########################################################################################################################################################################
      #parent wait till all children are done and then rebuild the CRISPR and the Statistics out of the temporary files
      ##########################################################################################################################################################################
      
      $pm->wait_all_children();
      foreach  my $cut (@cuts) {
            my $json = read_file( $temp_dir . "/" .$seq_obj->display_id .'_'. $cut . '.json', { binmode => ':raw' } );
            %finished_CRISPR_hash = ( %finished_CRISPR_hash, %{ decode_json $json } );
           unlink $temp_dir . "/" . $seq_obj->display_id .'_'. $cut . ".json";
            $json = read_file( $temp_dir . "/" . $seq_obj->display_id .'_'. $cut . 'stats.json', { binmode => ':raw' } );
            my %sechash=%{ decode_json $json };
            foreach  my $seckey (keys(%sechash)){
                        $tempstatistics{$seckey}+=$sechash{$seckey};
            }
           unlink $temp_dir . "/" . $seq_obj->display_id .'_'. $cut . "stats.json";
      }
      return (\%finished_CRISPR_hash,\%tempstatistics);
}




sub make_database{
        if (can_run('rsync') && can_run('wget') && (can_run('bowtie-build') && can_run('bowtie2-build'))) {                
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
		#unlink('temp.log');
	    system('
			rsync -av --progress '.$_[1].'gtf/'.$_[0].'/ .;
			rsync -av --progress --exclude "*abinitio*" '.$_[1].'fasta/'.$_[0].'/cdna/ .;
			rsync -av --progress --exclude "*primary_assembly*" --exclude "*dna_rm*" --exclude "*dna.chromosome*" --exclude "*dna_sm*" '.$_[1].'fasta/'.$_[0].'/dna/ .;
		');
	    print "All files were dowloaded\n";
        if (defined $something{"GUI"}) {$mw->update;};
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
            if (defined $something{"GUI"}) {$mw->update;};
			create_mygff($_[0].'.dna.toplevel.fa',$_[0].'.gtf',$_[0].'.all.dna.fa');
			print "All files were converted to gff\n";
            if (defined $something{"GUI"}) {$mw->update;};
			wrap_sequences($_[0].'.cdna.all.fa');
			cpg_for_all($_[0].'.dna.toplevel.fa'); #store CpG-island information in csv-files
			system('for f in *.fasta ; do rm $f ;done ;');
			print "The entire genome was checked for CPG islands\n";
            if (defined $something{"GUI"}) {$mw->update;};
			system('rm *all*.fasta;');
			system('rm *.flat;');
			system('rm *.gdx');
			system('rm \#*');
			system('rm *README*');
			system('rm *CHECKSUMS*');
			correct_cdna($_[0].".cdna.all.fa"); #change header of [organsim].cdna.all.fa in [organsim].cdna.all.facorrected.fa
			system('for f in *.cdna.all.fa ; do rm $f ;done ;');
			print "CDNA files were corrected\n";
            if (defined $something{"GUI"}) {$mw->update;};
			system("mv ".$_[0].".cdna.all.facorrected.fa ".$_[0].".cdna.all.fa;"); #rename [organsim].cdna.all.facorrected.fa
			system('for f in *.gff ; do rm $f ;done ;');
			 print "Annotation were formatted\n";
             if (defined $something{"GUI"}) {$mw->update;};
			include_cpg("."); #add CpG-island information to mygff-files
			system('for f in *.csv ; do rm $f ;done ;');
			system('for f in *.gtf ; do rm $f ;done ;');
			print "All prerequisites for building alignment indeces were built correctly.\nNow Bowtie indeces will be build, depending on the size of the target genome, this can take a while.\n";
			if (defined $something{"GUI"}) {$mw->update;};
            if(can_run('bowtie2-build')){system("bowtie2-build ".$_[0].".cdna.all.fa ".$_[0].".cdna & bowtie2-build ".$_[0].".all.dna.fa ".$_[0].".dna & bowtie2-build ".$_[0].".dna.toplevel.fa ".$_[0].".genome;");} #create files with bowtie2-indices
			if(can_run('bowtie-build')){system("bowtie-build ".$_[0].".cdna.all.fa ".$_[0].".cdna & bowtie-build ".$_[0].".all.dna.fa ".$_[0].".dna & bowtie-build ".$_[0].".dna.toplevel.fa ".$_[0].".genome;"); }#create files with bowtie-indices
			opendir my $curr_dir , ".";
			while (readdir($curr_dir)) {
				if (-z $_) {
					#unlink($_);
				}				
			}
			closedir($curr_dir);
            my $pwd = cwd();
            if (defined $something{"GUI"}) {$mw->update;};
			print "The database for the organism ".$_[0]." has been built in the following path:\n$pwd";
			if (defined $something{"GUI"}) {
                 $make_database->messageBox(
                        -icon => 'info',
                        -type => 'ok',
                        -title => 'Info',
                        -message => "The database for the organism ".$_[0]." has been built in the following path:\n$pwd",
                    );
            }
            
           
			}else{
                print "rsync, wget, bowtie and bowtie2 need to be installed and executable from the \$PATH variable.\n You can test this by running \"which wget\" and \"which bowtie\" from your terminal."
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
           print $chrom_file_gff "gene_".$id."::".$locus_tag."::".$line[6]."\t".$line[3]."\t".$line[4]."\n";
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
				#unlink($file);
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
		#unlink($filename)
}
#########################################################################################
#name:      comp
#function:  complement a DNA sequence without reversing it
#input:     string
#output:    string
#########################################################################################
sub comp {
      my $comp=$_[0];
      $comp=~ tr/ACGT/TGCA/;                                      
      return $comp;
}
#########################################################################################
#name:      translate_IUPAC
#function:  translate an DNA IUPAC code to ACGT basepaircode in regular expression, perl style
#input:     < string >
#output:    < string >
#########################################################################################
sub   translate_IUPAC {
      my  $reg_exp= $_[0];
      $reg_exp =~ s/U/T/g;
      $reg_exp =~ s/N/[ACGTU]/g;
      $reg_exp =~ s/K/[GTU]/g;
      $reg_exp =~ s/M/[AC]/g;
      $reg_exp =~ s/S/[CG]/g;
      $reg_exp=~ s/W/[ATU]/g;
      $reg_exp =~ s/R/[AG]/g;
      $reg_exp =~ s/Y/[CTU]/g;
      $reg_exp=~ s/B/[^A]/g;
      $reg_exp=~ s/D/[^C]/g;
      $reg_exp=~ s/H/[^G]/g;
      $reg_exp =~ s/V/[ACG]/g;
      $reg_exp=~ s/N/[ACGTU]/g;
      return $reg_exp;
}
#########################################################################################
#name:      from_pam_to_fasta_combis
#function:  find every possble ACGT sequence out of any given IUPAC sequence of any given length
#input:     string
#output:    string
#########################################################################################
sub from_pam_to_fasta_combis{
        my %translator;
        @{$translator{"U"}}="T";
        @{$translator{"A"}}="A";
        @{$translator{"G"}}="G";
        @{$translator{"T"}}="T";
        @{$translator{"C"}}="C";
        @{$translator{"N"}}=("A","C","G","T");
        @{$translator{"K"}}=("G","T");
        @{$translator{"M"}}=("A","C");
        @{$translator{"S"}}=("C","G");
        @{$translator{"W"}}=("A","T");
        @{$translator{"R"}}=("A","G");
        @{$translator{"Y"}}=("C","T");
        @{$translator{"B"}}=("G","C","T");
        @{$translator{"D"}}=("G","A","T");
        @{$translator{"H"}}=("C","A","T");
        @{$translator{"V"}}=("A","C","G");        
        my @old_words=($_[0]);
        my @words=();
        my @word_split=();
        my @tmp=();
        my $pos=0;
        while ($pos<length($_[0])) {
            @words=();
            foreach my $word (@old_words){
                @word_split=split("",$word);
                @tmp=();
                foreach my $translate (@{$translator{$word_split[$pos]}}){
                    @tmp=@word_split;
                    $tmp[$pos]=$translate;
                    push @words , join("",@tmp);
                }
            }
            @old_words=@words;
            $pos++;
        }
        return(@old_words);        
}
#########################################################################################
#name:      rev_com_IUPAC
#function:  reverse complement IUPAC nucleobases to IUPAC nulceobases
#input:     < string >
#output:    < string >
#########################################################################################
sub   rev_com_IUPAC {      
      my $rev = reverse $_[0] ;
      $rev =~ s/U/T/g ;
      $rev =~ tr/ACGTacgtNKMRYBVDH/TGCAtgcaNMKYRVBHD/ ;
      return $rev;
}

#########################################################################################
#name:      calc_doench_score
#function:  Calculate CRISPR Score after Doench et al. 2014 Rational design of highly active sgRNAs for CRISPR-Cas9Ðmediated gene inactivation
#input:     < string > #lengt 30 mandatory
#output:    <numeric double>
#########################################################################################
sub calc_doench_score{
    my $score;
    if (length($_[0])==30) {
    my %sing_nuc_hash = ('G2'=>-0.275377128,'A3'=>-0.323887456,'C3'=>0.172128871,'C4'=>-0.100666209,'C5'=>-0.20180294, 
                    'G5'=>0.245956633,'A6'=>0.036440041,'C6'=>0.098376835,'C7'=>-0.741181291,
                    'G7'=>-0.393264397,'A12'=>-0.466099015,'A15'=>0.085376945,'C15'=>-0.013813972,
                    'A16'=>0.272620512,'C16'=>-0.119022648,'T16'=>-0.285944222,'A17'=>0.097454592,
                    'G17'=>-0.17554617,'C18'=>-0.345795451,'G18'=>-0.678096426,'A19'=>0.22508903,
                    'C19'=>-0.507794051,'G20'=>-0.417373597,'T20'=>-0.054306959,'G21'=>0.379899366,
                    'T21'=>-0.090712644,'C22'=>0.057823319,'T22'=>-0.530567296,'T23'=>-0.877007428,
                    'C24'=>-0.876235846,'G24'=>0.278916259,'T24'=>-0.403102218,'A25'=>-0.077300704,
                    'C25'=>0.287935617,'T25'=>-0.221637217,'G28'=>-0.689016682,'T28'=>0.117877577,
                    'C29'=>-0.160445304,'G30'=>0.386342585);
    my %dinuc_hash = ('GT2'=>-0.625778696,'GC5'=>0.300043317,'AA6'=>-0.834836245,'TA6'=>0.760627772,'GG7'=>-0.490816749,
                      'GG12'=>-1.516907439,'TA12'=>0.7092612,'TC12'=>0.496298609,'TT12'=>-0.586873894,'GG13'=>-0.334563735,
                      'GA14'=>0.76384993,'GC14'=>-0.53702517,'TG17'=>-0.798146133,'GG19'=>-0.66680873,'TC19'=>0.353183252,
                      'CC20'=>0.748072092,'TG20'=>-0.367266772,'AC21'=>0.568209132,'CG21'=>0.329072074,'GA21'=>-0.836456755,
                      'GG21'=>-0.782207584,'TC22'=>-1.029692957,'CG23'=>0.856197823,'CT23'=>-0.463207679,'AA24'=>-0.579492389,
                      'AG24'=>0.649075537,'AG25'=>-0.077300704,'CG25'=>0.287935617,'TG25'=>-0.221637217,'GT27'=>0.117877577,
                      'GG29'=>-0.697740024);
    my $gc = ( substr($_[0],4,20) =~ tr/GC/GC/);
    if ($gc < 10){
        $score=0.597636154+(abs($gc-10)*-0.202625894)
    }else{
        $score=0.597636154+(($gc-10)*-0.166587752)
    }        
    foreach my $i (0..29){        
       my $key = substr($_[0],$i,1).($i+1);
       if ($sing_nuc_hash{$key}) {
        $score+=$sing_nuc_hash{$key};
       }
       if($i<29){
        $key =substr($_[0],$i,2).($i+1);
        if ($dinuc_hash{$key}){
                $score+=$dinuc_hash{$key};
        }
       }
    }
    return(1/(1+exp(-$score)))
      #code
    }else{
        return(0);
    }
}
#########################################################################################
#name:      calc_XU_score
#function:  Calculate CRISPR Score after XU et al. 2015 Sequence determinants of improved CRISPR sgRNA design
#input:     < string > #lengt 30 mandatory 20 Protspacer followed by 10 including NGG PAM
#output:    <numeric double>
#########################################################################################
sub calc_XU_score{
    my $score=0;  
    if (length($_[0])==30) {
        my %scoring_matrix;
        @{$scoring_matrix{'A'}}=(0,0,0,0,0.025840846,0,0,0,0.02156311,0.129118902,0.030483786,0.093646913,0,0,0.202820147,0.129158071,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
        @{$scoring_matrix{'C'}}=(0,0,-0.113781378,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.23502822,0,-0.125927965,0,0,0,0.179639101,0,0,0,0,0,0);
        @{$scoring_matrix{'G'}}=(0,0,0,0.080289971,0.072680697,0.100642827,0.082839514,0,0,0,0,0,0,-0.214271553,0,0,0.107523301,0,0.238517854,0.353047311,0,0,0,0,0,0,0,0,0,0);
        @{$scoring_matrix{'T'}}=(0,0,0,0,0,0,-0.070933894,0,0,0,-0.169986128,0,0,0.073750154,0,0,-0.349240474,-0.145493093,-0.300975354,-0.221752041,-0.155910373,0,0,0,0,0,-0.116646129,0,0,0);
        my $pos=0;
        while ( $_[0]=~m/(\w)/g) {
            $score+=@{$scoring_matrix{$1}}[$pos];
            $pos++;
        }
        return(($score-(-0.5))/(2.5));
    }else{
        return(0);
    }
}
#########################################################################################
#name:      rndStr
#function:  generate a random string of n characters from a array A
#input:     n < int >, A <char array>
#output:    <string>
#########################################################################################
sub rndStr{ join"", @_[ map{ rand @_ } 1 .. shift ] };

sub mock{};

1; 

