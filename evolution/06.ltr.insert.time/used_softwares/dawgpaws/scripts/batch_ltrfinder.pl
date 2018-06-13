#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_ltrfinder.pl - Run ltr_finder in batch mode         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 10/03/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the LTR_FINDER program in batch mode and extract     |
#  output to gff files, summary text files, fasta files and |
#  a tab delimited summary sequence file.                   |
#                                                           |
# VERSION: Release 1.0                                      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TO DO: -Reset vals to null and only print output when expected vals
#         are not null
#        -No gff output if no hits in ltr_finder
#        -Add possiblity to extract sequence data

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Input dir of fasta files
my $outdir;                    # Base output dir 
my $config_file;               # Configuration file

# GLOBAL SEQUENCE STRING OBJECTS
my $seq_in;                    # Full fasta file input object
my $qry_seq;                   # The full query sequence object

# LTR Finder Parameters
my @lf_params = ();            # 2d Array to hold the find_ltr parameters 

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff_convert = 1;
my $do_feat_seq = 0;           # Extract feature sequence data

# Counters/Index Vals
my $i = 0;                     # Array index val
my $file_num = 0;              # File number
my $proc_num = 0;              # Process number

# Vars needing 
my $name_root;
my $gff_dir;
my $ls_out;
my $fl_gff_outpath;

#-----------------------------+
# USER ENV OPTIONS            |
#-----------------------------+
# If variable not defined in the ENV then
# assume it is in the user's path and that the Oryza tRNA database
# is the preferred database
my $trna_db = $ENV{TRNA_DB} || 
    "Os-tRNAs.fa";
my $prosite_dir = $ENV{PROSITE_DIR} ||
    "ps_scan";
my $lf_path = $ENV{LTR_FINDER} ||
    "ltr_finder";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    "c|config=s"   => \$config_file,
		    # ADDITIONAL OPTIONS
		    "ltr-finder=s" => \$lf_path,     # LTR Finder Path
		    "s|trna-db=s"  => \$trna_db,     # s as per program
		    "a|prosite=s"  => \$prosite_dir, # a as per program
		    "f|feat-seq"   => \$do_feat_seq, 
		    "g|gff"        => \$do_gff_convert,
		    "q|quiet"      => \$quiet,
		    "verbose"      => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( ($show_usage) ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_version) {
    print "\nbatch_ltrfinder.pl:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) || (!$config_file) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory must be specified". 
	" at the command line.\n" if (!$indir);
    print STDERR "ERROR: An output directory must be specified". 
	" at the command line.\n" if (!$outdir);
    print STDERR "ERROR: A config file must be specified". 
	" at the command line.\n" if (!$config_file);
    print_help ("usage", $0 );
    exit;
}

#-----------------------------------------------------------+
# IF SEQ FEATURE STRINGS REQUESTED LOAD REQUIRED MODULES    |
#-----------------------------------------------------------+
# This will only attempt to use these objects when they are needed for parsing
# sequence strings.
if ($do_feat_seq) {
    use Bio::Location::Simple;    # Used to get revcom of substrings
    use Bio::Seq;                 # Fetch seq and do substrings
    use Bio::SeqIO;               # Seq IO
    use Bio::Location::Simple;    # Used to get reverse complement when needed
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         | 
#-----------------------------------------------------------+

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# LOAD THE CONFIG FILE        |
#-----------------------------+
$i=0;
my $config_line_num=0;

open (CONFIG, "<$config_file") ||
    die "ERROR Can not open the config file:\n $config_file";

while (<CONFIG>) {
    chomp;
    $config_line_num++;
    unless (m/^\#/) {
       	my @in_line = split (/\t/);           # Implicit split of $_ by tab
	my $num_in_line = @in_line; 
	
	# Can have just a name to run LTR_Finder with default settings
	# or can have two columns with additional parameter options
	# parameter options in the second columns.
	# I will currently stick with the two column config file
	# since there are so many options availabe with LTR_FINDER
	# that a multiple column config file would get messy.
	if ($num_in_line < 3) { 
	    $lf_params[$i][0] = $in_line[0] || "NULL";  # Name
	    $lf_params[$i][1] = $in_line[1] || "";      # Suffix
	    $i++;
	} # End of if $num_in_line is 10
	else {
	    print "\a";
	    print STDERR "WARNING: Unexpected number of line in config".
		" file line $config_line_num\n$config_file\n";
	}

   } # End of unless comment line
} # End of while CONFIG file
close CONFIG;

# Number of parameter sets specified in the config file
my $num_par_sets = $i;

if ($num_par_sets == 0) {
    print "\a";
    print STDERR "ERROR: No parameter sets were found in the config file:\n".
	"$config_file\n";
    exit;
}

#-----------------------------+
# GET THE FASTA FILES FROM THE|
# DIRECTORY PROVIDED BY THE   |
#$indir VARIABLE              |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $num_files = @fasta_files;

if ($num_files == 0) {
    print "\a";
    print "ERROR: No input fasta files were found in the directroy:\n$indir\n";
    print "Fasta files must end with the fasta or fa extension\n";
    exit;
}

my $num_proc_total = $num_files * $num_par_sets;

print STDERR "$num_proc_total find_ltr runs to process\n\n" if $verbose;

for my $ind_file (@fasta_files) {

    $file_num++;

    # Get root file name
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }

    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE LTR_FINDER OUTDIR    |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $ltrfinder_dir = $name_root_dir."ltr_finder/";
    unless (-e $ltrfinder_dir) {
	mkdir $ltrfinder_dir ||
	    die "Could not create ltr_finder out dir:\n$ltrfinder_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    if ($do_gff_convert) {
	$gff_dir = $name_root_dir."gff/";
	unless (-e $gff_dir) {
	    mkdir $gff_dir ||
		die "Could not create gff out dir:\n$gff_dir\n";
	}
    }

    #-----------------------------------------------------------+
    # LOAD SEQUENCE TO STRING IF WE WANT TO EXTRACT FEATURES    |
    #-----------------------------------------------------------+
    # TODO: This should also do the test to make sure that this is a fasta file
    #       with a single record
    if ($do_feat_seq) {

	my $lf_qryseqpath = $indir.$ind_file;

	$seq_in = Bio::SeqIO->new ( '-format' => 'fasta',
				    '-file' => "<$lf_qryseqpath") 
	    || die "Can not open sequence file $lf_qryseqpath";
	
	# THE FOLLOWING ERROR CHECK WOULD BE HELPFUL
	# If number of seqs in object is one move forward
	# otherwise throw error.
	# Try using next_seq outside of a while statement
	# This appears to work.
	$qry_seq = $seq_in->next_seq();

    }

    #-----------------------------------------------------------+
    # RUN LTR_FINDER FOR EACH SET OF PARAM VALS IN CONFIG FILE |
    #-----------------------------------------------------------+
    for ($i=0; $i<$num_par_sets; $i++) {
	
	$proc_num++;
	
	# Load parameters from array refs to name vars
	my $lf_gff_suffix = $lf_params[$i][0]; # Name of the parameter set
	my $lf_cmd_suffix = $lf_params[$i][1]; # Name of the parameter set

	# Status proces number
	my $stat_par_num = $i+1;

	print STDERR "\n#------------------------------------\n" if $verbose;
	print STDERR "# Processing $name_root $lf_gff_suffix\n" if $verbose;
	print STDERR "# Seq $file_num of $num_files\n" if $verbose;
	print STDERR "# Par $stat_par_num of $num_par_sets\n" if $verbose;
	print STDERR "# PROCESS $proc_num of $num_proc_total\n" if $verbose;
	print STDERR "#------------------------------------\n\n" if $verbose;

	# lf_out is path of the ltrfinder output file
	my $lf_inseqpath = $indir.$ind_file;
	my $lf_out = $ltrfinder_dir.$name_root."_".$lf_gff_suffix.
	    "_ltr_finder.txt";
	my $lf_gff_outpath = $gff_dir.$name_root."_".$lf_gff_suffix.
	    "_ltr_finder.gff";
	my $lf_summary_outpath = $gff_dir.$name_root."_".$lf_gff_suffix.
	    "_lf_summary";
	my $lf_cmd = "$lf_path $lf_inseqpath".
	    " -s $trna_db".
	    " -a $prosite_dir". 
	    " $lf_cmd_suffix".
	    " > $lf_out";

	print STDERR "\tCMD:$lf_cmd\n" if $verbose;

	system ($lf_cmd);

	#-----------------------------+
	# CONVERT OUTPUT TO GFF       |
	#-----------------------------+
	if ( (-e $lf_out) && ($do_gff_convert) ) {
	    ltrfinder2gff ( $name_root, $lf_out, $lf_gff_outpath, 
			    0, $lf_gff_suffix, $do_feat_seq,$ltrfinder_dir,
			    $lf_summary_outpath);
	}


	# Check the global aspect of the seq record
	#print STDERR "SEQ LEN: ".$qry_seq->length()."\n" if $verbose;

    }


} # End of for each individual fasta file



exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub print_help {
    my ($help_msg, $podfile) =  @_;
    # help_msg is the type of help msg to use (ie. help vs. usage)
    
    print "\n";

    #-----------------------------+
    # PIPE WITHIN PERL            |
    #-----------------------------+
    # This code made possible by:
    # http://www.perlmonks.org/index.pl?node_id=76409
    # Tie info developed on:
    # http://www.perlmonks.org/index.pl?node=perltie 
    #
    #my $podfile = $0;
    my $scalar = '';
    tie *STDOUT, 'IO::Scalar', \$scalar;
    
    if ($help_msg =~ "usage") {
	podselect({-sections => ["SYNOPSIS|MORE"]}, $0);
    }
    else {
	podselect({-sections => ["SYNOPSIS|ARGUMENTS|OPTIONS|MORE"]}, $0);
    }

    untie *STDOUT;
    # now $scalar contains the pod from $podfile you can see this below
    #print $scalar;

    my $pipe = IO::Pipe->new()
	or die "failed to create pipe: $!";
    
    my ($pid,$fd);

    if ( $pid = fork() ) { #parent
	open(TMPSTDIN, "<&STDIN")
	    or die "failed to dup stdin to tmp: $!";
	$pipe->reader();
	$fd = $pipe->fileno;
	open(STDIN, "<&=$fd")
	    or die "failed to dup \$fd to STDIN: $!";
	my $pod_txt = Pod::Text->new (sentence => 0, width => 78);
	$pod_txt->parse_from_filehandle;
	# END AT WORK HERE
	open(STDIN, "<&TMPSTDIN")
	    or die "failed to restore dup'ed stdin: $!";
    }
    else { #child
	$pipe->writer();
	$pipe->print($scalar);
	$pipe->close();	
	exit 0;
    }
    
    $pipe->close();
    close TMPSTDIN;

    print "\n";

    exit 0;
   
}

sub ltrfinder2gff {
    
    # seq_id could be extracted from the ltrfinder result
    # but passing it to the subfunction directly allows for cases
    # where the id assigned by ltrfinder differs from the way
    # the user is referring to the assembly

    # The gff outfile path is passed to the function
    my ($seq_id, $lf_infile, $gffout, $do_append, $gff_suffix,
	$do_feat_seq, $ltrfinder_dir,$sumout_root) = @_;
    
    # FILE PATHS
    my $sumout_file = $sumout_root."_short.txt";
    my $sumout_seq_file = $sumout_root."_seq.txt";
    
    # The gff src id
    my $gff_src = "ltr_finder:".$gff_suffix;
    my $gff_str_out;        # A single out string line of gff out

    # LTF FINDER COUTNERS/INDICES
    my $lf_id_num = 0;       # Incremented for every predicted model

    # LTR_FINDER BOOLEANS
    my $in_emp_details = 0;  # Exact match pairs details
    my $in_ltr_align = 0;    # Details of the LTR alignment
    my $in_pbs_align = 0;
    my $in_ppt;
    
    # By default assume ltr_retro found
    my $has_ltr_retro = 1;   # Does the sequence have an LTR retrotransposon

    
    #
    my $lf_prog_name;             # LTR Finder program name
    my $lf_seq_id;                # Seq id
    my $lf_seq_len;               # Sequence length
    my $lf_version;               # Version of LTR finder being parsed
    my $lf_trna_db;               # tRNA database used
    
    # Genearl LTR Retrotransposon Span variables
    my $lf_strand;                # Strand of the LTR Retro
    my $lf_span_start;            # Location start
    my $lf_span_end;              # Location end
    my $lf_length;                # Length
    my $lf_score;                 # Score
    my $lf_ltr_similarity;        # Similarity of the LTRs
    
    # Status strings
    my $has_5ltr_tg;                 # TG in 5' END of 5' LTR
    my $has_5ltr_ca;                 # CA in 3' END of 5' LTR
    my $has_3ltr_tg;                 # TG in 5' END of 3' LTR
    my $has_3ltr_ca;                 # CA in 3' END of 3' LTR
    my $has_tsr;                     # Has Target Site Replication
    my $has_pbs;                     # Has Primer Binding Site
    my $has_ppt;                     # Has Poly Purine Tract
    my $has_rt;                      # Has Reverse Transcriptase
    my $has_in_core;                 # Has Integrase Core
    my $has_in_cterm;                # Has Integrase C-term
    my $has_rh;                      # Has RNAseH
    
    my $lf_ltr_id;                   # LTR Retro Id number 

    # LTR COORDINATES
    my $lf_5ltr_start;
    my $lf_5ltr_end;
    my $lf_5ltr_len;
    my $lf_3ltr_start;
    my $lf_3ltr_end;
    my $lf_3ltr_len;

    # TSR COORDINATES
    my $lf_5tsr_start;             # Start of the 5' TSR
    my $lf_5tsr_end;               # End of the 5' TSR
    my $lf_3tsr_start;             # Start of the 3' TSR
    my $lf_3tsr_end;               # End of the 3' TSR 
    my $lf_tsr_string;             # String of bases in the TSR

    # BOUNDARY SHARPNESS
    my $lf_sharp_5;                # Sharpness of 5' Boundary
    my $lf_sharp_3;                # Sharpness of 3' Boundary
    
    # PBS
    my $lf_pbs_num_match;          # Number of matched bases in the PBS
    my $lf_pbs_aln_len;            # PBS alignment length
    my $lf_pbs_start;              # Start of the PBS signal
    my $lf_pbs_end;                # End of the PBS signal
    my $lf_pbs_trna;               # PBS tRNA type and anti-codon

    # PPT 
    my $lf_ppt_num_match;          # Number of matched based in the PPT
    my $lf_ppt_aln_len;            # PPT alignment length
    my $lf_ppt_start;              # Start of the PPT
    my $lf_ppt_end;                # End of the PPT
    
    #-----------------------------+
    # DOMAIN DATA                 |
    #-----------------------------+
    
    # GENERAL DOMAIN VARS
    my $lf_domain_dom_start;
    my $lf_domain_dom_end;
    my $lf_domain_orf_start;
    my $lf_domain_orf_end;
    my $lf_domain_name;            # Type of the domain
    
    # INTEGRASE CORE
    #my $has_in_core = 0;           # Boolean for has integrase core
    my $lf_in_core_dom_start;
    my $lf_in_core_dom_end;
    my $lf_in_core_orf_start;
    my $lf_in_core_orf_end;
    
    # INTEGRASE C-TERM
    #my $has_in_cterm = 0;
    my $lf_in_cterm_dom_start;
    my $lf_in_cterm_dom_end;
    my $lf_in_cterm_orf_start;
    my $lf_in_cterm_orf_end;

    # RNASEH
    #my $has_rh = 0;
    my $lf_rh_dom_start;
    my $lf_rh_dom_end;
    my $lf_rh_orf_start;
    my $lf_rh_orf_end;
    
    # RT
    #my $has_rt = 0;
    my $lf_rt_dom_start;
    my $lf_rt_dom_end;
    my $lf_rt_orf_start;
    my $lf_rt_orf_end;
    
    #-----------------------------+
    # OPEN SUMMARY OUTFILE        |
    #-----------------------------+
    open (SUMOUT, ">$sumout_file" ) ||
	die "ERROR: Can not open summary outfile:n\ $sumout_file\n";

    open (SUMSEQ, ">$sumout_seq_file" ) ||
	die "ERROR Can not open summary seq outfile:\n $sumout_seq_file\n";

    #-----------------------------+
    # OPEN GFF OUTFILE            |
    #-----------------------------+
    if ($do_append) {
	open (GFFOUT, ">>$gffout") ||
	    die "ERROR: Can not open gff outfile:\n $gffout\n";
    }
    else {
	open (GFFOUT,">$gffout") ||
	    die "ERROR: Can not open gff outfile:\n $gffout\n";
    }
    
    # Test of squence in the subfunction
    if ($do_feat_seq) {
	print "SEQ NAME: ".$qry_seq->primary_id()."\n" if $verbose;
    }
    
    #-----------------------------+
    # OPEN INPUT FILE             |
    #-----------------------------+
    open (INFILE, "<$lf_infile") ||
	die "ERROR: Can not open LTR_FINDER result file\n $lf_infile\n";
    
    while (<INFILE>) {
	chomp;
	# The following can be used to print LTR Finder output as
	# it is parsed
        #    print STDERR $_."\n";
	
	
	# CHECK BOOLEANS
	
	
	# 
	if (m/Nothing Header(.*)/) {
	    
	}

	elsif (m/No LTR Retrotransposons Found/) {
	    # No LTR Retrotransposons Found
	    # Set has_ltr_retro to false
	    $has_ltr_retro = 0;
	    print STDERR "\n\n\n\nHAS NO LTR RETROS\n\n\n\n" if $verbose;
	}
	
	#///////////////////////////////////////
	# PAY ATTENTION BELOW THIS MAY BE THE
	# PLACE TO PRINT OUT SAVED GFF DATA
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
	
	# IN NEW REC, GET ID
	elsif (m/^\[(.*)\]/) {
	    
	    
	    
	    #-----------------------------------------------------------+
	    # PRINT STORED GFF OUTPUT AND SUMMARY OUT                   |
	    #-----------------------------------------------------------+
	    unless ($1 == 1 ) {

		#-----------------------------+
		# SUMMARY OUTPUT              |
		#-----------------------------+
		print SUMOUT $seq_id."_LTRFinder_".$gff_suffix.$lf_ltr_id."\t".
		    $lf_span_start."\t".     # 2
		    $lf_span_end."\t".       # 3
		    $lf_length."\t".         # 4
		    $lf_score."\t".          # 5
		    # LTR INFO
		    $lf_ltr_similarity."\t". # 6
		    $lf_5ltr_len."\t".       # 7
		    $lf_3ltr_len."\t".       # 8
		    $lf_sharp_5."\t".        # 9
		    $lf_sharp_3."\t".        # 10
		    # BOOLEANS
		    $has_5ltr_tg."\t".       # 11
		    $has_5ltr_ca."\t".       # 12
		    $has_3ltr_tg."\t".       # 13
		    $has_3ltr_ca."\t".       # 14
		    $has_tsr."\t".           # 15
		    $has_pbs."\t".           # 16
		    $has_ppt."\t".           # 17
		    $has_rt."\t".            # 18
		    $has_in_core."\t".       # 19
		    $has_in_cterm."\t".      # 20
		    $has_rh.                 # 21
		    "\n";                    # NEWLINE
		    
		
		print SUMSEQ $seq_id."_LTRFinder_".$gff_suffix.$lf_ltr_id."\t";

		# Two alternatives here
		# keep appending to gff_lout or LTR Retro object 
		# variable and then print set
		# or print line at a time .. currently printing 
		# a line at a time. JCE 10/02/2007
		# Advantage to this approach is that it is not memory
		# limited. Extracted LTR sequences for entire chromosomes
		# could run into mememory limites using an object
		# oriented approach.

		#-----------------------------+
		# FULL SPAN OF LTR RETRO      |
		#-----------------------------+
		$gff_str_out = "$lf_seq_id\t".   # Seq ID
		    "$gff_src\t".                # Source
		    "LTR_retrotransposon\t".     # Data type
		    "$lf_span_start\t".          # Start
		    "$lf_span_end\t".            # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
		# Extract sequence
		if ($do_feat_seq) {
		    my $feat_strand;
		    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			$lf_ltr_id."_LTR_retrotransposon"; 
		    my $out_seq_path = $ltrfinder_dir.$seq_header.".fasta";
		    
		    print ">$seq_header\n";
		    
		    if ($lf_strand =~ "-") {
			$feat_strand="-1";
		    }
		    else {
			$feat_strand="1";
		    }
		    
		    # TEMP PRINT OF FEATURE STRAND
		    print STDERR "STRAND: $feat_strand\n" if $verbose;
		    
		    my $loc = Bio::Location::Simple->
			new(-start  => $lf_span_start,
			    -end    => $lf_span_end,
			    -strand => $feat_strand);
		    # Sequene of the feature
		    my $feat_seq = $qry_seq->subseq($loc);
		    
		    open (SEQOUT,">$out_seq_path") ||
			die "Can not open\n:$out_seq_path\n";
		    print SEQOUT ">$seq_header\n";
		    print SEQOUT "$feat_seq\n";
		    close SEQOUT;

		    print SUMSEQ "$feat_seq\t";
		}
		# END EXTRACT SEQUENCE
		
		#-----------------------------+
		# PRIMER BINDING SITE         |
		#-----------------------------+
		if ($has_pbs) {
		    $gff_str_out = "$lf_seq_id\t".  # Seq ID
			"$gff_src\t".               # Source
			"primer_binding_site\t".    # Data type
			"$lf_pbs_start\t" .         # Start
			"$lf_pbs_end\t".            # End
			"$lf_score\t".              # Score
			"$lf_strand\t".             # Strand
			".\t".                      # Frame
			"ltr_finder_$lf_ltr_id\n";  # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    # Extract sequence
		    if ($do_feat_seq) {
			my $feat_strand;
			my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			    $lf_ltr_id."_primer_binding_site"; 
			my $out_seq_path = $ltrfinder_dir.
			    $seq_header.".fasta";
			
			print ">$seq_header\n";
			
			if ($lf_strand =~ "-") {
			    $feat_strand="-1";
			}
			else {
			    $feat_strand="1";
			}
			
			my $loc = Bio::Location::Simple->
			    new(-start  => $lf_pbs_start,
				-end    => $lf_pbs_end,
				-strand => $feat_strand);
			
			# Sequene of the feature
			my $feat_seq = $qry_seq->subseq($loc);
			
			open (SEQOUT,">$out_seq_path") ||
			    die "Can not open\n:$out_seq_path\n";
			print SEQOUT ">$seq_header\n";
			print SEQOUT "$feat_seq\n";
			close SEQOUT;
			
			print SUMSEQ "$feat_seq\t"; 
			
		    }
		    # END EXTRACT SEQUENCE
		} #
		else {
		    # Case where we want to extract seq features but
		    # PBS is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		
		#-----------------------------+
		# POLYPURINE TRACT            |
		#-----------------------------+
		if ($has_ppt) {
		    $gff_str_out = "$lf_seq_id\t".  # Seq ID
			"$gff_src\t".               # Source
			"RR_tract\t".               # Data type
			"$lf_ppt_start\t".          # Start
			"$lf_ppt_end\t".            # End
			"$lf_score\t".              # Score
			"$lf_strand\t".             # Strand
			".\t".                      # Frame
			"ltr_finder_$lf_ltr_id\n";  # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    # Extract sequence
		    if ($do_feat_seq) {
			my $feat_strand;
			my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			    $lf_ltr_id."_RR_tract"; 
			my $out_seq_path = $ltrfinder_dir.
			    $seq_header.".fasta";
			
			print ">$seq_header\n";
			
			if ($lf_strand =~ "-") {
			    $feat_strand="-1";
			}
			else {
			    $feat_strand="1";
			}
			
			my $loc = Bio::Location::Simple->
			    new(-start  => $lf_ppt_start,
				-end    => $lf_ppt_end,
				-strand => $feat_strand);
			# Sequene of the feature
			my $feat_seq = $qry_seq->subseq($loc);
			
			open (SEQOUT,">$out_seq_path") ||
			    die "Can not open :\n$out_seq_path\n";
			print SEQOUT ">$seq_header\n";
			print SEQOUT "$feat_seq\n";
			close SEQOUT;
		    }
		    # END EXTRACT SEQUENCE
		}
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		 
		
		#-----------------------------+
		# 5' LTR                      |
		#-----------------------------+
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "five_prime_LTR\t".         # Data type
		    "$lf_5ltr_start\t".         # Start
		    "$lf_5ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
		
		# Extract sequence
		if ($do_feat_seq) {
		    my $feat_strand;
		    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			$lf_ltr_id."_five_prime_LTR"; 
		    my $out_seq_path = $ltrfinder_dir.$seq_header.".fasta";
		    
		    print ">$seq_header\n";
		    
		    if ($lf_strand =~ "-") {
			$feat_strand="-1";
		    }
		    else {
			$feat_strand="1";
		    }
		    
		    my $loc = Bio::Location::Simple->new
			(-start  => $lf_5ltr_start,
			 -end    => $lf_5ltr_end,
			 -strand => $feat_strand);
		    # Sequene of the feature
		    my $feat_seq = $qry_seq->subseq($loc);
		    
		    open (SEQOUT,">$out_seq_path") ||
			die "Can not open :\n$out_seq_path\n";
		    print SEQOUT ">$seq_header\n";
		    print SEQOUT "$feat_seq\n";
		    close SEQOUT;

		    # PRINT TO SUMMARY SEQUENCE FILE
		    print SUMSEQ "$feat_seq\t";
		}
		# END EXTRACT SEQUENCE
		
		#-----------------------------+
		# 3' LTR                      |
		#-----------------------------+
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "three_prime_LTR\t".        # Data type
		    "$lf_3ltr_start\t".         # Start
		    "$lf_3ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
		
		# Extract sequence
		if ($do_feat_seq) {
		     my $feat_strand;
		     my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			 $lf_ltr_id."_three_prime_LTR"; 
		     my $out_seq_path = $ltrfinder_dir.$seq_header.".fasta";
		     
		     print ">$seq_header\n";
		     
		     if ($lf_strand =~ "-") {
			 $feat_strand="-1";
		     }
		     else {
			 $feat_strand="1";
		     }
		     
		     my $loc = Bio::Location::Simple->
			 new(-start  => $lf_3ltr_start,
			     -end    => $lf_3ltr_end,
			     -strand => $feat_strand);
		     # Sequene of the feature
		     my $feat_seq = $qry_seq->subseq($loc);
		     
		     open (SEQOUT,">$out_seq_path") ||
			 die "Can not open:\n$out_seq_path\n";
		     print SEQOUT ">$seq_header\n";
		     print SEQOUT "$feat_seq\n";
		     close SEQOUT;

		     
		     # PRINT TO SUMMARY SEQUENCE FILE
		     print SUMSEQ "$feat_seq\t";
		 }
		# END EXTRACT SEQUENCE
		
		#-----------------------------+
		# TARGET SITE DUPLICATION     |
		#-----------------------------+
		if ($has_tsr) {
		    
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"target_site_duplication\t". # Data type
			"$lf_5tsr_start\t".          # Start
			"$lf_5tsr_end\t".            # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;	    
		    
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"target_site_duplication\t". # Data type
			"$lf_3tsr_start\t".          # Start
			"$lf_3tsr_end\t".            # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    # Extract sequence
		    if ($do_feat_seq) {
			my $feat_strand;
			my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			    $lf_ltr_id."_target_site_duplication"; 
			my $out_seq_path = $ltrfinder_dir.
			    $seq_header.".fasta";
			
			print ">$seq_header\n";
			
			if ($lf_strand =~ "-") {
			    $feat_strand="-1";
			}
			else {
			    $feat_strand="1";
			}
			
			my $loc = Bio::Location::Simple->
			    new(-start  => $lf_5tsr_start,
				-end    => $lf_5tsr_end,
				-strand => $feat_strand);
			 # Sequene of the feature
			 my $feat_seq = $qry_seq->subseq($loc);
			 
			 open (SEQOUT,">$out_seq_path") ||
			     die "Can not open\n:$out_seq_path\n";
			 print SEQOUT ">$seq_header\n";
			 print SEQOUT "$feat_seq\n";
			 close SEQOUT;
		     }
		     # END EXTRACT SEQUENCE
		 }
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		    
		    

		#-----------------------------+
		# INTEGRASE                   |
		#-----------------------------+
		if ($has_in_core) {
		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"integrase_core_domain\t".   # Data type
			"$lf_in_core_dom_start\t".   # Start
			"$lf_in_core_dom_end\t".     # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    
		    # ONLY DO THE FOLLOWING IF AND INTEGRASE
		    # CORE ORF WAS FOUND. It is possible that a
		    # core was found but a full orf can not be extracted
		    if ($lf_in_core_orf_start) {
			#/////////
			# NOT SONG
			#\\\\\\\\\
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "integrase_core_orf\t".      # Data type
			    "$lf_in_core_orf_start\t".   # Start
			    "$lf_in_core_orf_end\t".     # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
			
			# Extract sequence
			if ($do_feat_seq) {
			    my $feat_strand;
			    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
				$lf_ltr_id."_in_core_orf"; 
			    my $out_seq_path = $ltrfinder_dir.
				$seq_header.".fasta";
			    
			    print ">$seq_header\n";
			    
			    if ($lf_strand =~ "-") {
				$feat_strand="-1";
			    }
			    else {
				$feat_strand="1";
			    }
			    
			    my $loc = Bio::Location::Simple->new
				(-start  => $lf_in_core_orf_start,
				 -end    => $lf_in_core_orf_end,
				 -strand => $feat_strand);
			    # Sequene of the feature
			    my $feat_seq = $qry_seq->subseq($loc);
			    
			    open (SEQOUT,">$out_seq_path") ||
				die "Can not open:\n$out_seq_path\n";
			    print SEQOUT ">$seq_header\n";
			    print SEQOUT "$feat_seq\n";
			    close SEQOUT;
			}
			# END EXTRACT SEQUENCE
		    }
		    else {
			# Case where we want to extract seq features but
			# PPT is not present
			if ($do_feat_seq) {
			    print SUMSEQ ".\t"; 
			}
		    }
		   

		    
		} # End of has in_core
		
		
		#-----------------------------+
		# GAG C-TERMINAL END          |
		#-----------------------------+
		if ($has_in_cterm) {
		    
		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"integrase_cterm_domain\t".  # Data type
			"$lf_in_cterm_dom_start\t".  # Start
			"$lf_in_cterm_dom_end\t".    # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    if ($lf_in_cterm_orf_start) {
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "integrase_cterm_orf\t".     # Data type
			    "$lf_in_cterm_orf_start\t".  # Start
			    "$lf_in_cterm_orf_end\t".    # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
		    }

		} # End of has_in_cterm
		

		#-----------------------------+
		# RNASEH                      |
		#-----------------------------+
		if ($has_rh) {
		    
		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"rnaseh_domain\t".           # Data type
			"$lf_rh_dom_start\t".        # Start
			"$lf_rh_dom_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    if ($lf_rh_orf_start) {
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "rnaseh_orf\t".              # Data type
			    "$lf_rh_orf_start\t".        # Start
			    "$lf_rh_orf_end\t".          # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
			
			# Extract sequence
			if ($do_feat_seq) {
			    my $feat_strand;
			    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
				$lf_ltr_id."_rnaseh_orf"; 
			    my $out_seq_path = $ltrfinder_dir.
				$seq_header.".fasta";
			    
			    print ">$seq_header\n";
			    
			    if ($lf_strand =~ "-") {
				$feat_strand="-1";
			    }
			    else {
				$feat_strand="1";
			    }
			    
			    my $loc = Bio::Location::Simple->
				new(-start  => $lf_rh_orf_start,
				    -end    => $lf_rh_orf_end,
				    -strand => $feat_strand);
			    # Sequene of the feature
			    my $feat_seq = $qry_seq->subseq($loc);
			    
			    open (SEQOUT,">$out_seq_path") ||
				die "Can not open $out_seq_path\n";
			    print SEQOUT ">$seq_header\n";
			    print SEQOUT "$feat_seq\n";
			    close SEQOUT;
			}
			# END EXTRACT SEQUENCE
		    }
		}
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		   
		
		
		#-----------------------------+
		# REVERSE TRANSCRIPTASE       |
		#-----------------------------+
		if ($has_rt) {

		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"rt_domain\t".               # Data type
			"$lf_rt_dom_start\t".        # Start
			"$lf_rt_dom_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;


		    if ($lf_rt_orf_start) {
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "rt_orf\t".                  # Data type
			    "$lf_rt_orf_start\t".        # Start
			    "$lf_rt_orf_end\t".          # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
			
			# Extract sequence
			if ($do_feat_seq) {
			    my $feat_strand;
			    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
				$lf_ltr_id."_rt_orf"; 
			    my $out_seq_path = $ltrfinder_dir.
				$seq_header.".fasta";
			    
			    print ">$seq_header\n";
			    
			    if ($lf_strand =~ "-") {
				$feat_strand="-1";
			    }
			    else {
				$feat_strand="1";
			    }
			    
			    my $loc = Bio::Location::Simple->
				new(-start  => $lf_rt_orf_start,
				    -end    => $lf_rt_orf_end,
				    -strand => $feat_strand);
			    # Sequene of the feature
			    my $feat_seq = $qry_seq->subseq($loc);
			    
			    open (SEQOUT,">$out_seq_path") ||
				die "Can not open:\n$out_seq_path\n";
			    print SEQOUT ">$seq_header\n";
			    print SEQOUT "$feat_seq\n";
			    close SEQOUT;
			}
			# END EXTRACT SEQUENCE
		    } # End has ORF

		} # End of has_reverse_transcriptase
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		   


		# PRINT LINE ENDING TO SUMMARY OUTFILE
		# AND CLOSE THE FILE
		#print SUMOUT "\n";
		print SUMSEQ "\n";

	    } # END OF PRINT OUTPUT


	    #-----------------------------+
	    # RESET VARS TO NULL          |
	    #-----------------------------+
	    # May not need to reset vars since existence of the vars
	    # is designated by the $do_* variables.
	    # These vars are extracted from the LTR Finder output


	    #-----------------------------+
	    # LOAD ID VAR                 |
	    #-----------------------------+
	    $lf_ltr_id = $1;

	} # END OF PRINT STORED VARS
	
	# SEQ ID AND LENGTH
	elsif (m/>Sequence: (.*) Len:(.*)/){
	    $lf_seq_id = $1;
	    $lf_seq_len = $2;
	}
	
	# SPAN LOCATION, LENGTH, AND STRAND
	elsif (m/^Location : (\d*) - (\d*) Len: (\d*) Strand:(.)/){
	    $lf_span_start = $1;
	    $lf_span_end = $2;
	    $lf_length = $3;
	    $lf_strand = $4;
	}
	
	# SCORE SIMILARITY
	elsif (m/^Score    : (.*) \[LTR region similarity:(.*)\]/){
	    $lf_score = $1;
	    $lf_ltr_similarity = $2;
	}
	
	# STATUS SET
	elsif (m/^Status   : (\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)/){
	    # Since this is a binary string, it can be split as digits
	    # and used to load the $has_* booleans
	    $has_5ltr_tg = $1;
	    $has_5ltr_ca = $2;
	    $has_3ltr_tg = $3;
	    $has_3ltr_ca = $4;
	    $has_tsr = $5;
	    $has_pbs = $6;
	    $has_ppt = $7;
	    $has_rt = $8;
	    $has_in_core = $9;
	    $has_in_cterm = $10;
	    $has_rh = $11;

	}
	
	# 5' LTR
	elsif (m/^5\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	    $lf_5ltr_start = $1;
	    $lf_5ltr_end = $2;
	    $lf_5ltr_len = $3;
	}
	
	# 3' LTR
	elsif (m/^3\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	    $lf_3ltr_start = $1;
	    $lf_3ltr_end = $2;
	    $lf_3ltr_len = $3;
	}
	
    # TARGET SITE REPLICATION
	elsif (m/TSR      : (\d*) - (\d*) , (\d*) - (\d*) \[(.*)\]/){
	    $lf_5tsr_start = $1;
	    $lf_5tsr_end = $2;
	    $lf_3tsr_start = $3;
                   	    $lf_3tsr_end = $4;
	    $lf_tsr_string = $5;
	}
	
	# SHARPNESS METRIC
	elsif (m/^Sharpness: (.*),(.*)/){
	    $lf_sharp_5 = $1;
	    $lf_sharp_3 = $2;
	}
	
	# PBS
	elsif (m/PBS   : \[(\d*)\/(\d*)\] (\d*) - (\d*) \((.*)\)/) {
	    $lf_pbs_num_match = $1;
	    $lf_pbs_aln_len = $2;
	    $lf_pbs_start = $3;
	    $lf_pbs_end = $4;
	    $lf_pbs_trna = $5;
	}
	
	# PPT
	elsif (m/PPT   : \[(\d*)\/(\d*)\] (\d*) - (\d*)/) {
	    $lf_ppt_num_match = $1;
	    $lf_ppt_aln_len = $2;
	    $lf_ppt_start = $3;
	    $lf_ppt_end = $4;
	}
	
	# PROTEIN DOMAINS
	# This will need to be modified and checked after another run
	# using ps_scan to get the additional domains
	#
	#Domain: 56796 - 57326 [possible ORF:56259-59144, (IN (core))]
	elsif (m/Domain: (\d*) - (\d*) \[possible ORF:(\d*)-(\d*), \((.*)\)\]/) {
	    
	    $lf_domain_dom_start = $1;
	    $lf_domain_dom_end = $2;
	    $lf_domain_orf_start = $3;
	    $lf_domain_orf_end = $4;
	    $lf_domain_name = $5;
	    
	    # Temp while I work with this data
	    #print "DOMAIN:".$lf_domain_name."\n";
	    
	    if ($lf_domain_name =~ 'IN \(core\)') {
		
		$lf_in_core_dom_start = $lf_domain_dom_start;
		$lf_in_core_dom_end = $lf_domain_dom_end;
		$lf_in_core_orf_start = $lf_domain_orf_start;
		$lf_in_core_orf_end = $lf_domain_orf_end;
		
		# Temp while debug
		# This is to check the regexp vars can be fetched here
		#print "\tDom Start: $lf_in_core_dom_start\n";
		#print "\tDom End:   $lf_in_core_dom_end\n";
		#print "\tOrf Start: $lf_in_core_orf_start\n";
		#print "\tOrf End:   $lf_in_core_orf_end\n";
		
	    }
	    elsif ($lf_domain_name =~ 'IN \(c-term\)') {
		$lf_in_cterm_dom_start = $lf_domain_dom_start;
		$lf_in_cterm_dom_end = $lf_domain_dom_end;
		$lf_in_cterm_orf_start = $lf_domain_orf_start;
		$lf_in_cterm_orf_end = $lf_domain_orf_end;
	    }
	    elsif ($lf_domain_name =~ 'RH') {
		$lf_rh_dom_start = $lf_domain_dom_start;
		$lf_rh_dom_end = $lf_domain_dom_end;
		$lf_rh_orf_start = $lf_domain_orf_start;
		$lf_rh_orf_end = $lf_domain_orf_end;
	    }
	    elsif ($lf_domain_name =~ 'RT') {
		
		$lf_rt_dom_start = $lf_domain_dom_start;
		$lf_rt_dom_end = $lf_domain_dom_end;
		$lf_rt_orf_start = $lf_domain_orf_start;
		$lf_rt_orf_end = $lf_domain_orf_end;  
		
	    }
	    # RNASE H BY ANOTHER NAME
	    elsif ($lf_domain_name =~ 'RNase H') {
		$lf_rh_dom_start = $lf_domain_dom_start;
		$lf_rh_dom_end = $lf_domain_dom_end;
		$lf_rh_orf_start = $lf_domain_orf_start;
		$lf_rh_orf_end = $lf_domain_orf_end;
	    }
	    else {
		print "\a";
		print STDERR "Unknown domain type: $lf_domain_name\n";
	    }
	    
	    
	    
	} # End of elsif Domain
	
	#-----------------------------+
	# FILE HEADER INFORMATION     |
	#-----------------------------+
	
	# PROGRAM NAME
	elsif (m/^Program    : (.*)/) {
	    $lf_prog_name = $1;
	}
	
	# PROGRAM VERSION
	elsif (m/^Version    : (.*)/) {
	    $lf_version = $1;
	}
	
    }
    
    close INFILE;
    
    
    #//////////////////////////////
    # BEGIN OF PRINT LAST RECORD
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    #-----------------------------+
    # CHECK TO SEE IF IT HAD HITS |
    #-----------------------------+
    if ($has_ltr_retro) {
	
		#-----------------------------+
		# SUMMARY OUTPUT              |
		#-----------------------------+
		print SUMOUT $seq_id."_LTRFinder_".$gff_suffix.$lf_ltr_id."\t".
		    $lf_span_start."\t".     # 2
		    $lf_span_end."\t".       # 3
		    $lf_length."\t".         # 4
		    $lf_score."\t".          # 5
		    # LTR INFO
		    $lf_ltr_similarity."\t". # 6
		    $lf_5ltr_len."\t".       # 7
		    $lf_3ltr_len."\t".       # 8
		    $lf_sharp_5."\t".        # 9
		    $lf_sharp_3."\t".        # 10
		    # BOOLEANS
		    $has_5ltr_tg."\t".       # 11
		    $has_5ltr_ca."\t".       # 12
		    $has_3ltr_tg."\t".       # 13
		    $has_3ltr_ca."\t".       # 14
		    $has_tsr."\t".           # 15
		    $has_pbs."\t".           # 16
		    $has_ppt."\t".           # 17
		    $has_rt."\t".            # 18
		    $has_in_core."\t".       # 19
		    $has_in_cterm."\t".      # 20
		    $has_rh.                 # 21
		    "\n";                    # NEWLINE
		    
		
		print SUMSEQ $seq_id."_LTRFinder_".$gff_suffix.$lf_ltr_id."\t";

		# Two alternatives here
		# keep appending to gff_lout or LTR Retro object 
		# variable and then print set
		# or print line at a time .. currently printing 
		# a line at a time. JCE 10/02/2007
		# Advantage to this approach is that it is not memory
		# limited. Extracted LTR sequences for entire chromosomes
		# could run into mememory limites using an object
		# oriented approach.

		#-----------------------------+
		# FULL SPAN OF LTR RETRO      |
		#-----------------------------+
		$gff_str_out = "$lf_seq_id\t".   # Seq ID
		    "$gff_src\t".                # Source
		    "LTR_retrotransposon\t".     # Data type
		    "$lf_span_start\t".          # Start
		    "$lf_span_end\t".            # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
		# Extract sequence
		if ($do_feat_seq) {
		    my $feat_strand;
		    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			$lf_ltr_id."_LTR_retrotransposon"; 
		    my $out_seq_path = $ltrfinder_dir.$seq_header.".fasta";
		    
		    print ">$seq_header\n";
		    
		    if ($lf_strand =~ "-") {
			$feat_strand="-1";
		    }
		    else {
			$feat_strand="1";
		    }
		    
		    # TEMP PRINT OF FEATURE STRAND
		    print STDERR "STRAND: $feat_strand\n" if $verbose;
		    
		    my $loc = Bio::Location::Simple->
			new(-start  => $lf_span_start,
			    -end    => $lf_span_end,
			    -strand => $feat_strand);
		    # Sequene of the feature
		    my $feat_seq = $qry_seq->subseq($loc);
		    
		    open (SEQOUT,">$out_seq_path") ||
			die "Can not open\n:$out_seq_path\n";
		    print SEQOUT ">$seq_header\n";
		    print SEQOUT "$feat_seq\n";
		    close SEQOUT;

		    print SUMSEQ "$feat_seq\t";
		}
		# END EXTRACT SEQUENCE
		
		#-----------------------------+
		# PRIMER BINDING SITE         |
		#-----------------------------+
		if ($has_pbs) {
		    $gff_str_out = "$lf_seq_id\t".  # Seq ID
			"$gff_src\t".               # Source
			"primer_binding_site\t".    # Data type
			"$lf_pbs_start\t" .         # Start
			"$lf_pbs_end\t".            # End
			"$lf_score\t".              # Score
			"$lf_strand\t".             # Strand
			".\t".                      # Frame
			"ltr_finder_$lf_ltr_id\n";  # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    # Extract sequence
		    if ($do_feat_seq) {
			my $feat_strand;
			my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			    $lf_ltr_id."_primer_binding_site"; 
			my $out_seq_path = $ltrfinder_dir.
			    $seq_header.".fasta";
			
			print ">$seq_header\n";
			
			if ($lf_strand =~ "-") {
			    $feat_strand="-1";
			}
			else {
			    $feat_strand="1";
			}
			
			my $loc = Bio::Location::Simple->
			    new(-start  => $lf_pbs_start,
				-end    => $lf_pbs_end,
				-strand => $feat_strand);
			
			# Sequene of the feature
			my $feat_seq = $qry_seq->subseq($loc);
			
			open (SEQOUT,">$out_seq_path") ||
			    die "Can not open\n:$out_seq_path\n";
			print SEQOUT ">$seq_header\n";
			print SEQOUT "$feat_seq\n";
			close SEQOUT;
			
			print SUMSEQ "$feat_seq\t"; 
			
		    }
		    # END EXTRACT SEQUENCE
		} #
		else {
		    # Case where we want to extract seq features but
		    # PBS is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		
		#-----------------------------+
		# POLYPURINE TRACT            |
		#-----------------------------+
		if ($has_ppt) {
		    $gff_str_out = "$lf_seq_id\t".  # Seq ID
			"$gff_src\t".               # Source
			"RR_tract\t".               # Data type
			"$lf_ppt_start\t".          # Start
			"$lf_ppt_end\t".            # End
			"$lf_score\t".              # Score
			"$lf_strand\t".             # Strand
			".\t".                      # Frame
			"ltr_finder_$lf_ltr_id\n";  # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    # Extract sequence
		    if ($do_feat_seq) {
			my $feat_strand;
			my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			    $lf_ltr_id."_RR_tract"; 
			my $out_seq_path = $ltrfinder_dir.
			    $seq_header.".fasta";
			
			print ">$seq_header\n";
			
			if ($lf_strand =~ "-") {
			    $feat_strand="-1";
			}
			else {
			    $feat_strand="1";
			}
			
			my $loc = Bio::Location::Simple->
			    new(-start  => $lf_ppt_start,
				-end    => $lf_ppt_end,
				-strand => $feat_strand);
			# Sequene of the feature
			my $feat_seq = $qry_seq->subseq($loc);
			
			open (SEQOUT,">$out_seq_path") ||
			    die "Can not open :\n$out_seq_path\n";
			print SEQOUT ">$seq_header\n";
			print SEQOUT "$feat_seq\n";
			close SEQOUT;
		    }
		    # END EXTRACT SEQUENCE
		}
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		 
		
		#-----------------------------+
		# 5' LTR                      |
		#-----------------------------+
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "five_prime_LTR\t".         # Data type
		    "$lf_5ltr_start\t".         # Start
		    "$lf_5ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
		
		# Extract sequence
		if ($do_feat_seq) {
		    my $feat_strand;
		    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			$lf_ltr_id."_five_prime_LTR"; 
		    my $out_seq_path = $ltrfinder_dir.$seq_header.".fasta";
		    
		    print ">$seq_header\n";
		    
		    if ($lf_strand =~ "-") {
			$feat_strand="-1";
		    }
		    else {
			$feat_strand="1";
		    }
		    
		    my $loc = Bio::Location::Simple->new
			(-start  => $lf_5ltr_start,
			 -end    => $lf_5ltr_end,
			 -strand => $feat_strand);
		    # Sequene of the feature
		    my $feat_seq = $qry_seq->subseq($loc);
		    
		    open (SEQOUT,">$out_seq_path") ||
			die "Can not open :\n$out_seq_path\n";
		    print SEQOUT ">$seq_header\n";
		    print SEQOUT "$feat_seq\n";
		    close SEQOUT;

		    # PRINT TO SUMMARY SEQUENCE FILE
		    print SUMSEQ "$feat_seq\t";
		}
		# END EXTRACT SEQUENCE
		
		#-----------------------------+
		# 3' LTR                      |
		#-----------------------------+
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "three_prime_LTR\t".        # Data type
		    "$lf_3ltr_start\t".         # Start
		    "$lf_3ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
		
		# Extract sequence
		if ($do_feat_seq) {
		     my $feat_strand;
		     my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			 $lf_ltr_id."_three_prime_LTR"; 
		     my $out_seq_path = $ltrfinder_dir.$seq_header.".fasta";
		     
		     print ">$seq_header\n";
		     
		     if ($lf_strand =~ "-") {
			 $feat_strand="-1";
		     }
		     else {
			 $feat_strand="1";
		     }
		     
		     my $loc = Bio::Location::Simple->
			 new(-start  => $lf_3ltr_start,
			     -end    => $lf_3ltr_end,
			     -strand => $feat_strand);
		     # Sequene of the feature
		     my $feat_seq = $qry_seq->subseq($loc);
		     
		     open (SEQOUT,">$out_seq_path") ||
			 die "Can not open:\n$out_seq_path\n";
		     print SEQOUT ">$seq_header\n";
		     print SEQOUT "$feat_seq\n";
		     close SEQOUT;

		     
		     # PRINT TO SUMMARY SEQUENCE FILE
		     print SUMSEQ "$feat_seq\t";
		 }
		# END EXTRACT SEQUENCE
		
		#-----------------------------+
		# TARGET SITE DUPLICATION     |
		#-----------------------------+
		if ($has_tsr) {
		    
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"target_site_duplication\t". # Data type
			"$lf_5tsr_start\t".          # Start
			"$lf_5tsr_end\t".            # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;	    
		    
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"target_site_duplication\t". # Data type
			"$lf_3tsr_start\t".          # Start
			"$lf_3tsr_end\t".            # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    # Extract sequence
		    if ($do_feat_seq) {
			my $feat_strand;
			my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
			    $lf_ltr_id."_target_site_duplication"; 
			my $out_seq_path = $ltrfinder_dir.
			    $seq_header.".fasta";
			
			print ">$seq_header\n";
			
			if ($lf_strand =~ "-") {
			    $feat_strand="-1";
			}
			else {
			    $feat_strand="1";
			}
			
			my $loc = Bio::Location::Simple->
			    new(-start  => $lf_5tsr_start,
				-end    => $lf_5tsr_end,
				-strand => $feat_strand);
			 # Sequene of the feature
			 my $feat_seq = $qry_seq->subseq($loc);
			 
			 open (SEQOUT,">$out_seq_path") ||
			     die "Can not open\n:$out_seq_path\n";
			 print SEQOUT ">$seq_header\n";
			 print SEQOUT "$feat_seq\n";
			 close SEQOUT;
		     }
		     # END EXTRACT SEQUENCE
		 }
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		    
		    

		#-----------------------------+
		# INTEGRASE                   |
		#-----------------------------+
		if ($has_in_core) {
		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"integrase_core_domain\t".   # Data type
			"$lf_in_core_dom_start\t".   # Start
			"$lf_in_core_dom_end\t".     # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    
		    # ONLY DO THE FOLLOWING IF AND INTEGRASE
		    # CORE ORF WAS FOUND. It is possible that a
		    # core was found but a full orf can not be extracted
		    if ($lf_in_core_orf_start) {
			#/////////
			# NOT SONG
			#\\\\\\\\\
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "integrase_core_orf\t".      # Data type
			    "$lf_in_core_orf_start\t".   # Start
			    "$lf_in_core_orf_end\t".     # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
			
			# Extract sequence
			if ($do_feat_seq) {
			    my $feat_strand;
			    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
				$lf_ltr_id."_in_core_orf"; 
			    my $out_seq_path = $ltrfinder_dir.
				$seq_header.".fasta";
			    
			    print ">$seq_header\n";
			    
			    if ($lf_strand =~ "-") {
				$feat_strand="-1";
			    }
			    else {
				$feat_strand="1";
			    }
			    
			    my $loc = Bio::Location::Simple->new
				(-start  => $lf_in_core_orf_start,
				 -end    => $lf_in_core_orf_end,
				 -strand => $feat_strand);
			    # Sequene of the feature
			    my $feat_seq = $qry_seq->subseq($loc);
			    
			    open (SEQOUT,">$out_seq_path") ||
				die "Can not open:\n$out_seq_path\n";
			    print SEQOUT ">$seq_header\n";
			    print SEQOUT "$feat_seq\n";
			    close SEQOUT;
			}
			# END EXTRACT SEQUENCE
		    }
		    else {
			# Case where we want to extract seq features but
			# PPT is not present
			if ($do_feat_seq) {
			    print SUMSEQ ".\t"; 
			}
		    }
		   

		    
		} # End of has in_core
		
		
		#-----------------------------+
		# GAG C-TERMINAL END          |
		#-----------------------------+
		if ($has_in_cterm) {
		    
		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"integrase_cterm_domain\t".  # Data type
			"$lf_in_cterm_dom_start\t".  # Start
			"$lf_in_cterm_dom_end\t".    # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    if ($lf_in_cterm_orf_start) {
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "integrase_cterm_orf\t".     # Data type
			    "$lf_in_cterm_orf_start\t".  # Start
			    "$lf_in_cterm_orf_end\t".    # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
		    }

		} # End of has_in_cterm
		

		#-----------------------------+
		# RNASEH                      |
		#-----------------------------+
		if ($has_rh) {
		    
		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"rnaseh_domain\t".           # Data type
			"$lf_rh_dom_start\t".        # Start
			"$lf_rh_dom_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;
		    
		    if ($lf_rh_orf_start) {
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "rnaseh_orf\t".              # Data type
			    "$lf_rh_orf_start\t".        # Start
			    "$lf_rh_orf_end\t".          # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
			
			# Extract sequence
			if ($do_feat_seq) {
			    my $feat_strand;
			    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
				$lf_ltr_id."_rnaseh_orf"; 
			    my $out_seq_path = $ltrfinder_dir.
				$seq_header.".fasta";
			    
			    print ">$seq_header\n";
			    
			    if ($lf_strand =~ "-") {
				$feat_strand="-1";
			    }
			    else {
				$feat_strand="1";
			    }
			    
			    my $loc = Bio::Location::Simple->
				new(-start  => $lf_rh_orf_start,
				    -end    => $lf_rh_orf_end,
				    -strand => $feat_strand);
			    # Sequene of the feature
			    my $feat_seq = $qry_seq->subseq($loc);
			    
			    open (SEQOUT,">$out_seq_path") ||
				die "Can not open $out_seq_path\n";
			    print SEQOUT ">$seq_header\n";
			    print SEQOUT "$feat_seq\n";
			    close SEQOUT;
			}
			# END EXTRACT SEQUENCE
		    }
		}
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		
		#-----------------------------+
		# REVERSE TRANSCRIPTASE       |
		#-----------------------------+
		if ($has_rt) {

		    #/////////
		    # NOT SONG
		    #\\\\\\\\\
		    $gff_str_out = "$lf_seq_id\t".   # Seq ID
			"$gff_src\t".                # Source
			"rt_domain\t".               # Data type
			"$lf_rt_dom_start\t".        # Start
			"$lf_rt_dom_end\t".          # End
			"$lf_score\t".               # Score
			"$lf_strand\t".              # Strand
			".\t".                       # Frame
			"ltr_finder_$lf_ltr_id\n";   # Retro ID
		    print GFFOUT $gff_str_out;


		    if ($lf_rt_orf_start) {
			$gff_str_out = "$lf_seq_id\t".   # Seq ID
			    "$gff_src\t".                # Source
			    "rt_orf\t".                  # Data type
			    "$lf_rt_orf_start\t".        # Start
			    "$lf_rt_orf_end\t".          # End
			    "$lf_score\t".               # Score
			    "$lf_strand\t".              # Strand
			    ".\t".                       # Frame
			    "ltr_finder_$lf_ltr_id\n";   # Retro ID
			print GFFOUT $gff_str_out;
			
			# Extract sequence
			if ($do_feat_seq) {
			    my $feat_strand;
			    my $seq_header = $seq_id."_LTRFinder_".$gff_suffix.
				$lf_ltr_id."_rt_orf"; 
			    my $out_seq_path = $ltrfinder_dir.
				$seq_header.".fasta";
			    
			    print ">$seq_header\n";
			    
			    if ($lf_strand =~ "-") {
				$feat_strand="-1";
			    }
			    else {
				$feat_strand="1";
			    }
			    
			    my $loc = Bio::Location::Simple->
				new(-start  => $lf_rt_orf_start,
				    -end    => $lf_rt_orf_end,
				    -strand => $feat_strand);
			    # Sequene of the feature
			    my $feat_seq = $qry_seq->subseq($loc);
			    
			    open (SEQOUT,">$out_seq_path") ||
				die "Can not open:\n$out_seq_path\n";
			    print SEQOUT ">$seq_header\n";
			    print SEQOUT "$feat_seq\n";
			    close SEQOUT;
			}
			# END EXTRACT SEQUENCE
		    } # End has ORF

		} # End of has_reverse_transcriptase
		else {
		    # Case where we want to extract seq features but
		    # PPT is not present
		    if ($do_feat_seq) {
			print SUMSEQ ".\t"; 
		    }
		}
		   
		# PRINT LINE ENDING TO SUMMARY OUTFILE
		# AND CLOSE THE FILE
		#print SUMOUT "\n";
		print SUMSEQ "\n";

    } # End of if sequence has LTR_retro

    #////////////////////////////////
    # END OF PRINT LAST RECORD
    #\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    close (SUMOUT);
    close (SUMSEQ);
    
} # End of ltrfinder2gff subfunction


1;
__END__

# Old print_help subfunction

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"MyProg.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --quiet        # Run program with minimal output\n";
	
    if ($opt =~ "full") {
	print "\n$usage\n\n";
	print "$args\n\n";
    }
    else {
	print "\n$usage\n\n";
    }
    
    exit;
}

=head1 NAME

batch_ltrfinder.pl - Run the LTRFinder program in batch mode

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    batch_ltrfinder.pl -i InDir -o OutDir -c ConfigFile

=head2 Required Arguments

    -i,--indir    # Directory of fasta files to process
    -o,--outdir   # Path to the base output directory
    -c,--config   # Path to the config file

=head1 DESCRIPTION

This program will run the LTR_FINDER program for a set of fasta files.
You can set multiple parameter sets using a config file.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c,--config

Path to the batch_ltrfinder config file. This is a tab delimited text file 
made up of two columns. The first column indicates the name assigned to
the parameter set, while the second column contains the flags that
will be passed to the LTR_FINDER program.

=back

=head1 OPTIONS

=over 2

=item --ltr-finder

Path to the LTR_FINDER binary.

=item -s,--trna-db

Path to the tRNA database used by LTR_FINDER.
This file is part of the LTR_FINDER download.

=item -a,--prosite

Path to the prosite directory for use by LTR_FINDER.

=item -g,gff

Convert the outout to gff format.

=item -f,--feat-seq

Extract sequence string of features. These will be extracted as separate fasta
files for each sequence feature class. These fasta files will be stored in the
ltr_finder directory. Using rootname as the name of the contig and ltrid
as the number assigned by LTR_Finder, the fasta files that are created are:

=over

=item rootname_ltrid_ltr5.fasta

The sequence of the five prime LTR. This is not extracted from the alignment
returned by LTR_Finder, but is extracted from the original sequence using
the coordinates returned by LTR finder.

=item rootname_ltr3.fasta

The squence of the three primer LTR. This is also extracted from the 
original sequence using the coordinates returned by LTR finder.

=item rootname_ltr

=back

The name of the individual LTR retrotransposon identified will be 
consistant across these fasta files. Not all 

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item --verbose

Run the program in verbose mode.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

The list of error messages that can be generated,
explanation of the problem
one or more causes
suggested remedies
list exit status associated with each error

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

The location of the configuration file is indicated by the --config option
at the command line.
This is a tab delimited text file allowing you to run the LTR_FINDER
program with multiple parameter sets for each query sequence file.
This config file is a two column file. Lines beginning with the # symbol
are ignored and provide a way to add comments to your config file. This config
file is expected to have UNIX format line endings. You should therefore
avoid creating config files using Windows based programs such as MS Word.

=over 2

=item Col 1. Parameter set name

This is the name assigned to the parameter set defined on the current
line. This name must not have any spaces. This name will be appended to
the output gff file.

=item Col 2. LTR_FINDER Options

A string of command line options to send to the LTR finder program. For
a full set of possible variables, see the LTR_FINDER User Manual 
http://tlife.fudan.edu.cn/ltr_finder/help/single.html

=back

Example config file:

   # Simple batch_ltrfinder config file
   #
   def
   p_30	-p 30
   p_10	-p 10
   # END

=head2 User Environment

This program makes use of the following variables defined in the
user's environment.

=over 2

=item TRNA_DB

This is the path to TRNA_DB used by LTR_FINDER.

=item PROSITE_DIR

This is the path to the directory of Prosite models use by LTR_FINDER.

=item LTR_FINDER

This is the path to the LTR_FINDER binary.

=back

Example environment variables set in the bash shell:

   export TRNA_DB='/home/yourname/apps/LTR_Finder/tRNAdb/Os-tRNAs.fa'
   export PROSITE_DIR='/home/yourname/Apps/LTR_Finder/ps_scan'
   export LTR_FINDER='ltr_finder'

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * LTR_FINDER

The LTR_FINDER program can be obtained as a Linux binary by contacting
the author: xuzh <at> fudan.edu.cn

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the BLAST results.

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 SEE ALSO

The batch_ltrfinder.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Config file must use UNIX format line endings

The config file must have UNIX formatted line endings. Because of
this any config files that have been edited in programs such as
MS Word must be converted to a UNIX compatible text format before
being used with batch_blast.

=back

=head1 REFERENCE

A manuscript is being submitted describing the DAWGPAWS program. 
Until this manuscript is published, please refer to the DAWGPAWS 
SourceForge website when describing your use of this program:

JC Estill and JL Bennetzen. 2009. 
The DAWGPAWS Pipeline for the Annotation of Genes and Transposable 
Elements in Plant Genomes.
http://dawgpaws.sourceforge.net/

=head1 LICENSE

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 10/02/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 10/02/2007
# - Started program from general template
# - Added the ltrfinder2gff subfunction from the 
#   cnv_ltrfinder2gff program.
#
# 12/12/2007
# - Added SVN Rev to propset
# - Updated POD documentation
# - Updated print_help to subfunction that extracts
#   usage and help statements from POD documentation
# - Added an example config file to the config dir
# - Added some default locations of vars outside of base
#   base ENV vars. These assume that that the relevant
#   files have default names and are in the user's path.
#   
# 01/20/2008
# - Removed STDOUT print of gff strings
# - Adding fasta output of sequence features
#   This uses the -f,--feat-seq command line options
#
# 01/22/2008
# - Working version of extract sequence features added.
#   Does not work for the last predicted LTR retro in the
#   LTR Finder output, but will produce sequence for the
#   following features
#    ----------------------------------------------------
#    FEATURE                     SONG NAME
#    ----------------------------------------------------
#    LTR Retrotransposon span  .LTR_retrotransposon
#    Five Prime LTR ............five_prime_LTR
#    Three Prime LTR ...........three_prime_LTR
#    Target Site Duplication....target_site_duplication
#    Polypurine Tract ..........RR_tract
#    Primer Binding Site .......primer_binding_site
#    Integrase Core ORF ........in_core_orf
#   There are not SONG complient names for protein coding
#   sequence features or the C terminal end of gag.
#
# 01/23/2008
# - Cleaning up sequence extraction code
# - Adding ability to extract sequence for last seq in the
#
# 01/25/2008
# - Added check that sequence actually has LTR_retro before
#   trying to print GFF and extract sequences of features
# - It is possible that a CORE to the coding sequence was found
#   but a full length ORF could not be found. In those cases it 
#   is necessary to check for the existence of an ORF before
#   printing to output. Updating code to take this into account
# - Adding RNase H has a domain type
# - Can extract root_name name from *.masked.fasta
# - root_name default is full sequence name if other root not
#   found. This is better then setting value to UND
# - Added more information to status message:
#    - process number of total number of processes
#    - sequence number of total number of sequences
#    - parameter set number of total number of par sets
#
# 01/30/2008
# - Extracting additional data from LTR_Finder output
#    - lf_summary_path
#    - This is in the form of a tab deliminted text file that can
#      be parse into a database
#    - This will be created by default
#  
