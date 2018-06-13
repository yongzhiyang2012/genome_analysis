#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_ltrseq.pl - Run the LTR_seq program in batch mode   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/06/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

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
my $indir;
my $outdir;
my $file_config;

# LTR SEQ DIR

# This is the dir that will hold the config files

my $ltrseq_dir = $ENV{LTR_SEQ_DIR} || $ENV{HOME};

my $ltrseq_bin = $ENV{LTR_SEQ_BIN} || "LTR_seq";

# Booleans
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $verbose = 0;
my $do_test = 0; 

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions("i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # OPTIONS
		    "c|config=s"   => \$file_config,
		    "config-dir=s" => \$ltrseq_dir,
		    "ltrseq-bin=s" => \$ltrseq_bin,
		    "test"         => \$do_test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "verbose"      => \$verbose,
		    "q|quiet"      => \$quiet,);

#-----------------------------+
# SHOW REQUESTED HELP         |
#-----------------------------+
if ( $show_usage ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\nbatch_ltrseq.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}


#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print_help ("usage", $0 );
}


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

unless ($ltrseq_dir =~ /\/$/ ) {
    $ltrseq_dir = $ltrseq_dir."/";
}

#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+ 

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# GET INFO FROM CONFIG FILE   |
#-----------------------------+
my @config_opts;
if ($file_config) {

    my $i;
    # Load to the 2d config_opts array
    open (CONFIG, $file_config) ||
	die "Can't open the config file:\n$file_config";
    $i = 0;
    my $line_num = 0;
    
    while (<CONFIG>) {
	$line_num ++;
	
	next if m/^\#/;
	chomp;
	my @tmpary = split (/\t/);
	my $count_tmp = @tmpary;
	
	if ($count_tmp == 2) {
	    $config_opts[$i][0] = $tmpary[0];  # Config name
	    # The following assumes the full path is in the config file
	    #$config_opts[$i][1] = $tmpary[1];  # Config file
	    $config_opts[$i][1] = $ltrseq_dir.$tmpary[1];  # Config file path
	    $i++;
	} 
	elsif ($count_tmp == 1) {
	    # Problem, this may also happen if the config file has not
	    # been delimited by tabs.
	    # use default options for ltr_seq if only a name is given
	    $config_opts[$i][0] = $tmpary[0];  # Config name
	    $config_opts[$i][1] = "DEF";       # Tag as default config
	    $i++;   
	}
	else {
	    print STDERR "ERROR: Config file line number $line_num\n";
	    print STDERR "       $line_num columns of data were found\n";
	    print STDERR "       Two columns were expected\n";
	}

	
    } # End of while CONFIG
    close CONFIG;
}

#-----------------------------+
# If no config file use the   |
# default options tag         |
#-----------------------------+
else {
    $config_opts[0][0] = "default";   # Config name
    $config_opts[0][1] = "LTR.cfg";   # Set path to LTR.cfg
}

#-----------------------------+
# CHECK TO SEE IF THE LTRSEQ  |
# CONFIG FILE EXISTS          |
#-----------------------------+
for my $config_opt (@config_opts) {
    
    my $ltrseq_cfg_file = $config_opt->[1];
    
    unless ($ltrseq_cfg_file =~ "DEF") {
	unless (-e $ltrseq_cfg_file) {
	    print "\a";
	    die "LTR_Seq config file not found at:\n$ltrseq_cfg_file\n";
	} # End of check to see if the config files exists
    } # End of unless we are using the defualt parameter set
    
} # For each option set 

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

for my $ind_file (@fasta_files) {

    my $name_root;
    
    #-----------------------------+
    # Get the root name of the    |
    # file to mask                |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
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
	    die "ERROR: Could not create root name dir:\n$name_root_dir\n";
    }
    
    #-----------------------------+
    # CREATE THE LTRSEQ DIR       |
    #-----------------------------+
    my $ltrseq_dir = $name_root_dir."ltrseq/";
    unless (-e $ltrseq_dir) {
	mkdir $ltrseq_dir ||
	    die "Could not create ltrseq out dir:\n$ltrseq_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # This will hold the gff file modified from the
    # original gff fine
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "ERROR: Could not create GFF out dir:\n$gff_dir\n";
    }

    #-----------------------------+
    # RUN LTR SEQ                 |
    #-----------------------------+
    # Do this for every config set in the config file

    # Note: LTR_seq usage
    #LTR_seq {Fasta File} {#Sequences} [LTR.cfg file location optional] 
    
    for my $config_opt (@config_opts) {

	my $config_name = $config_opt->[0];
	my $ltrseq_cfg_file = $config_opt->[1];

	my $seq_file = "$indir$ind_file";
	
	# There should always be a single sequence in the input fasta file
	# therefore there is a one in the command syntax
	my $ltrseq_cmd = "$ltrseq_bin $seq_file 1";
	
	# Append config file to command if given
	unless ($ltrseq_cfg_file =~ "DEF") {
	    $ltrseq_cmd = $ltrseq_cmd." ".$ltrseq_cfg_file;
	}
	else {
	    $ltrseq_cmd = $ltrseq_cmd." LTR.cfg";
	}
	
	my $ltrseq_outfile = $ltrseq_dir.$name_root."_ltrseq_".
	    $config_name.".out";
	my $gff_outfile = $gff_dir.$name_root."_ltrseq_".
	    $config_name.".gff";
	
	$ltrseq_cmd = $ltrseq_cmd." > ".$ltrseq_outfile;

	print STDERR "Processing: $ltrseq_cmd\n" if $verbose;

	system ($ltrseq_cmd) unless $do_test;;
	
	#-----------------------------+
	# CONVERT OUTPUT TO GFF       |
	#-----------------------------+
        #ltrseq2gff ($ltrseq_in, $gff_out, $append_gff, $seqname);
	if (-e $ltrseq_outfile) {
	    ltrseq2gff ($ltrseq_outfile, $gff_outfile, 0, 
			$name_root, $config_name);
	}
	else {
	    # Report error, can not translate to gff
	}
    }
    
}

exit 0;

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
sub ltrseq2gff {

    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    my ($ltrseq_in, $gff_out, $append_gff, $seqname, $param_name) = @_;

    #-----------------------------+
    # LTRSEQ VARS                 |
    #-----------------------------+
    my $ltrseq_infile;              # Input file that was analyzed by LTR_seq
    my $ltrseq_numseqs;             # Number of seqs analyzed by LTR_seq
    my $ltrseq_msrt;                # Max Score Ration Threshold
    my $ltrseq_ltrmin;              # LTRmin
    my $ltrseq_ltr_minexact;        # LTR Min Exact Match
    my $ltrseq_dmin;                # Dmin
    my $ltrseq_dmax;                # Dmax

    # Postion of the LTR Retrotransposon features
    my $ltrseq_name;                # Unique name assigned to the prediction
    my $ltrspan_start;              # Start of the entire LTR Retrotransposon
    my $ltrspan_end;                # End of the entire LTR Retrotransposon
    my $ltrspan_len;                # Length of the LTR span
    my $ltr5_start;                 # Start of the 5' LTR
    my $ltr5_end;                   # End of the 5' LTR
    my $ltr3_start;                 # Start of the 3' LTR
    my $ltr3_end;                   # End of the 3' LTR
    my $ltr_len;                    # Length of the LTRs
    my $mid_start;                  # Start of the LTR Mid region
    my $mid_end;                    # End of the LTR Mid region
    my $ltr_diff;                   # Percent Difference in the LTRs
    my $ltr_conf;                   # Confidence score for the prediction
    my $ltr_tsr;                    # Target site rep

    my @in_split = ();              # Split of the infile line
    my $num_in;                     # Number of split vars in the infile

    # Initialize Counters
    my $ltrseq_num = 0;             # ID Number of putatitve LTR retro 
                                    # (Incremented Number)

    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+
    open (INFILE, "<$ltrseq_in") ||
	die "ERROR: Can not open input file:\n$ltrseq_in\n";

    if ($append_gff) {
	open (GFFOUT, ">>$gff_out") ||
	    die "ERROR: Can not open output file for appending\n$gff_out\n";
    }
    else {
	open (GFFOUT, ">$gff_out") ||
	    die "ERROR: Can not open output file for output\n$gff_out\n";
    } # End of if append_gff

    #-----------------------------+
    # PROCESS INFILE              |
    #-----------------------------+
    while (<INFILE>) {
	chomp;

	# Print the INFILE
	#print "$_\n";

	if (m/^Report:/) {

	    # If this is a Report line, accepted or rejected can be
	    # determined by counting the line number
	    # 15 is a rejected line
	    # 19 is an accepted line 
	    my @in_split = split;
	    my $num_in = @in_split;   

	    #print "Report Line: $num_in parts\n";
	    #print "\t$_\n";

	    if ($num_in == 15) {
		#print "\tREJECTED\n";
	    }
	    #-----------------------------+
	    # ACCEPTED LTR RETRO Parts    | 
	    #-----------------------------+
	    elsif ($num_in == 19) {
		#print "\tACCEPTED\n";
		#print "\n".$_."\n";
		$ltrseq_num++;              # Increment the count

		$ltr5_start = $in_split[4] || "NULL";
		$ltr5_end = $in_split[5] || "NULL";
		$ltr3_start = $in_split[7] || "NULL";
		$ltr3_end = $in_split[8] || "NULL";
		$mid_start = $ltr5_end + 1;
		$mid_end = $ltr3_start - 1;
		$ltrspan_start = $ltr5_start;
		$ltrspan_end = $ltr3_end;
		$ltrspan_len = $in_split[12];
		$ltr_len = $in_split[10];
		$ltr_diff = $in_split[14];
		$ltr_conf = $in_split[16];
		$ltr_tsr = $in_split[18];

		$ltrseq_name = $seqname."_ltrseq_".$ltrseq_num;

		#-----------------------------+
		# GET TSR LOCATIONS           |
		#-----------------------------+
		# The location of the math to generate these will
		# be strand dependent.
		my $has_tsr = 0;
		my $tsr5_start;
		my $tsr5_end;
		my $tsr3_start;
		my $tsr3_end;

		if ($ltr_tsr =~ m/\#(.*)\#/) {
		    $ltr_tsr = $1;
		    my $ltr_tsr_len = length($ltr_tsr);

		    if ($ltr_tsr_len > 0) {
			$has_tsr = 1;
			$tsr5_start = $ltr5_start - $ltr_tsr_len - 1;
			$tsr5_end = $ltr5_start - 1;
			$tsr3_start = $ltr3_end + 1;
			$tsr3_end = $ltr3_end + $ltr_tsr_len  + 1;
		    }
		}

		#-----------------------------+
		# SHOW INFO IF VERBOSE
		#-----------------------------+
		if ($verbose) {
		    print STDERR "$ltrseq_name\n";
		    print STDERR "\tSTART:    $ltrspan_start\n";
		    print STDERR "\tEND:      $ltrspan_end\n";
		    print STDERR "\tLTRLEN:   $ltr_len\n";
		    print STDERR "\tSPAN_LEN: $ltrspan_len\n";
		    print STDERR "\tLTSR:     $ltr_tsr\n";
		    print STDERR "\tLTRCON:   $ltr_conf\n";
		    print STDERR "\tLTRDIF:   $ltr_diff\n";
		} # End of if verbose

		#-----------------------------+
		# PRINT TO GFF OUTPUT FILE    |
		#-----------------------------+
		
		# LTR_seq predicts everything in the positive strand

		# SET THE PROGRAM SOURCE
		my $prog_src;
		if ($param_name) {
		    $prog_src = "LTR_seq:".$param_name;
		}
		else {
		    $prog_src = "LTR_seq";
		}

		# LTR RETRO PREDICTION SPAN
		print GFFOUT "$seqname\t".     # Name of sequence
		    "$prog_src\t".             # Source name
		    "LTR_retrotransposon\t".   # Feature, exon for Apollo
		    "$ltrspan_start\t".        # Start of the ltr span
		    "$ltrspan_end\t".          # End of the ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)
		
		# 5' LTR
		print GFFOUT "$seqname\t".     # Name of sequence
		    "$prog_src\t".             # Source name
		    "five_prime_LTR\t".        # Feature, exon for Apollo
		    "$ltr5_start\t".           # Start of the 5'ltr
		    "$ltr5_end\t".             # End of the 5' ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)

		# 3' LTR
		print GFFOUT "$seqname\t".     # Name of sequence
		    "$prog_src\t".             # Source name
		    "three_prime_LTR\t".       # Feature, exon for Apollo
		    "$ltr3_start\t".           # Start of the 3'ltr
		    "$ltr3_end\t".             # End of the 3' ltr span
		    "$ltr_conf\t".             # LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)

		# REPORT TSRs IF PRESENT AND KNOWN
		if ($has_tsr) {

		    # TSR 5 - 5' Target site duplication
		    print GFFOUT "$seqname\t".     # Name of sequence
			"$prog_src\t".             # Source name
			"target_site_duplication\t". # Feature
			"$tsr5_start\t".           # Start of the 3'ltr
			"$tsr5_end\t".             # End of the 3' ltr span
			"$ltr_conf\t".             # LTR Confidence Score
			".\t".                     # Strand
			".\t".                     # Frame
			"$ltrseq_name\n";          # Features (Name)

		    # TSR - 5' Target site duplication
		    print GFFOUT "$seqname\t".     # Name of sequence
			"$prog_src\t".             # Source name
			"target_site_duplication\t". # Feature
			"$tsr3_start\t".           # Start of the 3'ltr
			"$tsr3_end\t".             # End of the 3' ltr span
			"$ltr_conf\t".             # LTR Confidence Score
			".\t".                     # Strand
			".\t".                     # Frame
			"$ltrseq_name\n";          # Features (Name)
		} # End has_tsr

	    } # end of num parts is 19


	    # 15 Parts are Rejected models

	}
	# HEADER INFORMATION
	elsif (m/^Info:/) {
	    
	}
	# The input string that was passed
	elsif (m/^Input:/) {

	}

    } # End of while INFILE

}


1;
__END__

=head1 NAME

batch_ltrseq.pl - Run the LTR_seq program in batch mode.

=head1 VERSION

This documentation refers to program version Release 1.0.

=head1 SYNOPSIS

=head2 Usage

    batch_ltrseq.pl -i in_dir -o out_dir [-c config]

=head2 Required Options

    -i        # Intput directory of files to process
              # These must be in fasta file format
    -o        # Base output directory

=head1 DESCRIPTION

This program runs the LTR_seq program for each file in the input directory.
An optional configuartion file can also be used to run the LTR_seq program
for multiple LTR_seq parameter combinations. If no configuration file is
provided, batch_ltrseq.pl will run LTR_seq using default parameters.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item -c,--config

Path to a configuration file. This config file  will allow the 
batch_ltrseq.pl program to run the LTR_seq program for multiple parameter 
combinations for every fasta file in the input sequence directory. 
If no configuartion file is used, the LTR_seq program will be run with
default parameters.

The configuration file will
be a two column, tab delimited text file with the following options:

=over

=item * Col. 1 Configuration Set Name

This is the name that will be used in the gff file to tag the configuration
set. For example def for default options or old for an optional set to
try to finder LTR retrotransposon insertions that are old.

=item * Col. 2 Config File

LTR_seq config file describing the parameter set to use. This file
should follow the format used for LTR_Seq config files. If this value is
set to DEF, then batch_ltrseq.p will run the LTR_seq program with
default parameters.

The directory that this config file is located in will be specified by
the LTR_SEQ_DIR environment variable or defined in the command line
with --config-dir.

=back

An example configuration file is shown below:

    def     LTR.cfg
    long    long_ltr.cfg
    old     lod_ltr.cfg

This configuration file would run three sets of parameter files for every
sequence in the input directory.

=item --config-dir

The directory that the configuration files are located in. This can
also be specified by the LTR_SEQ_DIR environment variable.

=item --config-bin

The location of the LTR_seq binary. This can also be specified by the
LTR_SEQ_BIN environment variable.

=item --test

Run the batch_ltrseq.pl program in test mode. This will create the required
directory structure, test for existence of files and other program functions.
However, this will not actually run the LTR_seq program.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

=over 2

=item ERROR: Could not create root name dir

The program could not create the directory needed to hold the output for
the sequence file. It is possible that you not have the proper access to
create a directory in the path that you specified. It is also possible
that the directory that you are trying to create this subdirectory in 
does not exist.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

This program can optionally make use of a configuration file. The path to
configuration file is declared with the --config option at the command line.
This config file  will allow the 
batch_ltrseq.pl program to run the LTR_seq program for multiple parameter 
combinations for every fasta file in the input sequence directory. 
If no configuartion file is used, the LTR_seq program will be run with
default parameters.

The configuration file will
be a two column, tab delimited text file with the following options:

=over

=item * Col. 1 Configuration Set Name

This is the name that will be used in the gff file to tag the configuration
set. For example def for default options or old for an optional set to
try to finder LTR retrotransposon insertions that are old.

=item * Col. 2 Config File Path

The config file describing the parameter set to use. This file
should follow the format used for LTR_Seq config files. If this value is
set to DEF, then batch_ltrseq.p will run the LTR_seq program with
default parameters.

=back

An example configuration file is shown below:

    def     DEF
    long    long_ltr.cfg
    old     old_ltr.cfg

This configuration file would run three sets of LTR_seq configurations
for every sequence in the input directory. These configurations would
include (1) the default parameter set, (2) the parameter set name long
that is specified by the LTR_seq config file long_ltr.cfg, and (3) the
parameter set name old that is specified by the parameter set old_ltr.cfg.

Please see the LTR_seq documentation for specific information on the 
LTR_seq configuration files. The configuration files in LTR_seq allow you 
to specifiy the following parameters:

=over 2

=item * Dmin

The starting positions of the 5' and 3' LTRs should be separated by a minimum 
of this distance in base pairs.

=item * Dmax

The starting positions of the 5' and 3' LTRs should be separated by a maximum 
of this distance in base pairs.

=item * LTRmax

The 5' and 3' LTRs can individually span a maximum of this distance in base
pairs.

=item * LTRmin

The 5' and 3' LTRs can individually span a minimum of this distance in 
base pairs.

=item * LTRminExactMatch

The 5' and 3' LTRs should contain an exact match of this length in
base pairs.

=item * match

Dynamic programming score for a match (alignment).

=item * mismatch

Dynamic programming score for a mismatch (alignment).

=item * hgap

Dynamic programming score for a gap opening (alignment).

=item * gap

Dynamic programming score for a gap continuation (alignment).

=item * AlignmentWithN 

Dynamic programming score for aligning an `N' with any other base (alignment).

=item * MaxScoreRatioThreshold

Dynamic programming score threshold in percentage. This allows upto 15% 
difference in the scores of the optimal alignment vs. the best possible 
alignment given 100% identity.

=item * TSR_len

Length of the target site repeat.

=item * window

The LTR_seq documentation specifies that this parameter should not be
changed and is for internal use only.

=back

An example LTR_seq config file is shown below:

    window  12
    Dmin 600
    Dmax 15000
    LTRmax 2000
    match 2
    mismatch -5
    gap -1
    hgap -6
    AlignmentWithN -5
    MaxScoreRatioThreshold  15
    LTRmin 100
    TSR_len 6
    LTRminExactMatch 30

=head2 Environment

The batch_ltrseq.pl program can make use of the following variables
defined in the user environment.

=over 2

=item LTR_SEQ_DIR

This is the path to the directory that contains the LTR_seq configuration
files. This can also be specified at the command line using the 
--config-dir option. If this is not specified, the program will attempt
to look for the configuration files in the user's home directory.

=item LTR_SEQ_BIN

This is the path to the LTR_Seq binary. If this option is not specified, 
the batch_ltrseq.pl program will assume that LTR_seq is located in your
PATH.

=back

=head1 DEPENDENCIES

=head2 Required Software

The following software is required for the batch_ltrseq.pl program
to run properly.

=over

=item * ltr_seq

The ltr_seq is the latest incarnation of the ltr_par program
described by Kalyanaraman and Aluru J. Bioinform Comput Biol 4: 197-216.
A binary of this program is available upon contacuting the program
author: ananth <at> eecs.wsu.edu. The author's current homepage is:
http://www.eecs.wsu.edu/~ananth/

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

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
being used with batch_ltrseq.

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

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/06/2007

UPDATED: 03/24/2009

VERSION: $Rev: 728 $

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/06/2007
# - Program started, bulk of program written
# 
# 12/13/2007
# - Moved POD documentation to the end of the program
# - Added new print_help subfunction that extracts usage
#   and help statements from POD documentation
# - Updated POD documentation
#
# 01/26/2009
# - Updated POD documentation
# - Adding the ability to accept a configuration file
# - Added the ltrseq2gff subfunction
# - Added the ability to pass a configuration set name to the ltrseq2gff subfun
#
# 01/27/2009
# - Updated GFF output to use LTR feature names instead of exon
# - Added location of Target Site Duplications of GFF outfile
