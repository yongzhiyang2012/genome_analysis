#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_findltr.pl - Run the findltr program in batch mode  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/13/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the find_ltr.pl LTR finding program in batch mode.   |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
# VERSION: Release 1.0                                      |
#                                                           |
#-----------------------------------------------------------+

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
use File::Copy;
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
my $indir;                     # Base intput dir
my $outdir;                    # Base output dir
my $config_file;               # Configuration file
my $name_root;                 # Root name of the file being analyzed
my $errmsg;                    # Var to hold error message strings
my $fl_gff_outpath;            # Full output path of the gff output file

# Array of find_ltr parameters
my @fl_params = ();            # 2d Array to hold the find_ltr parameters 

# Path vars
my $fl_path = $ENV{FIND_LTR_PATH};  # Path to the find_ltr.pl program
my $gff_dir;                   # Dir to hold the gff output

# Counters/Index Vals
my $i = 0;                     # Array index val
my $file_num = 0;              # File number
my $proc_num = 0;              # Process number

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff_convert = 1;
my $do_test = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    "c|config=s"  => \$config_file,
		    # ADDITIONAL OPTIONS
		    "fl-path=s"   => \$fl_path,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "gff"         => \$do_gff_convert,
		    "test"        => \$do_test,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

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
    print "\nbatch_findltr.pl:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

# Throw error if required arguments not present
if ( (!$indir) || (!$outdir) || (!$config_file) ) {
    print "\a";
    print STDERR "ERROR: An input directory must be specified" if !$indir;
    print STDERR "ERROR: An output directory must be specified" if !$outdir;
    print STDERR "ERROR: A config file must be specified" if !$config_file;
    print_help ("help",  $0 );
    exit;
}

# If not path for find_ltr.pl exists as an ENV option, and the
# the path is not given at the command line, assume that find_ltr.pl
# exists somewhere that the user has included in their PATH var.
if ( !$fl_path ) {
    $fl_path = "find_ltr.pl";
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
       	my @in_line = split;           # Implicit split of $_ by whitespace
	my $num_in_line = @in_line; 
	
	if ($num_in_line == 10) { 
	    $fl_params[$i][0] = $in_line[0] || "NULL";  # Name
	    $fl_params[$i][1] = $in_line[1] || "NULL";  # MIN-MEM
	    $fl_params[$i][2] = $in_line[2] || "NULL";  # MIN-MEM-DIST
	    $fl_params[$i][3] = $in_line[3] || "NULL";  # MAX-MEM-DIST
	    $fl_params[$i][4] = $in_line[4] || "NULL";  # MAX-MEM-GAP
	    $fl_params[$i][5] = $in_line[5] || "NULL";  # MIN-LEN-LTR 
	    $fl_params[$i][6] = $in_line[6] || "NULL";  # MAX-LEN-LTR
	    $fl_params[$i][7] = $in_line[7] || "NULL";  # RANGE-BIN 
	    $fl_params[$i][8] = $in_line[8] || "NULL";  # MIN LEN ORF
	    $fl_params[$i][9] = $in_line[9] || "NULL";  # MAX E Value HMM
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

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $num_files = @fasta_files;

my $num_proc_total = $num_files * $num_par_sets;

print STDERR "$num_proc_total find_ltr runs to process\n";

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


    # The following added for temp work with gff output
    # 09/28/2007
#    if ($file_num == 5) {exit;}
#    print "Processing: $name_root\n";

    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE FIND_LTR OUTDIR      |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $findltr_dir = $name_root_dir."find_ltr/";
    unless (-e $findltr_dir) {
	mkdir $findltr_dir ||
	    die "Could not create genscan out dir:\n$findltr_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    if ($do_gff_convert) {
	$gff_dir = $name_root_dir."gff/";
	unless (-e $gff_dir) {
	    mkdir $gff_dir ||
		die "Could not create genscan out dir:\n$gff_dir\n";
	}
    }

    #-----------------------------------------------------------+
    # RUN find_ltr.pl FOR EACH SET OF PARAM VALS IN CONFIG FILE |
    #-----------------------------------------------------------+
    for ($i=0; $i<$num_par_sets; $i++) {
	
	# ASSUME THAT A CONFIG FILE IS USED
	# OTHERWISE COULD JUST USE DEFAULT VALS
	$proc_num++;
	
	# Load array vals to usefull short names
	my $fl_suffix       = $fl_params[$i][0]; # Name of the parameter set
	my $fl_min_mem      = $fl_params[$i][1];
	my $fl_min_mem_dist = $fl_params[$i][2];
	my $fl_max_mem_dist = $fl_params[$i][3];
	my $fl_max_mem_gap  = $fl_params[$i][4];
	my $fl_min_ltr      = $fl_params[$i][5];
	my $fl_max_ltr      = $fl_params[$i][6];
	my $fl_range_bin    = $fl_params[$i][7];
	my $fl_min_orf      = $fl_params[$i][8];
	my $fl_e_val        = $fl_params[$i][9];

	# Show processing information
	print "Processing: $name_root $fl_suffix\n";
	
	# The following name will need to be copied from the default 
	# named output from find_ltr
	my $fl_inseqpath = $indir.$ind_file;
	my $fl_ltrpos_out = $fl_inseqpath.".ltrpos";
	my $fl_out_copy = $findltr_dir.$name_root."_findltr_".$fl_suffix.".txt";
	my $fl_ltrseq_out = $fl_inseqpath.".ltrseq";
	my $fl_ltrseq_out_cp = $findltr_dir.$name_root."_findltr_".
	    $fl_suffix.".ltrseq";
	if ($do_gff_convert) {
	    $fl_gff_outpath = $gff_dir.$name_root."_findltr_".$fl_suffix.".gff";
	}
	
	#-----------------------------+
	# RUN THE find_ltr.pl program |
	#-----------------------------+
	my $fl_cmd = $fl_path. 
	    " --seq ".$fl_inseqpath.
	    " --min-mem ".$fl_min_mem.
	    " --min-dist ".$fl_min_mem_dist.
	    " --max-dist ".$fl_max_mem_dist.
	    " --max-gap ".$fl_max_mem_gap.
	    " --min-ltr ".$fl_min_ltr.
	    " --max-ltr ".$fl_max_ltr.
	    " --min-ltr ".$fl_min_ltr.
	    " --range-bin ".$fl_range_bin.
	    " --min-orf ".$fl_min_orf.
	    " --e-val ".$fl_e_val;

	if ($verbose) {
	    $fl_cmd = $fl_cmd." --verbose";
	}

	print "\n---------------------------------------+\n" if $verbose;
	print " Process $proc_num of $num_proc_total\n" if $verbose;
	print "---------------------------------------+\n" if $verbose;

	print STDERR "FL CMD: $fl_cmd\n\n" if $verbose;

	system ($fl_cmd) unless $do_test;

	#-----------------------------+
	# CONVERT OUTPUT TO GFF       |
	#-----------------------------+
	if ( (-e $fl_ltrpos_out) && ($do_gff_convert)) {
	    findltr2gff ( $fl_ltrpos_out, $fl_gff_outpath, 
			  0, $name_root, $fl_suffix);
	}


	#-----------------------------+
	# MOVE RESULTS FILES          |
	#-----------------------------+
	if (-e $fl_ltrpos_out) {
	    $errmsg = "Could not move\n $fl_ltrpos_out to\n $fl_out_copy\n";
	    move ($fl_ltrpos_out, $fl_out_copy)
		|| print STDERR "\a$errmsg";
	}

	if ($fl_ltrseq_out) {
	    $errmsg = "Could not move\n $fl_ltrseq_out to \n  $fl_ltrseq_out_cp\n";
	    move ($fl_ltrseq_out,  $fl_ltrseq_out_cp )
		|| print STDERR "\a$errmsg";
	}


    }

} # End of for each ind_file in @fasta_files

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


sub findltr2gff {

    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    # gff_suffix is the name appended to the end of the gff_source
    my ($findltr_in, $gff_out, $append_gff, $seqname, $gff_suffix) = @_;

    # find_ltr
    my $gff_source;                 # 
    my $findltr_id;                 # Id as assigned from find_ltr.pl
    my $findltr_name;               # Full name for the find_ltr prediction
    my $ltr5_start;                 # Start of the 5' LTR
    my $ltr5_end;                   # End of the 5' LTR
    my $ltr5_len;                   # Length of the 5' LTR
    my $ltr3_start;                 # Start of the 3' LTR
    my $ltr3_end;                   # End of the 3' LTR
    my $ltr3_len;                   # Length of the 3' LTR
    my $el_len;                     # Length of the entire element
    my $mid_start;                  # Start of the LTR Mid region
    my $mid_end;                    # End of the LTR Mid region
    my $ltr_similarity;             # Percent similarity between LTRs
    my $ltr_strand;                 # Strand of the LTR

    my @in_split = ();              # Split of the infile line
    my $num_in;                     # Number of split vars in the infile

     # Initialize Counters
    my $findltr_num = 0;            # ID Number of putatitve LTR retro

    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+
    open (INFILE, "<$findltr_in") ||
	die "Can not open input file:\n$findltr_in\n";

    if ($append_gff) {
	open (GFFOUT, ">>$gff_out") ||
	    die "Could not open output file for appending\n$gff_out\n";
    }
    else {
	open (GFFOUT, ">$gff_out") ||
	    die "Could not open output file for output\n$gff_out\n";
    } # End of if append_gff
    
    #-----------------------------+
    # PROCESS INFILE              |
    #-----------------------------+
    while (<INFILE>) {
	chomp;

	my @in_split = split;
	my $num_in = @in_split;   
	
	# Load split data to vars if expected number of columns found
	if ($num_in == 10) {

	    $findltr_num++;

	    $findltr_id = $in_split[0];
	    $ltr5_start = $in_split[1];
	    $ltr5_end = $in_split[2];
	    $ltr3_start = $in_split[3];
	    $ltr3_end = $in_split[4];
	    $ltr_strand = $in_split[5];
	    $ltr5_len = $in_split[6];	    
	    $ltr3_len = $in_split[7];
	    $el_len = $in_split[8];
	    $ltr_similarity = $in_split[9];

	    $mid_start = $ltr5_end + 1;
	    $mid_end = $ltr3_start - 1;   

	    $findltr_name = $seqname."_findltr_"."".$findltr_id;
	    $gff_source = "find_ltr:".$gff_suffix;

	    #-----------------------------+
	    # PRINT GFF OUTPUT            |
	    #-----------------------------+
	    # Data type follows SONG
	    # http://song.cvs.sourceforge.net/*checkout*/song/ontology/so.obo

	    # Full span of LTR Retrotransposon
	    print GFFOUT "$seqname\t".  # Name of sequence
		"$gff_source\t".        # Source
		"LTR_retrotransposon\t".# Features, exon for Apollo
		"$ltr5_start\t".        # Feature start
		"$ltr3_end\t".	        # Feature end
		".\t".                  # Score, Could use $ltr_similarity
		"$ltr_strand\t".        # Strand
		".\t".                  # Frame
		"$findltr_name\n";      # Features (name)

	    # 5'LTR
	    print GFFOUT "$seqname\t". # Name of sequence
		"$gff_source\t".       # Source
		"five_prime_LTR\t".    # Features, exon for Apollo
		"$ltr5_start\t".       # Feature start
		"$ltr5_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".       # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)

#	    # MID
#	    # This is not a SONG complient feature type name
#	    print GFFOUT "$seqname\t". # Name of sequence
#		"$gff_source\t".       # Source
#		"mid\t".               # Features, exon for Apollo
#		"$mid_start\t".        # Feature start
#		"$mid_end\t".	       # Feature end
#		".\t".                 # Score, Could use $ltr_similarity
#		"$ltr_strand\t".       # Strand
#		".\t".                 # Frame
#		"$findltr_name\n";     # Features (name)
	    
	    # 3'LTR
	    print GFFOUT "$seqname\t". # Name of sequence
		"$gff_source\t".       # Source
		"three_prime_LTR\t".   # Features, exon for Apollo
		"$ltr3_start\t".       # Feature start
		"$ltr3_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".       # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)

	} # End of if num_in is 10

    } # End of while INFILE


} # End of findltr2gff

1;
__END__

=head1 NAME

batch_findltr.pl - Run the find_ltr.pl program in batch mode.

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    batch_findltr.pl -i InDir -o OutDir -c Config.cfg [--gff]

=head2 Required Arguments

    --indir         # Path to the input directory of fasta files
    --outdir        # Path to the base output directory
    --config        # Config file containg batch_findltr.pl paramaters
    --gff           # Produce GFF formatted output

=head1 DESCRIPTION

Runs the program find_ltr.pl in batch mode. This makes use of a modified
version of the find_ltr.pl program that takes changes to the LTR finding
parameters at the command line.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the input directory. This is the directory that contains all
of the fasta files to anlayze. The fasta files should all end with
the I<fasta> extension to recognized. 

=item -o,--outdir

Path of the output directory. This is the base directory that will 
hold all of the batch_findltr.pl output

=item -c, --config 

Path of the config file that contains the model options for running
find_ltr. This config file is a white space delimited text file that 
should be in the following format.

  #---------------------------------------------------------------+
  #1   2     3     4      5     6     7      8     9   10         |
  #---------------------------------------------------------------+
  Def  40  1100  16000    40   100   1000   500  700   0.0000000001
  Alt  40  1100  1800     40   100   1000   500  400   0.00001 

More information about this file is available under configuration and 
environment heading below.

=back

=head1 OPTIONS

=over 2

=item --fl-path

Location of the find_ltr.pl program. This option can also be set in the
users envrionment. See Configuration and Environment below.

=item -q,--quiet

Run the program with minimal output.

=item -v, --verbose

Run the program in verbose mode.

=item --gff

Produce gff formatted output of the results.

=item --test

Run the program in test mode. The find_ltr.pl program will not be run, but
the location of source files, binaries, will be checked and the 
outupt directories will be created.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

The configuration file in batch_findltr.pl specifies the options for
running the find_ltr.pl program. This is a white space delimited text
file. All lines starting with the # symbol will be treated as comments.
An example of a config file is below:

  #---------------------------------------------------------------+
  #1   2     3     4      5     6     7      8     9   10         |
  #---------------------------------------------------------------+
  Def  40  1100  16000    40   100   1000   500  700   0.0000000001
  Alt  40  1100  1800     40   100   1000   500  400   0.00001 

These 10 columns represents the following information:

=over 2

=item Col 1. 

Base_name for the parameter set. This set name will be used to
name the output file, and will be added to the output of
the gff output file. B<DO NOT INCLUDE SPACES IN NAMES>

=item Col 2. 

Minimum Length MEM

=item Col 3. 

Mimimum distance between MEMs

=item Col 4. 

Maximum distance between MEMs

=item Col 5. 

Maximu gap between MEMs

=item Col 6. 

Minimum length of the LTR 

=item Col 7. 

Maximum length of the LTR

=item Col 8. 

Range Bin

=item Col 9. 

Minimum length of ORF

=item Col 10. 

Mac E value of HMM Hit

=back

=head2 FIND_LTR_PATH Environment

As an alternative to specifying the full path of the find_ltr program
with the --fl-path option,
the path of the find_ltr program can be specified in the users environment.
For example in bash shell, add the following line to your .bashrc

  export FIND_LTR_PATH='/usr/local/genome/find_ltr.pl'

assuming that the find_ltr.pl program is in the /usr/local/genome/
directory.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * find_ltr.pl>

A modified version of the find_ltr.pl program is required.

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the BLAST results.

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

=item * Modified version of findltr required

The batch_finltr.pl program requies a modified version of the batch_ltr.pl
program that accepts parameters from the command line.

=item * Config file must use UNIX format line endings

The config file must have UNIX formatted line endings. Because of
this any config files that have been edited in programs such as
MS Word must be converted to a UNIX compatible text format before
being used with batch_blast.

=back

=head1 SEE ALSO

The batch_findltr.pl program is part of the DAWG-PAWS package of genome
annotation programs. See the DAWG-PAWS web page 
( http://dawgpaws.sourceforge.net/ )
or the Sourceforge project page 
( http://sourceforge.net/projects/dawgpaws ) 
for additional information about this package.

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

STARTED: 09/13/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 09/13/2007
# - Basic input and output working
# 09/14/2007
# - Added findltr2gff from cnv_findltr2gff.pl 
# 09/28/2007
# - Making gff output type column SONG complient
# 12/10/2007
# - Updated POD documentation
# - Converted to new help subfunction that extracts help
#   and usage message from the POD documentation
