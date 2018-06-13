#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_trf.pl - Run Tandem Repeat Finder in batch mode     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/10/2008                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  Run the Tandem Repeat Finder program in batch mode.      |
#                                                           |
# VERSION: Release 1.0                                      |
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
use Cwd;                       # Get the current working directory
use File::Copy;                # Copy files

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;

# Allowing TRF bin to be defined in environment
my $trf_bin = $ENV{TRF_BIN} || "trf400.linux.exe";  # The TRF binary file

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    # ADDITIONAL OPTIONS
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$do_test,
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

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\nbatch_trf.pl:\n".
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

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
my $start_dir = getcwd();

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
	    die "Could not create dir:\n$name_root_dir\n";
    }
    
    #-----------------------------+
    # CREATE THE TRF DIR          |
    #-----------------------------+
    my $trf_dir = $name_root_dir."trf/";
    unless (-e $trf_dir) {
	mkdir $trf_dir ||
	    die "Could not create trf out dir:\n$trf_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # This will hold the gff file modified from the
    # original gff fine
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create genscan out dir:\n$gff_dir\n";
    }


    #-----------------------------+
    # RUN THE TRF PROGRAM         |
    #-----------------------------+
    # TRF will place the results in the dir from which the program was called

    # Change to the fasta folder as the working directory
    # This will place the trf dat file in the fasta dir
    chdir $indir;
    
    # The parameter set used on 
    #trf400.linux.exe HEX0014K09.masked.fasta 2 2 7 80 10 30 5
    my $trf_match = "2";        # Matching weight
    my $trf_mismatch = "2";     # Mismatch penalty
    my $trf_delta = "7";        # Indel probability
    my $trf_pm = "80";          # Match probability
    my $trf_pi = "10";          # Indel probability
    my $trf_min = "30";         # TRF minimum alignment score
    my $trf_max = "5";          # Maximum period size

    # h switch to supress html output
    # d switch to produce a *.dat file
    my $trf_cmd = "$trf_bin".
	" $ind_file".
	" $trf_match".
	" $trf_mismatch".
	" $trf_delta".
	" $trf_pm".
	" $trf_pi".
	" $trf_min".
	" $trf_max".
	" -h -d";                   

    print "TRF CMD: $trf_cmd\n" if $verbose;


    system ($trf_cmd);
    

    # change back to the starting working dir
    # This will make sure that relative paths still work

    #-----------------------------+
    # MOVE TRF RESULT FROM FASTA  |
    # DIR TO THE TRF DIR          |
    #-----------------------------+
    chdir $start_dir;

    # The base file name for the trf result
    my $trf_result_file = "$ind_file.$trf_match.$trf_mismatch.$trf_delta.".
	"$trf_pm.$trf_pi.$trf_min.$trf_max.dat";
    my $trf_result_path = "$indir$trf_result_file";  # Original TRF placement
    my $dp_result_path = "$trf_dir$trf_result_file"; # DawgPaws placement

    move( $trf_result_path, $dp_result_path ) ||
	die "Can not move:\n $trf_result_path TO\n$dp_result_path\n";
    

    #-----------------------------+
    # TRANSLATE TRF TO GFF        |
    #-----------------------------+
    open (TRFIN, "<$dp_result_path") ||
	die "Can not open TRF dat file:\n$dp_result_path\n";

    my $gff_out = $gff_dir.$name_root."_trf.gff";
    open (TRFGFF, ">$gff_out") ||
	die "Can not open GFF output file $gff_out";
    
    my $trf_count = 0; # Trf feature count
    while (<TRFIN>) {

	# Print the trf result
	#print STDERR"$_\n" if $verbose;
	
	my @trf_parts = split;
	my $num_trf_parts = @trf_parts;
	
	if ($num_trf_parts == 15) {

	    $trf_count++;  # Increment feature count
	    my $pad_len = 4;
	    my $trf_count_pad = sprintf("%0${pad_len}d", $trf_count);


	    # Format count for id name
	    my $start = $trf_parts[0];
	    my $end = $trf_parts[1];
	    my $repeat_word = $trf_parts[13];

	    #print STDERR"\tSTART: $start\n\tEND: $end\n" if $verbose;
	    #print STDERR "\tWORD: $repeat_word\n" if $verbose;

	    # PRINT THE GFF OUTPUT

	    # The gff string
	    my $gff_str = "$name_root\t".  # 1 Sequence name
		"TRFv4.00\t".             # 2 Source program
		"tandem_repeat\t".        # 3 Feature type
		"$start\t".               # 4 Feature Start
		"$end\t".                 # 5 Feature End
		".\t".                    # 6 Score
		".\t".                    # 7 Strand, could be +
		".\t".                    # 8 Frame
		"TRF$trf_count_pad\n";    # 9 Feature ID
	    print TRFGFF $gff_str;

	    print STDERR $gff_str if $verbose;

	} # END of if correct number of  TRF parts
    
    } # End of while TRFIN

    close (TRFIN);
    close (TRFGFF);


}

# The following required to make the test harness work
print STDOUT "\n";

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

1;
__END__

=head1 NAME

batch_trf.pl - Run Tandem Repeat Finder in batch mode

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    batch_trf.pl -i InDir -o OutDir

=head2 Required Arguments

    --indir         # Path to the input directory
    --outdir        # Path to the output directory

=head1 DESCRIPTION

This program runs the Tandem Repeat Finder program (TRF) in batch mode.
Given a directory of fasta files, this will run the TRF program for all
files in the directory.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -c, --config

Path to a config file. This is a tab delimited text file
indicating the required information for each of the databases to blast
against. Lines beginning with # are ignored.

=back

=head1 OPTIONS

=over 2

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

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

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

This program does not make use of a configuration file.

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * Tandem Repeats Finder

The batch_trf.pl program requires the Tandem Repeats Finder program. This
program is available from:
http://tandem.bu.edu/trf/trf.download.html

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited testing on versions of TRF

This program has been tested and is known to work with trf400.linux.exe

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
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

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 11/10/2008

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 11/10/2008
# - Main body of program written
# 
# 01/20/2009
# - Updated POD documentation
# - Added Rev to properties
