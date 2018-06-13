#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_genscan.pl - Run genscan gene prediction program    |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill at gmail.com                         |
# STARTED: 07/31/2007                                       |
# UPDATED: 04/28/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the genscan gene prediction program in batch mode.   |
#  Runs genscan as well as converts output to gff format.   |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use File::Copy;
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
# VARS USING ENV              |
#-----------------------------+
# If not provided here, use the default assumtion that
# the file location are in the user's path and will
# use Maize as the default library.
my $genscan_path = $ENV{DP_GENSCAN_BIN} ||
    "genscan";
my $lib_path = $ENV{DP_GENSCAN_LIB} ||
    "Maize.smat";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file

my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $name_root;                 # Root name to be used for output etc

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;
my $verbose = 0;

# COUNTERS
my $num_proc = 1;              # Number of processors to use

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # Optional strings
		    "genscan-path" => \$genscan_path,
		    "lib-path"     => \$lib_path,
		    "logfile=s"    => \$logfile,
		    # Booleans
		    "apollo"       => \$apollo,
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);

my $proc_num = 0;

#//////////////////////
my $file_num_max = 2;
my $file_num = 0;
#\\\\\\\\\\\\\\\\\\\\\\

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
    print "\nbatch_genscan.pl:\n".
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
# OPEN THE LOG FILE           |
#-----------------------------+
if ($logfile) {
    # Open file for appending
    open ( LOG, ">>$logfile" ) ||
	die "Can not open logfile:\n$logfile\n";
    my $time_now = time;
    print LOG "==================================\n";
    print LOG "  batch_mask.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
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
my $num_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($num_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# RUN REPEAT MAKSER AND PARSE |
# RESULTS FOR EACH SEQ IN THE |
# fasta_files ARRAY FOR EACH  |
# REPEAT LIBRARY IN THE       |
# RepLibs ARRAY               |
#-----------------------------+

for my $ind_file (@fasta_files) {
    
    $proc_num++;
    $file_num++;
    #if ($file_num == $file_num_max){exit;}

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FASTA FILE                  |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
	# This is hard masked fasta files
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.masked\.fasta$/) {
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

    my $infile_path = $indir.$ind_file;
    

    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root;
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE GENSCAN OUTDIR       |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $genscan_dir = $outdir.$name_root."/gene/";
    unless (-e $genscan_dir) {
	mkdir $genscan_dir ||
	    die "Could not create genscan out dir:\n$genscan_dir\n";
    }
    
    my $gff_dir = $outdir.$name_root."/gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create gff out dir:\n$gff_dir\n";
    }

    my $out_path = $genscan_dir.$name_root.".genscan.out";
    #my $gff_path = $genscan_dir.$name_root.".genscan.gff";
    my $gff_path = $gff_dir.$name_root.".genscan.gff";

    my $genscan_cmd = "$genscan_path $lib_path $infile_path -v > $out_path";

    print "=======================================\n" if $verbose;
    print "Running Genscan for $name_root\n" if $verbose;
    print " File $file_num of $num_files\n" if $verbose;
    print "=======================================\n" if $verbose;
    print "$genscan_cmd\n" if $verbose;
    system($genscan_cmd) unless $test;

    print "Converting Genscan output\n";
    if (-e $out_path) {
	genscan_2_gff($out_path, $gff_path, $name_root) unless $test;
    }
    else {
	print "ERROR: Could not find genscan output at:\n$out_path\n"
    }

} # End of for each file in the input folder

close LOG if $logfile;

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub genscan_2_gff 
{
# Convert a Genscan format output to the GFF format
# Modified from
# http://nucleus.cshl.org/agsa/presentations.html
# as Downloaded by JCE 02/22/2007

    my ($infile, $outfile, $name_root) = @_;

    my %exon_type = ('Sngl', 'Single Exon',
		     'Init', 'Initial Exon',
		     'Intr', 'Internal Exon',
		     'Term', 'Terminal Exon');
    
    open (IN, "<".$infile);
    open (OUT, ">".$outfile);
    
    while (<IN>) {
	
	# Last line before predictions contains nothing but spaces and dashes
	if (/^\s*-[-\s]+$/)  {
	    while (<IN>)  {
		my %feature; 
		if (/init|term|sngl|intr/i) {
		    
		    my @f  = split;
		    
		    my ($gene, $exon) = split (/\./, $f[0]); 
		    
                    #name must be a number
		    $feature {name} = $gene + ($exon/1000); 

		    #arrange numbers so that start is always < end
		    if ($f[2] eq '+') {
			$feature {'start'}  = $f[3];
			$feature {'end'}    = $f[4];
			$feature {'strand'} = "+";
		    } elsif ($f[2] eq '-') {
			$feature {'start'}  = $f[4];
			$feature {'end'}    = $f[3];
			$feature {'strand'} = "-";
		    }
		    
		    $feature {'score'}   = $f[12];
		    $feature {'p'}       = $f[11];
		    $feature {'type'}    = $exon_type{$f[1]};
		    $feature {'program'} = 'Genscan';
		    $feature {'program_version'} = '1.0';
		    $feature {'primary'} = 'prediction';
		    $feature {'source'}  = 'genscan';
		    
		    # Pring GFF format output
		    print OUT 
			"$name_root\t" .            # SeqName
			"Genscan:maize\t".          # Source
			"exon\t".                     # Feature
			$feature{start}."\t" .      # Start
			$feature{end}."\t" .        # End
			$feature{p}."\t" .          # Score
			$feature{strand} . "\t".    # Strand
			".\t" .                     # Frame
			"gene_".$gene."\n";               # Attribute
		} elsif (/predicted peptide/i) {
		    last;   
		}
	    } # End of second while statement
	} # End of if seach command
    } # End of first while statment

} #End of genscan_2_gff subfunction


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


# The deprecated print_help subfunction

sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;
    
    my $usage = "USAGE:\n".
	"  batch_genscan.pl -i DirToProcess -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directory containing the sequences\n".
	"                 # to process. The files must have one of the\n".
	"                 # following file extensions:\n".
	"                 # [fasta|fa]\n".
	"  --outdir       # Path to the output directory\n".
	"\n".
	"OPTIONS:\n".
	"  --gencan-path  # Full path to the genscan binary\n".
	"  --lib-path     # Full path to the prediction library:\n".
	"  --logfile      # Path to file to use for logfile\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
	"  --test         # Run the program in test mode\n".
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

batch_genscan.pl - Run genscan and parse results to a gff format file. 

=head1 VERSION

This documentation refers to batch_genscan.pl version $Rev: 782 $

=head1 SYNOPSIS

=head2 Usage

    batch_genscan.pl -i DirToProcess -o OutDir

=head2 Required Variables

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Run the GENSCAN gene prediction program in batch mode.
This will run genscan as well as convert the output to gff format.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --genscan-path

The full path to the genscan binary.

=item --lib-path

The full path to the library file.

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

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

=item --test

Run the program without doing the system commands.

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

The batch_genscan.pl program currently does not make use of any 
external configuration files.

=head2 User Environment

The following environment variables may be used to designate the full
path to the GENSCAN binary and the default library path to use.

=over 

=item DP_GENSCAN_BIN

The path to the genscan binary file.

=item DP_GENSCAN_LIB

The path to the default library to use for gene predictions using genscan.

=back

The following example shows the lines that would need to be added to your
.bashrc file in the bash shell to specify the location of the genscan
binary and the Maize library file in your home directory.

    export DP_GENSCAN_BIN='/home/yourname/apps/genscan/genscan'
    export DP_GENSCAN_LIB='/home/yourname/apps/genscan/Maize.smat'

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * GENSCAN

This program obviously requires the GENSCAN gene model prediction program.
Information for downloading the GENSCAN program is available from:
http://genes.mit.edu/GENSCANinfo.html

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the results.

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
being used with batch_blast.

=back

=head1 SEE ALSO

The batch_genscan.pl program is part of the DAWG-PAWS package of genome
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

STARTED: 07/31/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/31/2007
# - Program started and bulk of functions added and
#   working copy of the program finished
# 
# 12/11/2007
# - Added SVN tracking of Rev
# - Moved POD documentation to the end of the file
# - Added use strict
# - Added print_help subfunction that extracts help
#   and usage statements from the POD documentation
# - Update POD documentation
# - Changed the path to the genscan binary and the default
#   libary to variables set in the user environment
#
# - Changed Version reference to $0 instead of hard coded
#   program name
# - Additional updates to POD documentation
#
# 04/28/2009
# - Putting output in gff dir
