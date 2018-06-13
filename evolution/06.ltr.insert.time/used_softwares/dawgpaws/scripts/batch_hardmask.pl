#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_hardmask.pl - Hardmask a batch of softmasked files  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/19/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a directory of softmasked fasta files, this will   |
#  hardmask the files by replacing the lowecase letters     |
#  with an uppercase letter set by the user.                |
#                                                           |
# USAGE:                                                    |
#  batch_hardmask.pl -i InDir -o OutDir -m X                |
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
# VARIABLE SCOPE              |
#-----------------------------+

# VARS WITH DEFAULT VALUES
my $out_ext = ".hard.fasta";   # Default outfile extension
my $mask_char = "N";           # Default character to mask with

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode

# PACKAGE LEVEL SCOPE
my $file_to_mask;              # Path to the file to mask
my $hard_mask_out;             # Path to the hardmasked file output
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file
my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $name_root;                 # Root name to be used for output etc


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # Optional strings
		    "m|mask=s"     => \$mask_char,
		    "ext=s"        => \$out_ext,
		    "logfile=s"    => \$logfile,
		    # Booleans
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);


my $bac_parent_dir = $outdir;  

my ( $ind_lib , $RepMaskCmd, $MakeGffDbCmd, $MakeGffElCmd );
my ( @RepLibs );
my $ProcNum = 0;

#//////////////////////
my $file_num_max = 5;
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
    print "\nbatch_hardmask.pl:\n".
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
    print LOG "  batch_hardmask.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
}


#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# Do not do the check if the indir is actually a single file
unless (-f $indir) {
    unless ($indir =~ /\/$/ ) {
	$indir = $indir."/";
    }
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
my @fasta_files;
if (-f $indir) {
    push (@fasta_files, $indir);
}
else {
    opendir( DIR, $indir ) || 
	die "Can't open directory:\n$indir"; 
    @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
    closedir( DIR );
}

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

print "Creating output dir ...\n" unless $quiet;
unless (-e $outdir) {
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
my $NumRecs = 1;
for my $ind_file (@fasta_files)
{
    
    $file_num++;

    #-----------------------------+
    # Get the root name of the    |
    # file to mask                |
    #-----------------------------+
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

    if (-f $indir) {
	$file_to_mask = $indir;
    }
    else {
	$file_to_mask = $indir.$ind_file;
    }

    $hard_mask_out = $outdir.$name_root.$out_ext;
    
    print "Converting:\n".
	"\t$file_to_mask TO\n".
	"\t$hard_mask_out\n" unless $quiet;


    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+
    open (IN, "<".$file_to_mask) ||
	die "Can not open input file:\n$file_to_mask\n";
    
    open (OUT, ">".$hard_mask_out) ||
	die "Can not open output file:\n$hard_mask_out\n";
    
    #-----------------------------+
    # HARD MASK FILE              |
    #-----------------------------+
    # The tr regexp does not appear to accept variables
    # therefore I have to write this a bit convoluted with
    # if then statements for acceptable MaskCharacters
    while (<IN>)
    {
	unless (m/^>/)        # Do not mask header lines 
	{
	    # Mask with the selected character
	    if ($mask_char =~ "N"){
		tr/[a-z]/N/;
	    } elsif ($mask_char =~ "X"){
		tr/[a-z]/X/;
	    } elsif ($mask_char =~ "x"){
		tr/[a-z]/x/;
	    } elsif ($mask_char =~ "n"){
		tr/[a-z]/n/;
	    } else {
		$msg = "\aERROR: A valid mask character was not selected\n";
	    }# End of select mask character
		
		# Print masked string to the outfile
		print OUT $_;
	    
	} else {
	    print OUT $_;
	    $NumRecs++;       # For headers increment NumRecs
	}
    } # End of while IN
    
    #-----------------------------+
    # CLOSE FILES                 |
    #-----------------------------+
    close IN;
    close OUT;
    
#//////////////////////////////////
# MAY WANT TO LEAVE THE FOLLOWING
# TO MAKE A COPY IN THE GENERAL
# BAC DIRECTORY
#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
#
##    #-----------------------------+
#    # MAKE OUTPUT DIR             |
##    #-----------------------------+
#    # the $bac_out_dir is the summary directory that
#    # contains all of the data related to a seq
#    $bac_out_dir = $outdir.$name_root."/";
#    mkdir $bac_out_dir, 0777 unless (-e $bac_out_dir); 


#    # TEMP EXIT FOR DEBUG, WIll JUST RUN FIRST FILE TO BE MASKED
#    if ($file_num > $file_num_max ) {
#	print "\nDebug run finished\n\n";
#	exit;
#    }


} # End of for each file in the input folder

close LOG if $logfile;

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


1;
__END__

=head1 NAME

batch_hardmask.pl - Hardmask a directory of softmasked fasta files. 

=head1 VERSION

This documentation refers to batch_hardmask version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    batch_hardmask.pl -i DirToProcess -o OutDir

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Given a directory of FASTA formatted sequence files that have been
softmasked, this will replace lowercase masked characters (a,g,c,t) 
with a masked character. You can choose from a limited set of
characters to mask with (x,X,n,N).

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --m,mask

Single letter to mask with. Valid options are: [ N | n | X | x ]

=item --ext

The new outfile extension to use. Default value is .hard.fasta

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

This program does not currently depend on any external configuration
files or and variables set in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * No external software required

The batch_hardmask.pl program is not dependent on any external
software.

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

=item * Limited fasta file extensins

Input files will only be recognized as fasta files if they have the
*.fa or *.fasta extension.

=back

=head1 SEE ALSO

The batch_hardmask.pl program is part of the DAWG-PAWS package of genome
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

STARTED: 07/19/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/19/2007
# - Program started
# - Imported functions from batch_mask.pl as well as
#   HardMask.pl
# 12/11/2007
# - Moved POD documentation to the end of the file
# - Added svn Rev 
# - Fixed strings that were set as booleans in option list
