#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_gff2game.pl - Convert a gff file to game xml          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 02/09/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts gff data tracks from gff format to the game     |
#  xml format for use in the Apollo Genome Annotation       |
#  Curation program.                                        |
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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLES
#-----------------------------+

# CONSTANTS
#my $VERSION = "1.0";

my $infile_fasta;
my $infile_gff;
my $outfile;

# Options with default values
my $ap_path = $ENV{DP_APOLLO_BIN} || "apollo";

# Booleans
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED
		    "i|infile=s"  => \$infile_fasta,
		    "g|gff=s"     => \$infile_gff,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "ap-path=s"   => \$ap_path,
		    # BOOLEANS
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
		    "q|quiet"     => \$quiet,);

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
    print "\ncnv_gff2game.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK FOR REQUIRED VARS     |
#-----------------------------+
if ( (!$infile_fasta) || (!$infile_gff) || (!$outfile) ) {
    print_help("full");
}

# Convert each the gff file to game xml using the apollo_convert command
&apollo_convert ($infile_gff, "gff", $outfile, 
		 "game", $infile_fasta, "NULL");

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


sub apollo_convert {
#-----------------------------+
# CONVERT AMONG FILE FORMATS  |
# USING THE APOLLO PROGRAM    |
#-----------------------------+
# Converts among the various data formats that can be used 
# from the command line in tbe Apollo program. For example
# can convert GFF format files into the game XML format.
# NOTES:
#  - Currently assumes that the input file is in the correct
#    coordinate system.
#  - GFF files will require a sequence file
#  - ChadoDB format will require a db password


    # ApPath - the path of dir with the Apollo binary
    #          Specifying the path will allow for cases
    #          where the program is not in the PATHS
    # ApCmd  - the apollo commands to run

    my ($InFile,$InForm,$OutFile,$OutForm,$SeqFile,$DbPass) = @_;

    #InFile = $_[0];        # Input file path
    #InForm = $_[1];        # Output file format:
    #                         # game|gff|gb|chadoxml|backup
    #$OutFile = $_[2];       # Output file path
    #$OutForm = $_[3];       # Ouput file foramt
    #                           # chadoDB|game|chadoxml|genbank|gff|backup
    #$SeqFile = $_[4];       # The path of the sequence file
    #                           # This is only required for GFF foramt files
    #                           # When not required this can be passed as na
    #$DbPass = $_[5];        # Database password for logging on to the 
    #                           # chado database for reading or writing.
    my $ApCmd;

    # Set the base command line. More may need to be added for different
    # formats. For example, GFF in requires a sequence file and the CHADO
    # format will require a database password.
    $ApCmd = $ap_path." -i ".$InForm." -f ".$InFile.
	" -o ".$OutForm." -w ".$OutFile;

    # Make sure that that input output formats are in lowercase
    # may need to add something here to avoid converting chadoDB
    $InForm = lc($InForm);
    $OutForm = lc($OutForm);
    
    # Determine the proper command to use based on the input format
    # since GFF file also require a sequence file
    if ($InForm =~ "gff" ) {
	$ApCmd = $ApCmd." -s ".$SeqFile;
    }
    
    if ($InForm =~ "chadodb") {
	$ApCmd = $ApCmd." -D ".$DbPass;
    }

    # Do the apollo command
    system ( $ApCmd );

}

1;
__END__

=head1 NAME

cnv_gff2game.pl - Convert a gff file to game xml

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_gff2game.pl -i InFile.fasta -g GffFile.gff -o OutFile.game.xml

=head2 Required Arguments

    -i              # Path to the fasta file the gff refers to
    -g              # Path to the gff file to convert
    -o              # Path to the output game.xml file

=head1 DESCRIPTION

Converts gff foramt files to the game.xml format. This program uses the
Apollo genome annotation program to do this conversion, so you must have
Apollo installed for this script to work

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input fasta filefile.

=item -g, --gff

Path of the input gff file

=item -o,--outfile

Path of the output game xml file.
 
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

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

=over

=item Can not open input file

The input file path provided by the -i, --infile switch is not valid.

=item Can not open output file

The output file path provided by the -o, --outfile witch is not valid.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head 2 Configuration

This program does not make use a configuration file.

=head 2 Environment

The following variables can be defined in the user environment:

=over

=item * DP_APOLLO_BIN

The location of the apollo binary file. This is the path used to
lauch the Apollo genome annotation program. If not specified in the
user environment, this will attempt to call 'apollo'. The path may
also be specified using the --ap-path option at the command line.

=back

=head1 DEPENDENCIES

=head2 Required Software

=over 2 

=item Apollo

This program requires that apollo be installed on the local machine
since apollo is being used as the engine to do the
conversion between formats. Apollo is available from:
http://apollo.berkeleybop.org/current/index.html

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

=item * Limited Apollo Testing

New versions of Apollo come out on a regular basis, and older versions
of Apollo are not arcived. I am therefore limited to testing this program
to versions of Apollo that are extant at the time this script was written.
This is known to work on Apollo version 1.6.5.

The command line interface from Apollo does not appear to work on
text only terminals.

=item * Does not support STDIN/STDOUT

In general the cnv scripts are designed to support input from STDIN and
output from STDOUT. The gff to game converter relies on Apollo, and currently
does not support the STDIN and STOUT. :(

=back

=head1 SEE ALSO

This program is part of the DAWG-PAWS package of genome
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

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 02/09/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 01/19/2009
# - Cleanup of code
# - Adding the ability to read from STDIN and write to 
#   STDOUT
# - Added the new print_help subfunction
# - Added package namespace
