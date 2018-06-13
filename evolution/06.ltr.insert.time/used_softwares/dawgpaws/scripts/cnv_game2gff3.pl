#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_game2gff2.pl - Convert game.xml to gff3 format        |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 02/17/2009                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert a game xml format file to the gff3 format. This  |
#  program absolutely requires a local installation of      |
#  Apollo since it uses Apollo to do the conversion.        |
#                                                           |
# USAGE:                                                    |
#  cnv_game2gff.pl -i HEX001.game.xml -o HEX001.gff         |
#                                                           |
# VERSION: Release 1.0                                      |
#                                                           |
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
my $infile;
my $outfile;
my $apollo_bin = $ENV{DP_APOLLO_BIN} || "apollo";

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $test=0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "ap-path=s"   => \$apollo_bin,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$test,
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
    print "\ncnv_game2gff3.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

# check if input is actually a directory
# If it is a file, convert the single file
if (-f $infile) {
    apollo_convert ($apollo_bin, $infile, "game", $outfile, "gff3");
}
elsif (-d $infile) {
    my $indir = $infile;
    # Read the xml files in the directory
    unless ($indir =~ /\/$/ ) {
	$indir = $indir."/";
    }

    opendir( DIR, $indir ) || 
	die "Can't open directory:\n$indir"; 
    my @xml_files = grep /\.xml$/, readdir DIR ;
    closedir( DIR );
    
    my $count_files = @xml_files;
    
    # Report error if no xml files in the directory
    if ($count_files == 0) {
	print STDERR "\a";
	print STDERR "\nERROR: No xml files were found in the input".
	    " directory\n $indir\n".
	    "Game xml files must have the xml extension to be recognized.\n\n";
	exit;
    }

    #-----------------------------+
    # CONVERT EACH FILE           |
    #-----------------------------+
    for my $ind_file (@xml_files) {

	apollo_convert ($apollo_bin, $indir.$ind_file, "game", 
			$indir.$ind_file.".gff", "gff3");
	
    }

}

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

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

    my ($ap_path, $InFile,$InForm,$OutFile,$OutForm,$SeqFile,$DbPass) = @_;

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

cnv_game2gff3.pl - Convert game.xml to gff3 format  

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_game2gff3.pl -i infile.xml -o outfile.gff

=head2 Required Arguments

    --infile        # Path to the input file or direcotry
    --outfie        # Path to the output file

=head1 DESCRIPTION

This program can covert game.xml files to the gff3 format. This script
is essentially a wrapper around the Apollo genome annotation editor
program, so Apollo is require for this script to work. It can convert
a single file at a time, or all game.xml files in a directory.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file or directory. If the input path is a file, the single
file will be converted to the gff3 with a file name specified by the
--outfile option. If an entire directory of files is being specified by
the -i option, all of the files ending in xml will be be converted to
gff format. The gff file will be generated in the input directory, and
the files will have the existing name with gff appended.

=item -o,--outfile

Path of the output file. If the input path is a directory, then this
argument is not needed. 

=back

=head1 OPTIONS

=over 2

=item --ap-bin

Path to the apollo binary. By default, the program assumes that the directory
containing the apollo program is included in your PATH. If this is not
the case, you will need to specify the location of the apollo program
using the -ap-bin variable.

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

=back

=head1 EXAMPLES

=head2 Convert a Single File

To convert a single game.xml file to gff3 format, you would use the following
command:

    cnv_game2gff3.pl -i HEX001.game.xml -o HEX001.gff

This will result in a single file name HEX001.gff in the path specified
by the -o argument.

=head2 Convert an Entire Directory of Files

To convert an entire directory of game.xml files to the gff3 format,
simply specify the directory containing the game.xml files with
the -i argument:

    cnv_game2gff3.pl -i annotation_dir/

This will find all files with names ending in xml, and will convert
those files to the gff3 format. The .gff extension will be added on
to the existing file name. For example, the following list shows
the name of game files and the name of the gff files that would
be created:

    Game File             GFF File Created
    HEX001.game.xml       HEX001.game.xml.gff
    HEX002.game.xml       HEX001.game.xml.gff
    HEX003.game.xml       HEX003.game.xml.gff

The gff files will be created in the same directory as the game.xml files.
The original files will not be overwritten in this process.

=head1 DIAGNOSTICS

The following information includes error message you may encounter
when using this program, and information on possible solutions
to these errors.

=head2 ERROR: No xml files were found in the input directory.

This program that the files you want to convert will be named with 
the xml extension, and that all files that end with xml are game.xml files.
If you have files in this directory that you need to convert, they must
all end with xml.

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

=over

=item Apollo Genome Annotation Curation Tool

This program uses the Apollo program to convert from game.xml to gff3
format. The Apollo program can be obtained at:
http://apollo.berkeleybop.org/current/index.html

=back

=head2 Perl Modules

This program does not make use of Perl modules beyond the standard modules
included in basic installations of Perl.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over 2

=item Apollo Required

This program is completely dependent on a local installation of the Apollo
program since it depends on Apollo to do the conversion.

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

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 02/17/2009

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 02/17/2009
# - Main body of program written
# 02/20/2009
# - POD documentation finished
