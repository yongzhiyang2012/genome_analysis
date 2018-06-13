#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_repseek2gff.pl - Convert repseek output to gff format |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/24/2008                                       |
# UPDATED: 03/30/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert repseek output to a GFF format compatible with   |
#  the Apollo genome annotation program.                    |
#                                                           |
# USAGE:                                                    |
#  cnv_repseek2gff.pl -i repseek_out.txt -o repseek.gff     |
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
my $seq_id = "seq";        # Default seqname is seq
my $parameter_set = 0;     # The parameter set used, set to default

my $program = "repseek";

# BOOLEANS
my $test = 0;
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
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "s|seqname=s" => \$seq_id,
		    "program=s"   => \$program,
		    "p|param=s"   => \$parameter_set,
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
    print "\ncnv_repseek2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
# Just use the repseek2gff subfunction to do the conversion

&repseek2gff ($seq_id, $parameter_set, $infile, $outfile, $program);

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub repseek2gff {

    # Need to convert this to using hashes and use STDIN,STOUT
    # when file path arguments are not passed
    
    # $seq_id is the id of the query sequence
    # $repin is the path to the file to convert to gff
    # $repout is the path to the gff output file

    my ($seqname, $param_set, $repin, $repout, $source) = @_;

    #-----------------------------+
    # OPEN FILE HANDLES           |
    #-----------------------------+
    # DEFAULT TO STDIN
    if ($repin) {
	open (REPIN, "<$repin") ||
	    die "Can not open input file:\n$repin\n";
    } 
    else {
	print STDERR "Expecting input from STDIN\n";
	open (REPIN, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }

    # DEFAULT TO STDOUT
    if ($repout) {
	open (GFFOUT, ">$repout") ||
	    die "Can not open output file:\n$repout\n";
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    my $palindrome_id = 0; # Running palindrome id number
    my $tandem_id = 0;
    my $close_id = 0;
    my $overlap_id = 0;
    my $interseq_id = 0;
    my $repseek_id = 0;

    my $inv_id = 0;        # Counter/id for inverted repeats
    my $dir_id = 0;        # Counter/id for direct repeats

#    my $source;            # The data for the source column
#    my $source = "repseek";

    if ($param_set) {
	$source = $source.":".$param_set;
    }


    # MAY CONSIDER JUST RETURNING THE TANDEM REPEATS ?
    while (<REPIN>) {
	chomp;

	$repseek_id = $repseek_id + 1;
	my @rep_parts = split;
	my @type_parts = split(/\./, $rep_parts[0] );

	my $repeat_type = $type_parts[0];
	my $repeat_direction = $type_parts[1];

	my $copy1_start = $rep_parts[1];
	my $copy2_start = $rep_parts[2];
	my $copy1_end = $copy1_start + $rep_parts[3];
	my $copy2_end = $copy2_start + $rep_parts[4];
	my $alignment_score = $rep_parts[8];

	my $attribute = "repseek".$repseek_id.":$repeat_direction";

	# Will split the repseek source into inverted repeats 
	# vs direct repeats .. putative TEs

	# Source sets
	# - Overlap
        # - Palindromes
	my $feature;
	#$source = "repseek";



	# Translate repeat direction
	if ($repeat_direction =~ "inv") {
	    $repeat_direction = "inverted";
	}
	elsif ($repeat_direction =~ "dir") {
	    $repeat_direction = "direct";
	}

	if ($repeat_type =~ "Overlap") {
	    $feature = "overlapping_"."$repeat_direction"."_repeat";
	}
	else {
	    $feature = "$repeat_direction"."_repeat";
	}

	if ($repeat_type =~ "Palindrome") {
	    $feature = "palindromic_repeat";
	}

	# parameter set will go in the attribute column

	#-----------------------------+
	# PRINT OUTPUT TO GFF         |
	#-----------------------------+
	print GFFOUT "$seqname\t".       # Seqname
	    "$source\t".                 # Source
	    "$feature\t".                    # Features: Apollo kluge
	    "$copy1_start\t".            # Start
	    "$copy1_end\t".              # End
	    "$alignment_score\t".        # Score
	    "+\t".                       # Strand
	    ".\t".                       # Frame
	    "$attribute".                # Attribute
	    "\n";

	print GFFOUT "$seqname\t".       # Seqname
	    "$source\t".                 # Source
	    "$feature\t".                    # Features: Apollo kluge
	    "$copy2_start\t".            # Start
	    "$copy2_end\t".              # End
	    "$alignment_score\t".        # Score
	    "+\t".                       # Strand
	    ".\t".                       # Frame
	    "$attribute".                # Attribute
	    "\n";

    }

    close (REPIN);
    close (GFFOUT);


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

cnv_repseek2gff.pl - Convert repseek output to gff format

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_repseek2gff.pl -i InFile -o OutFile

=head2 Required Options

    --infile        # Path to the repseek output file
                    # Assumes STDIN if not givien
    --outfile       # Path to the output gff file
                    # Assumes STDOUT if not given

=head1 DESCRIPTION

Convert repseek output to a GFF format compatible with
the Apollo genome annotation program.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If an input file is not provided, the program
will expect input from STDIN.

=item -o,--outfile

Path of the output file. If an output path is not provided, the program
will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item -s,--seqname

The name of the sequence contig that is being annotated. This will be used
for the first column in the gff file. If this option is not specified
the name will default to 'seq'.

=item -p, --param

The name of the parameter set used in repseek. This allows the user to
define multiple parameter sets in repseek, and then draw them as
separate tracks in annotation visualization programs.

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

=item -v, --verbose

Run the program in verbose mode.

=back

=head1 DIAGNOSTICS

=over

=item Can not open input file

The input file path provided by the -i, --infile switch is not valid.

=item Can not open output file

The output file path provided by the -o, --outfile witch is not valid.

=item Expecting input from STDIN

When an input file path is not specified, the program expects input 
to come from STDIN. This default behaviour allows output from repseek
to be piped directly into the cnv_repseek2gff.pl program.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of configuration files or variables
set in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item RepSeek

This program is designed to parse output from the RepSeek program. RepSeek
is available from:
http://wwwabi.snv.jussieu.fr/~public/RepSeek/

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

=item * Limited RepSeek testing

This program has been tested and know to work with RepSeek version 10May2007.
Other versions have not been tested for compatibility.

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

STARTED: 11/25/2008

UPDATED: 03/30/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 11/24/2008
# - The earlier version of this program was lost so
#   starting this program de novo.
# - This program accepts input from STDIN and will output 
#   to STDOUT.
#
# 01/20/2009
# - Added svn revision tracking
# 
# 03/30/2009
# - changed source to not include repeat direction
# - changed feature name from exon to sequence ontology
#   complient names for repeats
#      -inv -> inverted_repeat
#      -dir -> direct_repeat
#      -overlapping_direct_repeat   
#         -- not SeqOntology complient
#      -overlapping_inverted_repeat
#         -- not SeqOntology complient
#      -palindromic_repeat
#  - Added support for --param
#  - Added support for --program
