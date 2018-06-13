#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnf_findltr2gff.pl - Converts find_ltr output to gff      |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/13/2007                                       |
# UPDATED: 01/29/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts output from the find_ltr.pl program to gff      |
#  format for easy import into apollo.                      |
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
my $infile;
my $outfile;
my $inseqname = "seq";         # Default name is seq
my $findltr_suffix;            # Parameter name
my $program = "FindLTR";

# Booleans
my $do_gff_append = 0;
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "s|seqname=s" => \$inseqname,
		    # ADDITIONAL OPTIONS
		    "program=s"   => \$program,
		    "p|param=s"   => \$findltr_suffix,
		    "append"      => \$do_gff_append,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);


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
    print "\ncnv_findltr2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# MAIN PROGRAM BODY           |
#-----------------------------+
if ($do_gff_append) {
    findltr2gff ( $program, $infile, $outfile, 1, $inseqname, $findltr_suffix);
}
else {
    findltr2gff ( $program, $infile, $outfile, 0, $inseqname, $findltr_suffix);
}

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
    my ($gff_source, 
	$findltr_in, $gff_out, $append_gff, $seqname, $gff_suffix) = @_;

    # find_ltr
    #my $gff_source;                 # 
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
    if ($findltr_in) {
	open (INFILE, "<$findltr_in") ||
	    die "Can not open input file:\n$findltr_in\n";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }


    if ($gff_out) {
	if ($append_gff) {
	    open (GFFOUT, ">>$gff_out") ||
		die "Could not open output file for appending\n$gff_out\n";
	}
	else {
	    open (GFFOUT, ">$gff_out") ||
		die "Could not open output file for output\n$gff_out\n";
	} # End of if append_gff
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    # Append parameter name if passed
    if ($gff_suffix) {
	$gff_source = $gff_source.":".$gff_suffix;
    }

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
	    
	    # FULL LTR Retrotransposon Span
	    print GFFOUT "$seqname\t". # Name of sequence
		"$gff_source\t".       # Source
		"LTR_retrotransposon\t".
		"$ltr5_start\t".       # Feature start
		"$ltr5_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".         # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)

	    # 5'LTR
	    print GFFOUT "$seqname\t". # Name of sequence
		"$gff_source\t".       # Source
		"five_prime_LTR\t".
		"$ltr5_start\t".       # Feature start
		"$ltr5_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".         # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)

	    # 3'LTR
	    print GFFOUT "$seqname\t". # Name of sequence
		"$gff_source\t".       # Source
		"three_prime_LTR\t".
		"$ltr3_start\t".       # Feature start
		"$ltr3_end\t".	       # Feature end
		".\t".                 # Score, Could use $ltr_similarity
		"$ltr_strand\t".         # Strand
		".\t".                 # Frame
		"$findltr_name\n";     # Features (name)

	} # End of if num_in is 10

    } # End of while INFILE


} # End of findltr2gff

=head1 NAME

cnv_findltr2gff.pl - Convert output from find_ltr.pl to gff format.

=head1 VERSION

This documentation refers to version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_findltr2gff.pl -i infile.ltrpos -o outfile.gff [--seqname HEX451]

=head2 Required Arguments

    -i         # Path to the ltrpos input file 
               # Defaults to STDIN if not specified
    -o         # Path to the gff format output file
               # Defaults to STDOUT if not specified
    -s         # Identifier for the annotated sequence
    -p         # Name of use as a suffix in the gff source column
               # This is generally the parameter set used

=head1 DESCRIPTION

Converts output from the find_ltr.pl ltr annotation program to 
the standard gff format. The find_ltr.pl program is a component of
the LTR DeNovo package 
( http://darwin.informatics.indiana.edu/evolution/LTR.tar.gz ) 
described in Rho, M., J. H. Choi, et al. (2007). BMC Genomics 8: 90.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. This is a computational result from the find_ltr.pl
program.  If an input file is not provided, the program
will expect input from STDIN.

=item -o,--outfile

Path of the output file.  If an output path is not provided,
the program will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item --program

The program name used to produce the computational prediction. By default
this is set to be 'FindLTR'.

=item -p, --param

Name to use as the suffix when naming the source in the gff output. This
can be used to distinguish among find_ltr runs that used different paramater
sets. This will be append to findltr: to produce a name that can be
unique for each parameter set result. For example the default parameter set
could be given a I<def> suffix while and alternate paramater set could
be given the suffix I<alt>. This will produce gff source lines of
I<findltr:def> and I<findltr:alt>.

=item -s,--seqname

The name of the sequence that is being annotated. This is placed in the 
first field of the gff output file.

=item --apend

Append the GFF output to the gff file indicate by the --outfile option. This
allows you to append results from multiple programs to a single GFF file, but
it is generally recommended to create separate GFF files and concatenate
them at a later stage.

=item -q,--quiet

Run the program with minimal output.

=item --verbose

Run the program in verbose mode.

=back

=head1 DIAGNOSTICS

Error messages that you may encounter and possible solutions are listed below:

=over 2

=item Expecting input from STDIN

If a file is not specified by the -i or --infile option, the program will
expect to receive intput from standard input.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuration file or variables set
in the user's environment.

=head1 DEPENDENCIES

=head2 Software

The following software is required for this program:

=over 2 

=item *ltr_finder Program

This program is designed to parse output from the ltr_finder.pl
program. http://darwin.informatics.indiana.edu/evolution/LTR.tar.gz

M Rho et al. 2007. I<'De novo identificatio of LTR retrotransposons
in eukaryotic genomes'> BMC Genomics 8:90.

=back

=head2 Perl Modules

This program does not make use of perl modules beyond those installed
with the basic Perl package. If you discover a dependency that is not
documented here, please email the author or file a bug report.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over 2

=item * Works with 2008 version of find_ltr.pl

This program can parse results generated by the find_ltr.p program 
as downloaded from
http://darwin.informatics.indiana.edu/evolution/LTR.tar.gz
in 2008.

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

=head1 HISTORY

STARTED: 09/13/2007

UPDATED: 01/29/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 09/13/2007
# - Base program written
# 09/14/2007
# - Added the suffix to the gff_source variable this allows
#   for different parameter sets to recieve unique names in
#   the gff output file.
# 01/29/2009
# -Updating POD documentation
# -Adding input from STDIN output to STDOUT
# -Added help function to extract help from POD doc
# -Dropped required variables (infile,outfile,seqname
#  replaced with STDIN,STDOUT and seq
# -Replaced exon if gff names with names of objects

