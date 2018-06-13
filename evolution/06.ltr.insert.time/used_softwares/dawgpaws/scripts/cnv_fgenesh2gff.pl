#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_fgenesh2gff.pl - Convert fgenesh output to gff        |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 01/31/2009                                       |
# UPDATED: 03/27/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert output from fgenesh to the gff format.           |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
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
use Bio::Tools::Fgenesh;

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;
my $outfile;
my $seqname;
my $param;
my $tmp_file_path;            # A temp file stripped of html tags
my $prog;                     # The program used to generate the data

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode
my $append = 0;
my $test = 0;
my $strip_html = 0;               # Attempt to strip html

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "p|param=s"   => \$param,
		    "program=s"   => \$prog,
		    # Allow name in addition to seqname
		    "s|name|seqname=s" => \$seqname,
		    "html"        => \$strip_html,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "append"      => \$append,
		    # ADDITIONAL INFORMATION
		    "usage"       => \$show_usage,
		    "test"        => \$test,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,);

# MAY NEED TO ITERATE ACROSS THE FILE AND GET RID OF COPYWRITE STATEMENT

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
    print "\ncnv_fgenesh2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# TAKE A LOOK AT THE FILE TO  |
# SEE IF IT IS HTML           |
#-----------------------------+
unless ($strip_html) {
# OPEN INPUT FILE HANDLE
    if ($infile) {
	open (TMPIN, "<$infile") ||
	    die "Can not open input file ";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (TMPIN, "<STDIN") ||
	    die ""
    }
    
    while (<TMPIN>) {
	if ($_ =~ m/DOCTYPE HTML PUBLIC/) {
	    print STDERR "------------------------------------------------\n";
	    print STDERR " WARNING: The input file appears to be HTML\n";
	    print STDERR " Attempting to strip HTML from text file\n";
	    print STDERR "------------------------------------------------\n";
	    $strip_html = 1;
	}
    }
    close TMPIN;
}

#-----------------------------+
# STRIP HTML                  |
#-----------------------------+
# A very simple attempt to strip HTML from the output
if ($strip_html) {

    # OPEN INPUT FILE HANDLE
    if ($infile) {
	open (TMPIN, "<$infile") ||
	    die "Can not open input file ";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (TMPIN, "<STDIN") ||
	    die ""
    }

    # OPEN OUTPUT FILE HANDLE
    if ($infile) {
	$tmp_file_path = "$infile.strip.tmp";
    }
    else {
	$tmp_file_path = "fgenesh.strip.tmp";
    }
    open (TMPOUT, ">$tmp_file_path") ||
	die "Can not write to temp file:\n$tmp_file_path\n";

    # STRIP HTML
    while (<TMPIN>) {
	next if m/^\</;              # Remove lines starting with <
	next if m/www\.softberry/;   # Remove copywrite
	s /\&gt\;/\>/;               # Replace &gt; with >
	#print STDERR $_;
	print TMPOUT $_;
    }
    
    close TMPIN;
    close TMPOUT;
}


#-----------------------------+
# DO THE CONVERSION           |
#-----------------------------+
#my $prog = "fgenesh";
if ($strip_html) {
    fgenesh2gff ($prog, $tmp_file_path, $outfile, $seqname, $param, $append);
}
else {
    fgenesh2gff ($prog, $infile, $outfile, $seqname, $param, $append);
}

exit 0;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub fgenesh2gff {
    # fgnesh_in  - path to the fgenesh program
    # gff_out    - path to the gff output file
    # seq_id     - id of the source sequence
    # src_suffix - parameter id for fgenesh run

    my ($source, $fgenesh_in, $gff_out, $seq_id, $src_suffix, $do_append ) = @_;

    #-----------------------------+
    # OPEN THE FGENESH INFILE     |
    #-----------------------------+
    my $fgenesh_result;
    if ($fgenesh_in) {
	$fgenesh_result = Bio::Tools::Fgenesh->new(-file => $fgenesh_in);
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	$fgenesh_result = Bio::Tools::Fgenesh->new( -fh  => \*STDIN );
    }

    #-----------------------------+
    # OPEN THE GFF OUTFILE        |
    #-----------------------------+
     # Default to STDOUT if no argument given
    if ($gff_out) {
	if ($do_append) {
	    open (GFFOUT, ">>$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	}
	else {
	    open (GFFOUT,">$gff_out") ||
		die "ERROR: Can not open gff outfile:\n $gff_out\n";
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    #-----------------------------+
    # SET PROGRAM SOURCE          |
    #-----------------------------+
    unless ($source) {
	$source = "fgenesh";
    }
    if ($src_suffix) {
	$source = $source.":".$src_suffix;
    }

    my $gene_num = 0;
    while (my $gene = $fgenesh_result->next_prediction()) {

	$gene_num++;
	#-----------------------------+
	# SET SEQUENCE ID             |
	#-----------------------------+
	unless ($seq_id) {
	    if ($gene->seq_id()) {
		$seq_id = $gene->seq_id();
	    }
	    else {
		$seq_id = "seq";
	    }
	}

	# $gene is an instance of Bio::Tools::Prediction::Gene, which inherits
	# off Bio::SeqFeature::Gene::Transcript.
	#
	# $gene->exons() returns an array of 
	# Bio::Tools::Prediction::Exon objects
	# all exons:
	my @exon_arr = $gene->exons();
	

	foreach my $ind_exon (@exon_arr) {
	    #print STDERR $ind_exon;
#	    if ($ind_exon->is_coding()) {
#		print STDERR "coding\t";
#	    }
#	    else {
#	    }

	    #-----------------------------+
	    # FORMAT STRAND               |
	    #-----------------------------+
	    my $strand = $ind_exon->strand()."\t";
	    if ($strand =~ "-1") {
		$strand = "-";
	    }
	    else {
		$strand = "+";
	    }

	    #-----------------------------+
	    # GET START AND END           |
	    #-----------------------------+
	    my $start = $ind_exon->start();
	    my $end = $ind_exon->end();

	    if ($start > $end) {
		$end =  $ind_exon->start();
		$start = $ind_exon->end();
	    }

	    #-----------------------------+
	    # PRINT GFF OUTPUT            |
	    #-----------------------------+
	    # The following prints one line at a time
	    #print GFFOUT $seq_id."\t";     # Seqname
	    #print GFFOUT $source."\t";     # Source
	    #print GFFOUT "exon\t";         #feature
	    #print GFFOUT $start."\t";      # start
	    #print GFFOUT $end."\t";        # end
	    #print GFFOUT $ind_exon->score()."\t";  # score
	    #print GFFOUT $strand."\t";       # strand
	    #print GFFOUT ".\t";              # frame
	    #print GFFOUT "gene_".$gene_num;  # attribute
	    #print GFFOUT "\n";

	    print GFFOUT $seq_id."\t".     # Seqname
		$source."\t".     # Source
		"exon\t".         #feature
		$start."\t".      # start
		$end."\t".        # end
		$ind_exon->score()."\t".  # score
		$strand."\t".       # strand
		".\t".              # frame
		"gene_".$gene_num.  # attribute
		"\n";


	    # The following does work
	    #print GFFOUT $ind_exon->primary_tag()."\t";


	    #///////////////////////////////
	    # The following do not work
	    #///////////////////////////////
	    #print GFFOUT $gene->cds()."\n";
	    #print GFFOUT $ind_exon->significance()."\t";
	    #print $ind_exon->predicted_cds();
	    #print GFFOUT $ind_exon->coding_signal_score()."\t";
	    #print GFFOUT $ind_exon->seq_id()."\t";
	    #print $ind_exon->significance()."\n";
	    # Get the CDS of the sequence
	    #print STDERR $ind_exon->cds()."\n";

	}
	
#       # initial exons only
#       @init_exons = $gene->exons('Initial');
#       # internal exons only
#       @intrl_exons = $gene->exons('Internal');
#       # terminal exons only
#       @term_exons = $gene->exons('Terminal');
#       # singleton exons: 
#       ($single_exon) = $gene->exons();


   }

    # CLOSE FGENESH
    $fgenesh_result->close();
    

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

cnv_fgenesh2gff.pl - Convert fgenesh gene predictions to gff format

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_fgenesh2gff.pl -i infile.txt -o outfile.gff

=head2 Required Arguments

    --infile        # Path to fgenesh result to convert
    --outfie        # Path to the gff format output

=head1 DESCRIPTION

This program converts output from the fgenesh program to the gff format. If 
the fgenesh output file appears to be saved from the web, the program
will attempt to first strip the HTML tags from the text before converting
to the GFF format.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. This should a text file of the result of the fgenesh
gene prediction program. If an input file is not specified, then the program
will expect input from STDIN.

=item -o,--outfile

Path of the gff file that is produced by the program. If an output file
is not specified, the program will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item --html

Use this to convert the output from the softberry website if you
saved the text in html format.

=item -p,--param

The label used to describe the parameter set used for the the annotation
program. This identifier will be appended the source column (col 2)
in the GFF output.

=item -s,--seqname

This is the name of the sequence that was annotated. This will be used
in the source column (col 1) of the gff output file. By default, the program
will use the name of the sequence as specified in the fgenesh output file.

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

=head1 EXAMPLES

=over 2

=item Typical Use

Typically you will be using this program to convert the fgenesh annotation 
output for an individual sequence file to the gff format.

  cnv_fgenesh2gff.pl -i fgenesh_result.txt -o fgenesh_result.gff

This will result in a GFF result similar to the following:

 HEX3045G05  fgenesh   exon	961	1456	12.02	+     .	gene_1
 HEX3045G05  fgenesh   exon	1702	2725	0.04	+     .	gene_1
 HEX3045G05  fgenesh   exon	3619	3982	10.41	+     .	gene_1
 HEX3045G05  fgenesh   exon	6960	7273	13.70	+     .	gene_2
 HEX3045G05  fgenesh   exon	7435	7789	21.29	+     .	gene_2
 HEX3045G05  fgenesh   exon	7904	8091	1.14	+     .	gene_2
 HEX3045G05  fgenesh   exon	8248	9163	13.79	+     .	gene_2
 HEX3045G05  fgenesh   exon	9206	9587	8.00	+     .	gene_2
 ...

=item Specify the Sequence ID

Generally the cvn_fgenesh2gff.pl program will use the label for the sequence
as reported in the fgenesh report file. Otherwise, you can specify the
source sequence name using the -n or --name flag. For example:

  cnv_fgenesh2gff.pl -i result.txt -o result.gff -n wheat_1

Will result in a gff file like the following:

 wheat_1    fgenesh	exon	961	1456	12.02	+     .	gene_1
 wheat_1    fgenesh	exon	1702	2725	0.04	+     .	gene_1
 wheat_1    fgenesh	exon	3619	3982	10.41	+     .	gene_1
 wheat_1    fgenesh	exon	6960	7273	13.70	+     .	gene_2
 wheat_1    fgenesh	exon	7435	7789	21.29	+     .	gene_2
 wheat_1    fgenesh	exon	7904	8091	1.14	+     .	gene_2
 wheat_1    fgenesh	exon	8248	9163	13.79	+     .	gene_2
 wheat_1    fgenesh	exon	9206	9587	8.00	+     .	gene_2
 ...

This option allows you to change the name of the sequence source without
having to run the fgenesh program again.

=item Specify the Parameter Set

It is often useful to run a program using different parameter sets. The
cnv_fgenesh2gff.pl program therefore allows you to specify the label
for a set of parameters to be able to distinguish multiple prediction results
from the same program using different parameter combinations. This
parameter set label will be added to the second column of the gff
output file.

For example running the program with parameter set one:

  cnv_fgenesh2gff.pl -i result.txt -o result.gff -p set_1

This will result in a GFF file like the following:

 HEX3045G05  fgenesh:set_1  exon   961    1456	12.02	+    .	gene_1
 HEX3045G05  fgenesh:set_1  exon   1702   2725	0.04	+    .	gene_1
 HEX3045G05  fgenesh:set_1  exon   3619   3982	10.41	+    .	gene_1
 HEX3045G05  fgenesh:set_1  exon   6960   7273	13.70	+    .	gene_2
 HEX3045G05  fgenesh:set_1  exon   7435   7789	21.29	+    .	gene_2
 HEX3045G05  fgenesh:set_1  exon   7904   8091	1.14	+    .	gene_2
 HEX3045G05  fgenesh:set_1  exon   8248   9163	13.79	+    .	gene_2
 HEX3045G05  fgenesh:set_1  exon   9206   9587	8.00	+    .	gene_2
 ...

Then running the program wit parameter set two:

  cnv_fgenesh2gff.pl -i result.txt -o result.gff -p set_2

This will result in a GFF file like the following:

 HEX3045G05  fgenesh:set_2  exon   961    1456	12.02	+    .	gene_1
 HEX3045G05  fgenesh:set_2  exon   1702   2725	0.04	+    .	gene_1
 HEX3045G05  fgenesh:set_2  exon   3619   3982	10.41	+    .	gene_1
 ...

This will allow you to later distinguish between the result for parameter
set one and the parameter set two results.

=item Accepting Input from STDIN

It is often useful in working at the unix command line to pipe the output 
from one program to another. For that reason, the cnv_fgenesh2gff.pl program
can accept input from STDIN. For example, given a text file named result.txt.
You can send the result to cnv_fgenesh2gff.pl using the cat command and
then the pipe '|':

  cat result.txt | cnv_fgenesh2gff.pl

Since an output file is not specified, the result will be printed to
STDOUT and will appear on the screen.

=item Writing Output to STDOUT

Since the program can write output to STDOUT, it is possible to
directly load the GFF file to your database. For example, if you
have a script called load_gff2mydb.pl, you can pipe the GFF results
to this program directly:

  cnv_fgenesh2gff.pl -i result.txt | load_gff2mydb.pl

This will load the result to your database without generating a copy
of the GFF file on your hard drive.

=item Removing Text That Throws Warnings

Saving the output from the fgenesh webpage will included the copywrite statement
from. You will get the following warning:

  --------------------- WARNING ---------------------
  MSG: seq doesn't validate, mismatch is ?1999,2009,<,://,/>
  ---------------------------------------------------

This warning is only written to STDERR, and should not affect the gff output
of the program. However, you can remove the offending line of fgenesh
output using the grep command before piping the text to the
cnv_fgenesh2gff.pl program.

  grep -v 'www.softberry.com' fgenesh.txt | cnv_fgenesh2gff.pl 

=item Strip HTML Tags

It is also possible to parse output from the softberry website if
it was saved in html text format using the --html option. This
will attempt to strip the html and save a local tmp copy that
is in plain text file that will then be parsed:

  cnv_fgenesh2gff.pl -i infile.txt --html

=back

=head1 DIAGNOSTICS

The following lists some typical error messages and solutions:

=over 2

=item * MSG: seq doesn't validate, mismatch is ?1999,2009,<,://,/>

This generally will be seen when the fgenesh text file includes the
copywrite statement from the softberry web site.

=item * MSG: seq doesn't validate, mismatch is &,;,:[,]13,(,)961,39821884,,,+,

You may see something like this if you are trying to parse a result
you saved from the softberry web site in the html format. The solution 
to this problem is to save the program as text. You can strip the 
html from the program using the --html option.

=item *  WARNING: The input file appears to be HTML

You will see this message if the program detects that the fgenesh output
you are trying to parse is in HTML format. If this is the case, 
cnv_fgenesh2gff.pl will attempt to save a copy of the fgenesh
result as a normal text file before converting to GFF format.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuartion file or varaibles
defined in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over 2

=item * Fgenesh

This program is designed to parse ab initio gene annotation results generated
by the Fgenesh program. These results can be generated from a local copy
of the Fgenesh program, or can be results obtained by the Fgenesh web
service provided by softberry
http://linux1.softberry.com/berry.phtml

=back

=head2 Required Perl Modules

=over 2

=item * Bio::Tools::Fgenesh

This program requires the perl module Bio::Tools::Fgenesh. This module is
part of the bioperl package

=back

Other modules or software that the program is dependent on.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over 2

=item * Not Tested on Fgenesh Binary

This progarm has only been tested with output from the softberry
website and has not been tested with the Fgenesh binary. If you find
that this program does work with the standalone program, please
contact the author and let me know.

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

STARTED: 01/31/2009

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 01/31/2009
# - Program started, basic conversion to gff using the
#   bioperl module for parsing fgenesh programs
# 02/02/2009
# - Updated POD documentation
# - Added the --html flag to strip html if needed
# - Added autodetect of html format
#
# 03/27/2009
# - Working on adding code to deal with (c) statement
# - Added --program option
# - Renamed --name option to --seqname
