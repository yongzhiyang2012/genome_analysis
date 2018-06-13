#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_repmask2gff.pl - Convert repeatmasker output to gff   |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_sourceforge.net                   |
# STARTED: 10/30/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts the repeat masker *.out file into the GFF       |
#  format.                                                  |
#                                                           |
# USAGE:                                                    |
#  cnv_repmask2gff.pl -i infile.out -o outfile.gff          |
#                                                           |
# VERSION: Release 1.0                                      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TO DO: Better Sequence ontology feature names
#        currently just using exon for Apollo rendering

package DAWGPAWS;

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

#//////////////////////
my $file_num_max = 1;
#\\\\\\\\\\\\\\\\\\\\\\

#-----------------------------------------------------------+
# VARIABLE SCOPE                                            |
#-----------------------------------------------------------+
my $infile;                    # Input path of *.out file to convert
my $outfile;                   # Output gff file path
my $param;                     # The parameter set, db name used
my $seqname;                   # The id of the sequence annoated
my $prog = "repeatmasker";     # The program used 

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $apollo = 0;                # Path to apollo and apollo variables
my $test = 0;
my $verbose = 0;
my $debug = 0;                 # Run the program in debug mode 
my $append = 0;
my $pos_strand = 0;            # Put all results in positive strand
#my $neg_strand = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED ARGUMENTS
		    "i|infile=s"   => \$infile,
                    "o|outfile=s"  => \$outfile,
		    # ADDITIONAL OPTIONS
		    "program=s"      => \$prog,
		    "p|param=s"      => \$param,
		    "n|s|seqname|name=s"  => \$seqname,
		    # BOOLEANS
		    "plus"         => \$pos_strand,
		    "append"       => \$append,
		    "verbose"      => \$verbose,
		    "debug"        => \$debug,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);

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
    print "\ncnv_repmask2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

# DO CONVERSION
rmout2gff ($prog, $infile, $outfile, $seqname, $param, $append, $pos_strand );

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+
sub rmout2gff {

    # $source = program name source for the annotation
    # $rm_file = path to the repeat masker file
    # $gff_file = path to the gff output file
    # $seq_id   = identifier of the sequence being processed
    # $src_suffix = parameter name/ db name suffix for source
    # $do_append  = append results to an existing gff file
    # $all_plus   = place all results on the plus strand
    # the default will be to concatenate when the prefix
    # variable can not be undertood
    my ( $source, $rm_file, $gff_out, $seq_id, $src_suffix, $do_append, 
	 $all_plus) = @_;

    my $strand;

    # Set the has_seq_id flag to true is a over ride sequence id was passed
    my $has_seq_id;
    if ($seq_id) {
	$has_seq_id = 1;
    }
    else {
	$has_seq_id = 0;
    }
    
    #-----------------------------+
    # OPEN RM FILE                |
    #-----------------------------+
    # Default is to expect intput from STDIN if infile path not given
    if ($rm_file) {
	open ( RM_IN, "<".$rm_file ) ||
	    die "Could not open the RM infile:\n $rm_file\n";
    }
    else {
	print STDERR "Expecting input from STDIN";
	open ( RM_IN, "<&STDIN" ) ||
	    die "Could not open STDIN for input\n";
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
    if ($src_suffix) {
	$source = $source.":".$src_suffix;
    }

    #-----------------------------------------------------------+
    # REPEATMASKER OUT FILE CONTAINS
    #-----------------------------------------------------------+
    # 0 Smith-Waterman score of the match
    # 1 % substitutions in matching region compared to the consensus
    # 2 % of bases opposite a gap in the query sequence (deleted bp)
    # 3 % of bases opposite a gap in the repeat consensus (inserted bp)
    # 4 name of query sequence
    # 5 starting position of match in query sequence
    # 6 ending position of match in query sequence
    # 7 no. of bases in query sequence past the ending position of match
    #   in parenthesis
    # 8 match is with the Complement of the consensus sequence in the database
    #   C
    # 9 name of the matching interspersed repeat
    #10 class of the repeat, in this case a DNA transposon 
    #11 no. of bases in (complement of) the repeat consensus sequence 
    #   in parenthesis
    #12 starting position of match in database sequence 
    #   (using top-strand numbering)
    #13 ending position of match in database sequence

    while (<RM_IN>) {
	chomp;
	my @rmout = split;

	my $line_len = @rmout;
	# Skip lines related to database used etc
	next if ($line_len < 13);
	
	my $cur_strand = $rmout[8];

	#-----------------------------+
	# GET STRAND                  |
	#-----------------------------+
	if ($cur_strand =~ "C" ) {
	    $strand = "-";
	}
	else {
	    $strand = "+";
	}

	if ($all_plus) {
	    $strand = "+";
	}

	#-----------------------------+
	# SET SEQUENCE ID             |
	#-----------------------------+
	unless ($has_seq_id) {
	    if ($rmout[4]) {
		$seq_id = $rmout[4];
	    }
	    else {
		$seq_id = "seq";
	    }
	}

	#-----------------------------+
	# PRINT GFF OUT               |
	#-----------------------------+
	# TO DO: Allow for attribute name at cmd line
	# May also want to place everything in positive strand
	print GFFOUT "$seq_id\t".              # qry sequence name
	    "$source\t".              # software used
	    "exon\t".                 # attribute name
	    "$rmout[5]\t".            # start
	    "$rmout[6]\t".            # stop
	    "$rmout[0]\t".            # smith waterman score"
	    "$strand\t".              # Postive strand
	    ".\t".                    # frame
	    "$rmout[9]".              # attribute
	    "\n";
	
    }
    
    close RM_IN;
    close GFFOUT;
    
    return;
    
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

batch_mask.pl - Convert RepeatMasker output to the gff format.

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_repmask2gff.pl -i infile.out -o outfile.gff

=head2 Required Variables

    -i    # Path to the repeatmasker file to convert
    -o    # Path to the gff output file

=head1 DESCRIPTION

This program will convert the output from RepeatMasker to the standard
gff file format.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path to the intput file that contains the repeatmasker out file to
convert to the gff format. If this option is not specified, the program
will expect input from STDIN.

=item -o,--oufile

Path to the gff output file. If and outfile file path is not specified, the
program will write the output to STDERR.

=back

=head1 OPTIONS

=over 2

=item -p,--param

The parameter name to append to the source program name. This information will
be appended to the second column of the gff output file. This is used to 
specify if
you masked with the same database using a different parameter set or if you
used a different database.

=item --program

The program name to use. This is the data in the second column of the 
gff output file. Be default, this is set to 'repeatmasker'. This option
allows you to specify other program names if desired.

=item -s,--seqname

Identifier for the sequence file that was masked with repeatmasker. The
out file from repeatmasker may have truncated your original file name,
and this option allows you to use the full sequence name.

=item --plus

Write all of the annotations to be in the positive (plus) strand. By default
the program will interpret strand results reported as 'C' to be in the
negative strand orientation. Setting the --plus flag will report all
RepeatMasker hit results as occurring in the plus strand orientation.

=item --append

Append the results to an existing gff file. This must be used in
conjunctions with the --outfile option.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands.

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

=head1 EXAMPLES

=over 2

=item Typical Use

The typical use of this program will be to convert a repeat masker
output file to the gff format.

  cnv_repmask2gff.pl -i HEX3045G05_TREP.out -o HEX3045G05_TREP.gff

This will produce a GFF file similar to the following:

 HEX3045G05  repeatmasker    exon     469	493	25	-     .	AT_rich
 HEX3045G05  repeatmasker    exon     716	754	25	+     .	AT_rich
 HEX3045G05  repeatmasker    exon     1764	2069	469	+     .	TREP20
 HEX3045G05  repeatmasker    exon     1816	2105	507	+     .	TREP58
 HEX3045G05  repeatmasker    exon     1920	2248	450	+     .	TREP214
 ...

=item Identify the Database Used

It may be also be useful to run repeatmasker against a number of different
databases. You would therefore want to specify the database used in your
gff output file. This can be specified using the --param option, to specify
the database in the parameter tag. For example, if you used the TREP
database as your database for masking:

  cnv_repmask2gff.pl -i rm_result.out -o rm_resout.gff --param TREP

This will append the parameter tag 'TREP' to the source column (col 2) of 
the gff output file and 
will produce a GFF file simlar to the following:

 HEX3045G05  repeatmasker:TREP	exon   469     493     25     -	.      AT_rich
 HEX3045G05  repeatmasker:TREP	exon   716     754     25     +	.      AT_rich
 HEX3045G05  repeatmasker:TREP	exon   1764    2069    469    +	.      TREP20
 HEX3045G05  repeatmasker:TREP	exon   1816    2105    507    +	.      TREP58
 HEX3045G05  repeatmasker:TREP	exon   1920    2248    450    +	.      TREP214
 ...

=item Program Source

It may also be useful for you to specify a different program source name
depending on the needs of your individual pipeline. You can do this
using the --program option. For example, to shorten the full name
repeatmasker to 'RM', you could use the following command

  cnv_repmask2gff.pl -i rm_result.out -o rm_result.gff --program RM

This will changed the source id in the second column of the output to RM
and will result in a GFF output file similar to the following:

  HEX3045G05	RM     exon	469	493	25	-     .	AT_rich
  HEX3045G05	RM     exon	716	754	25	+     .	AT_rich
  HEX3045G05	RM     exon	1764	2069	469	+     .	TREP20
  HEX3045G05	RM     exon	1816	2105	507	+     .	TREP58
  HEX3045G05	RM     exon	1920	2248	450	+     .	TREP214
  ...

This can also be used in conjunction with the param tag:

  cnv_repmask2gff.pl -i result.out -o result.gff --program RM --param TREP

This will result in a GFF file similar to the following

 HEX3045G05   RM:TREP	exon	469	493	25	-     .	AT_rich
 HEX3045G05   RM:TREP	exon	716	754	25	+     .	AT_rich
 HEX3045G05   RM:TREP	exon	1764	2069	469	+     .	TREP20
 HEX3045G05   RM:TREP	exon	1816	2105	507	+     .	TREP58
 ...

=item Specify the Sequence Source

It is possible that the repeatmasker out file will truncate the name of
your source sequence. You can restore this to the original name using the
--name option. For example, if your full name was HEX3045G05_A001 you could
specify this as:

  cnv_repmask2gff.pl -i result.out -o result.gff --name HEX3045G05_A001

This will result in a GFF output file similar to the following:

 HEX3045G05_A001   repeatmasker	exon	469	493	25	-    .	AT_rich
 HEX3045G05_A001   repeatmasker	exon	716	754	25	+    .	AT_rich
 HEX3045G05_A001   repeatmasker	exon	1764	2069	469	+    .	TREP20
 HEX3045G05_A001   repeatmasker	exon	1816	2105	507	+    .	TREP58
 HEX3045G05_A001   repeatmasker	exon	1920	2248	450	+    .	TREP214
 ...

=back

=head1 DIAGNOSTICS

The error messages that can be generated will be listed here.

=over 2

=item * Expecting Input from STDIN

You will see this message if you did not specify an input file with the -i
or --input options.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuration file or variables set in the
user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * RepeatMasker

This program is designed to parse output file produce by the RepeatMasker
program. RepeatMaske is available for download from:
http://www.repeatmasker.org/

=back

=head2 Required Perl Modules

This program does not make use of Perl modules outside of the normal suite
of modules present in a typical installation of perl.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited RepeatMasker Version Testing.

This program has been tested with RepeatMasker v  3.1.6. If you find that this
program is not compatible with other versions of RepeatMasker, please file
a BUG report an include an output file that is not able to be parsed. 

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

STARTED: 04/10/2006

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 4/10/2006
# - Program started
# - Parsing of RepeatMasker *.out file to an Apollo 
#   compatible *.GFF file done. Working copy for both
#   a single track for each repeat class, and a data
#   track with all of the results from a Repeat Library
#   grouped toether.
# - Two dimensional array containing the repeat library
#   database informaiton.
#
# 4/11/2006
# - Additional code comments and reformat
# - Code to cycle through a set of FASTA files and store
#   the output data in a separate named folder for each
#   of the FASTA flies.
# - Added the -fixed to the command to run RepeatMasker
#   to make sure that the text format is the same for all
# - Since Apollo does not allow additional GFF files to be
#   added later, I added code to create a GFF file that has
#   the outcome for all of the RepeatMask runs in one 
#   database.
#
# 4/12/2006
# - Added the AllRepeats to the dataset to make sure that 
#   that the repeat characterization is cumulative
# - ConvertGFF2Chado subfunction
# - Adding Smith Waterman score from RepeatMasker to the
#   GFF output file. This should be the first column in the 
#   *.out file and the 6th column in the GFF file.
#
# 4/16/2005
# - Added array for list of fasta files to process and testing
#   with small set of sequences.
# - The fasta file was manually edited to use just the GB
#   ID in the FASTA header. This could be done using the 
#   the seq object PERL module from bioperl.
# - The fasta files should be read from an input directory
#   and then the output should be copied to an output dir.
#
# 4/17/2006
# - Working out variables names
#
# 4/25/2006
# - Adding ability to get the input set of FASTA files
#   from a directory 
#
# 4/26/2006
# - Working out the copy of the relevant output files to
#   a central directory for the FASTA file (ie. all 
#   programatic output from the FASTA file goes to 
#   a dir named for the file.) 
# - Added use of File::Copy module for cp commands
#
# 5/24/2006
# - Added process log to the program. This will allow
#   me to launch the program at the lab and then
#   monitor the process at home.
#
# 5/31/2006
# - Made local copy dir a variable that can be set at the top
#
# 09/26/2006
# - A few changes made to the base code to make this
#   work on the altix. Using databases at
#   /scratch/jestill/repmask/
#
# 09/28/2006
# - Changes make to get this to work on the altix
# - Changed the format of the repeat databases list
#   This is currently a two-d array
#
# 07/13/2007
# - Added POD documentation.
# - Renamed to automask.pl
# - Made this the official 1.0 release
# - Adding command line variables
# - Added print_help subfunction
# - Modified ApolloConvert subfunction to 
#   apollo_convert
# - Getting rid of the die commands and replacing
#   them with the write to logfile commands. This
#   should prvent the program from dying when
#   working on really large batch files.
# - Added cmd line options for:
#   - Number of processors
#   - Engine to use in repeatmasker
#   - apollo variable at cmd line
#     initially will be used as boolean but can be
#     used to pass apollo version, path and
#     desired output (ie. chado, game.xml etc)
#
# 07/15/2007
# - Added the ability to show help when the
#   required options $indir and $outdir are not
#   present.
#
# 07/16/2007
# - Added error message when no fasta files were
#   found in the input directory.
# - Added the ability to append a slash to the end
#   of the $indir and $outdir variables if they
#   are not already present
#
# 07/17/2007
# - Added the ability to get the search_name directly
#   from the fasta file header
# - Getting root name for the output files and folders
#   by stripping off the fasta extension with regexp
# - Changed output dir to just the root name
# - Adding rmout_to_gff subfunction to convert the repeat
#   masker output to gff format. This will get rid
#   of any dependence on grep and awk and will provide
#   a more stable program
# 
# 07/18/2007
# - rmout_to_gff working now ... strangely had to make
#   the type exon to draw properly in apollo
# - Got rid of previous code that used grep and awk
#   to do this  conversion
# - Added command line variable to specify the full
#   path of the RepeatMasker program. This should
#   help this to run machines without having to 
#   make any changes to the user's environment.
#   (ie. don't have to add RepeatMaker to user's path)
# - Renamed program again to batch_mask.pl
# 
# 09/07/2007
# - Moving POD documentation to the end of the program
# - Changing to use a config file instead of internal 2-d array
# - Modified all variable names to lowercase
# - Added use strict
#
# 09/11/2007
# - Getting rid of LOG, all printing to STDERR
# - Dropped attempts to move *.tbl file.
#
# 02/02/2009
# - Dropped the Apollo converstion subfunction
# - Updating POD documentation
# - Added print_help subfunction to extract help from 
#   POD documentation
# - Dropped log file in favor of STDERR
#
# 03/30/2009
# - Added support for --seqname in addition to --name
