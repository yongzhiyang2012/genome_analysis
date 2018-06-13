#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ltrseq2gff.pl - Convert LTR_seq output to gff format. |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/03/2007                                       |
# UPDATED: 03/30/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts LTR_seq output to gff file format.              |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TODO: Extract only the unique LTR Predictions 

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
my $infile;                    # Path of the LTR_seq file to convert
my $outfile;                   # Path to the gff format file created
my $inseqname;                 # Name of the sequence analyzed
my $parameter;                 # The parameter set used
my $program="LTR_seq";         # The program name                     
# Optional
my $infasta;                   # Path to the FASTA file analyzed

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff_append = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# Required options
                    "i|infile=s"    => \$infile,
                    "o|outfile=s"   => \$outfile,
		    "s|seqname=s"   => \$inseqname,
		    # Additional Optons
		    "program=s"     => \$program,
		    "param=s"       => \$parameter,
		    #"f|fastafile=s" => \$infasta,
		    "append"        => \$do_gff_append,
		    "q|quiet"       => \$quiet,
		    "verbose"       => \$verbose,
		    # Additional information
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

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
    print "\ncnv_ltrseq2gff.pl:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$inseqname) ) {
    print "\a";
    print STDERR "ERROR: A sequence name was not specified at the".
	" command line\n" if (!$inseqname);
    print_help ("usage", $0 );
}

#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

if ($do_gff_append) {
    ltrseq2gff ( $infile, $outfile, 1, $inseqname, $program, $parameter);
} 
else {
    ltrseq2gff ( $infile, $outfile, 0, $inseqname, $program, $parameter);
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

sub ltrseq2gff {

    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    # $source is the program
    my ($ltrseq_in, $gff_out, $append_gff, $seqname,
	$source, $param) = @_;

    #-----------------------------+
    # LTRSEQ VARS                 |
    #-----------------------------+
    my $ltrseq_infile;              # Input file that was analyzed by LTR_seq
    my $ltrseq_numseqs;             # Number of seqs analyzed by LTR_seq
    my $ltrseq_msrt;                # Max Score Ration Threshold
    my $ltrseq_ltrmin;              # LTRmin
    my $ltrseq_ltr_minexact;        # LTR Min Exact Match
    my $ltrseq_dmin;                # Dmin
    my $ltrseq_dmax;                # Dmax

    # Postion of the LTR Retrotransposon features
    my $ltrseq_name;                # Unique name assigned to the prediction
    my $ltrspan_start;              # Start of the entire LTR Retrotransposon
    my $ltrspan_end;                # End of the entire LTR Retrotransposon
    my $ltrspan_len;                # Length of the LTR span
    my $ltr5_start;                 # Start of the 5' LTR
    my $ltr5_end;                   # End of the 5' LTR
    my $ltr3_start;                 # Start of the 3' LTR
    my $ltr3_end;                   # End of the 3' LTR
    my $ltr_len;                    # Length of the LTRs
    my $mid_start;                  # Start of the LTR Mid region
    my $mid_end;                    # End of the LTR Mid region
    my $ltr_diff;                   # Percent Difference in the LTRs
    my $ltr_conf;                   # Confidence score for the prediction
    my $ltr_tsr;                    # Target site rep

    my @in_split = ();              # Split of the infile line
    my $num_in;                     # Number of split vars in the infile

    # Initialize Counters
    my $ltrseq_num = 0;             # ID Number of putatitve LTR retro (Incremented Number)

    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+
    
    if ($ltrseq_in) {
	open (INFILE, "<$ltrseq_in") ||
	    die "ERROR: Can not open input file:\n$ltrseq_in\n";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }
    
    if ($gff_out) {
	if ($append_gff) {
	    open (GFFOUT, ">>$gff_out") ||
		die "ERROR: Can not open output file for appending\n$gff_out\n";
	}
	else {
	    open (GFFOUT, ">$gff_out") ||
		die "ERROR: Can not open output file for output\n$gff_out\n";
	} # End of if append_gff
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    if ($param) {
	$source = $source.":".$param
    }

    #-----------------------------+
    # PROCESS INFILE              |
    #-----------------------------+
    while (<INFILE>) {
	chomp;

	# Print the INFILE
	#print "$_\n";

	if (m/^Report:/) {

	    # If this is a Report line, accepted or rejected can be
	    # determined by counting the line number
	    # 15 is a rejected line
	    # 19 is an accepted line 
	    my @in_split = split;
	    my $num_in = @in_split;   

	    #print "Report Line: $num_in parts\n";
	    #print "\t$_\n";

	    if ($num_in == 15) {
		#print "\tREJECTED\n";
	    }
	    #-----------------------------+
	    # ACCEPTED LTR RETRO Parts    | 
	    #-----------------------------+
	    elsif ($num_in == 19) {
		#print "\tACCEPTED\n";
		print STDERR "\n".$_."\n" if $verbose;

		$ltrseq_num++;              # Increment the count

		$ltr5_start = $in_split[4] || "0";   # Replaced with zero
		$ltr5_end = $in_split[5] || "NULL";
		$ltr3_start = $in_split[7] || "NULL";
		$ltr3_end = $in_split[8] || "NULL";
		$mid_start = $ltr5_end + 1;
		$mid_end = $ltr3_start - 1;
		$ltrspan_start = $ltr5_start;
		$ltrspan_end = $ltr3_end;
		$ltrspan_len = $in_split[12];
		$ltr_len = $in_split[10];
		$ltr_diff = $in_split[14];
		$ltr_conf = $in_split[16];
		$ltr_tsr = $in_split[18];

		$ltrseq_name = "ltrseq_".$ltrseq_num;

		if ($ltr_tsr =~ m/\#(.*)\#/) {
		    $ltr_tsr = $1;
		}

		# SHOW INFO IF VERBOSE
		
		if ($verbose) {
		    print STDERR "$ltrseq_name\n";
		    print STDERR "\tSTART:    $ltrspan_start\n";
		    print STDERR "\tEND:      $ltrspan_end\n";
		    print STDERR "\tLTRLEN:   $ltr_len\n";
		    print STDERR "\tSPAN_LEN: $ltrspan_len\n";
		    print STDERR "\tLTSR:     $ltr_tsr\n";
		    print STDERR "\tLTRCON:   $ltr_conf\n";
		    print STDERR "\tLTRDIF:   $ltr_diff\n";
		} # End of if verbose

		#-----------------------------+
		# PRINT TO GFF OUTPUT FILE    |
		#-----------------------------+

		# LTR RETRO PREDICTION SPAN
		print GFFOUT "$seqname\t".     # Name of sequence
		    "$source\t".               # Source
		    "LTR_retrotransposon\t".   # Feature, exon for Apollo
		    "$ltrspan_start\t".        # Start of the ltr span
		    "$ltrspan_end\t".          # End of the ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)
		
		# 5' LTR
		print GFFOUT "$seqname\t".     # Name of sequence
		    "$source\t".               # Source
		    "five_prime_LTR\t".        # Feature, exon for Apollo
		    "$ltr5_start\t".           # Start of the 5'ltr
		    "$ltr5_end\t".             # End of the 5' ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)

		# 3' LTR
		print GFFOUT "$seqname\t".     # Name of sequence
		    "$source\t".               # Source
		    "three_prime_LTR\t".       # Feature, exon for Apollo
		    "$ltr3_start\t".           # Start of the 3'ltr
		    "$ltr3_end\t".             # End of the 3' ltr span
		    "$ltr_conf\t".             # Score, LTR Confidence Score
		    ".\t".                     # Strand
		    ".\t".                     # Frame
		    "$ltrseq_name\n";          # Features (Name)
	    }


	    # 15 Parts is Rejected

	}
	# HEADER INFORMATION
	elsif (m/^Info:/) {
	    
	}
	# The input string that was passed
	elsif (m/^Input:/) {

	}

    } # End of while INFILE

}


1;
__END__

# Deprecated print_help
sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"MyProg.pl -i InFile -o OutFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --infile       # Path to the input file\n".
	"  --outfile      # Path to the output file\n".
	"\n".
	"OPTIONS::\n".
	"  --version      # Show the program version\n".     
	"  --usage        # Show program usage\n".
	"  --help         # Show this help message\n".
	"  --man          # Open full program manual\n".
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

cnv_ltrseq2gff.pl - Convert LTR_seq output to GFF format.

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_ltrseq2gff.pl -i InFile -o OutFile -s SeqName

=head2 Required Arguments

    -i,--infile   # Directory of fasta files to process
    -o,--outfile  # Path to the output file
    -s,--seqname  # ID of the sequence record analyzed

=head1 DESCRIPTION

Converts program output from LTR_seq LTR prediction program to GFF format.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file which is the prediction results from LTR_seq.
If an input file is not specified, the program will expect intput from
STDIN.

=item -o,--outfile

Path of the output file which will be the GFF output from cnv_ltrseq2gff.pl.
If an input file is not specified, the program will print output
to STDOUT.

=item -s,--seqname

The sequence ID of the query sequence that was analyzed. This will be
used to write to the first column of the GFF output file.

=back

=head1 OPTIONS

=over 2

=item --append

Append the GFF results to an existing GFF file. The path to the existing
GFF file is indicated by the --outfile option.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

=item --version

Show program version.

=item --man

Show the full program manual. This uses the perldoc command to print the 
POD documentation for the program.

=item -v,--verbose

Run the program in verbose mode.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: Can not open input file

The path to the input file you provided by the --infile option
may be incorrect, or you do not have permission to read the file
that exists at that location.

=item ERROR: Can not open output file

The directory that you specified at the --outfile option may not
exist, or you may not have permission to write files in that 
directory.

=back

=head1 CONFIGURATION AND ENVIRONMENT

The cnv_ltrseq2gff.pl program does not required an external configuration
file or make use of variables defined in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * LTR_seq Program

This program requires output from LTR_seq. LTR_seq is available upon email
request from the author, Ananth Kalyanaraman: 
http://www.eecs.wsu.edu/~ananth/contact.htm.

Also see the original publication for the LTR_par program :

I<A. Kalyanaraman, S. Aluru. 2006. Efficient algorithms and software for 
detection of full-length LTR retrotransposons. Journal of 
Bioinformatics and Computational Biology (JBCB), 4(2):197-216.>

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

=item * Single Record FASTA Files

This program is currently limited to analyzing results from FASTA files
that contain a single record.

=back 

=head1 SEE ALSO

The batch_blast.pl program is part of the DAWG-PAWS package of genome
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

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/03/2007

UPDATED: 03/30/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/03/2007
# - Main body of program written
#
# 12/13/2007
# - Moved POD documentation to the end of the program
# - Updated POD documentation
# - Added print_help subfunction that extracts help
#   and usage messages from the POD documentation
# 03/30/2009
# - Added param
# - Added program
# - Added STDIN support
# - Added STOUT support
# - Added proper feature names for LTR_retrotransposon,
#   three_prime_LTR, and five_prime_LTR
# - Fixed error where vals starting at zero were 
#   set to null
