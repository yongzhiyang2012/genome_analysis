#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# seq_oligocount.pl - Count oligos of size k
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 10/11/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Short Program Description                                |
#                                                           |
# USAGE:                                                    |
#  seq_oligocount.pl 
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
use Bio::SeqIO;                # Seq IO used to to work with the input file
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
my $infile;                    # Path for the input query sequence
my $outdir;                    # Path for the output
my $index;                     # Path to the sequence index file
my $kmer_len = 20;             # Oligomer length, default is 20
my $win_len = 50;              # Window length
my $seq_name;                  # Sequence name used in the output

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff_kmer = 0;           # Produce GFF formatted output for kmers
my $do_gff_win = 0;            # Produce GFF formatted output for win summary
my $do_wig_kmer = 0;           # Produce wiggle output for kmers
my $do_wig_win = 1;            # Produce wiggle output for win summary

my $param;                     # The parameter set
my $program="vmatch";          # The source program the result is derived from

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
		    "d|db=s"      => \$index,
                    "o|outdir=s"  => \$outdir,
		    "n|s|seqname|name=s" => \$seq_name,
		    # ADDITIONAL OPTIONS
		    "param=s"     => \$param,
		    "program-s"   => \$program,
		    "k|kmer=s"    => \$kmer_len,
		    "l|len=s"     => \$win_len,
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
if ( ($show_usage) ) {
#    print_help ("usage", File::Spec->rel2abs($0) );
    print_help ("usage", $0 );
}

if ( ($show_help) || (!$ok) ) {
#    print_help ("help",  File::Spec->rel2abs($0) );
    print_help ("help",  $0 );
}

if ($show_version) {
    print "\nseq_oligocount.pl:\nVersion: $VERSION\n\n";
    exit;
}

if ($show_man) {
    # User perldoc to generate the man documentation.
    system("perldoc $0");
    exit($ok ? 0 : 2);
}

# Throw error if required options not present
if ( (!$infile) || (!$outdir) || (!$index) ) {

    print "\a";
    print "\n";
    print "ERROR: Input file path required\n" if (!$infile);
    print "ERROR: Output directory required\n" if (!$outdir);
    print "ERROR: Index file path required\n" if (!$index);
    print_help("");

}


#-----------------------------+
# CREATE OUTDIR IF NOT PRESENT|
#-----------------------------+
mkdir $outdir, 0777 unless (-e $outdir);

#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+ 

seq_kmer_count ($infile, $outdir, $index, $kmer_len, $win_len, $seq_name);

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

sub seq_kmer_count {

    my ($fasta_in, $outdir, $vmatch_index, $k, $l, $seq_name ) = @_;
    
    # Array of threshold values
    my $thresh = 50;
    my $in_seq_num = 0;

    my $inseq;
    if ($fasta_in) {
	$inseq = Bio::SeqIO->new( -file => "<$fasta_in",
				  -format => 'fasta');
    }
    else {
	$inseq = Bio::SeqIO->new( -fh => \*STDIN,
				  -format => 'fasta');
    }

    # Counts of the number of occurences of each oligo
    my @counts = ();    
    my @vdat;   # temp array to hold the vmatch data split by tab
    my $pos = 0;

    my $start_pos;
    my $end_pos;
    my $i;         # Array index value

    # FILE PATHS
    # These will need to be modified later to more useful names
    my $temp_fasta = $outdir."split_seq_".$k."_mer.fasta";
    my $vmatch_out = $outdir."vmatch_out.txt";
    my $gff_count_out = $outdir."vmatch_out.gff";

    while (my $seq = $inseq->next_seq) {

	$in_seq_num++;
	if ($in_seq_num == 2) {
	    print "\a";
	    die "Input file should be a single sequence record\n";
	}

	# Calculate base cooridate data
	my $seq_len = $seq->length();
	my $max_start = $seq->length() - $k;
	
	# Print some summary data
	print STDERR "\n==============================\n" if $verbose;
	print STDERR "SEQ LEN: $seq_len\n" if $verbose;
	print STDERR "MAX START: $max_start\n" if $verbose;
	print STDERR "==============================\n" if $verbose;
	
	# CREATE FASTA FILE OF ALL K LENGTH OLIGOS
	# IN THE INPUT SEQUENCE
	print STDERR "Creating oligo fasta file\n" if $verbose;
	open (FASTAOUT, ">$temp_fasta") ||
	    die "Can not open temp fasta file:\n $temp_fasta\n";

	for ($i=0; $i<=$max_start; $i++) {

	    $start_pos = $i + 1;
	    $end_pos = $start_pos + $k - 1;

	    my $oligo = $seq->subseq($start_pos, $end_pos);

	    # Set counts array to zero
	    $counts[$i] = 0;
	    
	    print FASTAOUT ">$start_pos\n";
	    print FASTAOUT "$oligo\n";
	    
	}

	close (FASTAOUT);

	# QUERY OLIGO FASTA FILE AGAINST VMATCH INDEX
	my $vmatch_cmd = "vmatch -q $temp_fasta -complete".
	    " $vmatch_index > $vmatch_out";
#	# OUTPUT TO STDOUT FOR DEBUG
#	my $vmatch_cmd = "vmatch -q $temp_fasta -complete".
#	    " $vmatch_index";
	print STDERR "\nVmatch cmd:\n$vmatch_cmd\n" if $verbose;
	system ($vmatch_cmd);

	# PARSE VMATCH OUTPUT FILE
	# increment
	unless (-e $vmatch_out) {
	    print STDERR "Can not file the expected vmatch output file\n";
	    die;
	}

	# PARSE THE VMATCH OUTPUT FILE AND INCREMENT
	# THE COUNTS ARRAY
	open (VMATCH, "<$vmatch_out") ||
	    die "Can not open vmatch output file:\n$vmatch_out\n";

	# Count oligos hits in index file
	print STDERR "\nCounting oligos ...\n" if $verbose;
	while (<VMATCH>) {
	    # ignore comment lines
	    unless (m/^\#/) {
		chomp;
		@vdat = split;
		

		my $len_vdat = @vdat;
		# print the following for debu
		#print "VDAT LEN: $len_vdat\n";
		#print $vdat[5]."\n";
		
		# Get the seq file in the new sub se
		# Counts index starts at zero
		# It may be possible to increment without needing
		# to add one
		#$counts[$vdat[5]]++;
		$counts[$vdat[5]] = $counts[$vdat[5]] + 1;

	    }
	}
	close (VMATCH);


	#///////////////////////////////////
	# NEED SEGMENTATION STEP HERE
	# THIS WILL JOIN INDIVIDUAL PIPS
	# FOR A GIVEN THRESHOLD COVERAGE
	# This could also be done with a 
	# separate script operating on the
	# gff output from this program.
	#\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

	# Output oligo counts to gff file
	
	# This is the pip file, a single pip with depth coverage data
	# for eack oligo

	# Get summary information ... across window

	#-----------------------------+
	# GET AVERAGE OF KMERS THAT   |
	# EACH BASE IS A MEMBER OF    |
	#-----------------------------+
	# Is this necessary
	# Could just smooth for all values that 
	# are in 
#	for ($i=0; $i<=$max_start; $i++) {
#
#
#	}

	# OPEN GFF OUT FILE, STDOUT IF NO FILE PATH USED
	if ($gff_count_out) {
	    open (GFFCOUNT, ">$gff_count_out") ||
		die "Can not open gff out file:\n$gff_count_out\n";
	}	
	else {
	    open (GFFCOUNT, ">&STDOUT") ||
		die "Can not print to STDOUT\n";
	}

	# SET VARIABLES FOR GFF OUT FILE
	if ($param) {
	    $program = $program.":".$param;
	}

	print STDERR "\nCreating gff output files ...\n" if $verbose;

	for ($i=0; $i<=$max_start; $i++) {

	    $start_pos = $i + 1;
	    $end_pos = $start_pos + $k - 1;
	    
	    # The count will be placed in the score position
	    my $seq_str = $seq->subseq($start_pos, $end_pos);
	    print GFFCOUNT "$seq_name\t".  # Ref sequence
		"$program\t".              # Source
 		"oligo_count\t".           # Type
		"$start_pos\t".            # Start
		"$end_pos\t".              # End
		$counts[$i]."\t".          # Score
		".\t".                     # Strand
		".\t".                     # Phase
		"Vmatch ".$k."mer\n";    # Group

# THE FOLLOWING REMOVED 03/26/2009
#	    # PRINT OUT SEQS EXCEEDING THRESHOLD VALUES
#	    if ($counts[$i] > $thresh) {
#		my $thresh_seq = $seq->subseq($start_pos, $end_pos);
#		print "$start_pos\t".$counts[$i]."\t".
#		    "$thresh_seq\n";
#	    }

	}
	
	close (GFFCOUNT);

    } # End of while seq object

    # May want to make a single fasta file of all oligos for the
    # qry fasta_sequence and parse the results from vmatch from that

}

1;
__END__

=head1 NAME

seq_oligocount.pl - Count oligos from an input sequence

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    seq_oligocount.pl -i InFile -o OutDir -db index.fasta
                      -n SeqName -k 20

=head2 Required Arguments

    -i,--infile   # Path to the input fasta file
    -d,--db       # Path to the mkvtree index file
    -o,--outdir   # Path to the base output directory
    -n,--name     # Name to assign to the sequence file
    -k,--kmer     # Oligomer query length

=head1 DESCRIPTION

The seq_oligocount program will take a query sequence, break it
into subsequences of size k and query it against an persistent index
created by the mkvtree program. It produces a GFF output file describing
the number of copies of every oligomer in the query sequence in the 
subject index database.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file.

=item -o,--outdir

Path of the output file.

=item -n,--name

Name to assign to the sequence file

=item -d,--db

Path to the fasta file that was indexed with the mkvtree program.

=back

=head1 OPTIONS

=over 2

=item -k,--kmer

Length of the kmer to index. The default value of this variable is 20.

=item -l, --len

Length of window for summarizing kmer index counts. This will be a non 
overlapping window and can range from 1 to the length of the sequence.

=item -s,--seqname

The name of the sequence being annotated. This is the first column of data
in the GFF output file.

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

Error messages generated by this program and possible solutions are listed
below.

=over 2

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=back

=head1 CONFIGURATION AND ENVIRONMENT

An external configuration file is not required for this program, and
it does not make use of any variables set in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over 2

=item Vmatch

This program requires the Vmatch package of programs.
http://www.vmatch.de . This software is availabe at no cost for
noncommercial academic use. See program web page for details.

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the BLAST results.

=item * Getopt::Long

This module is required to accept options at the command line.

=item * Bio:SeqIO

The Bio:SeqIO module is a component of the BioPerl package

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

=item * Limit to query sequence file

Since this program will generate a large file of your original sequence
broken into oligomers, you are somewhat limited to the size of query
sequence that is feasible to use this program with. 

=item * Output limited to GFF file

The output for this program is currently limited to GFF file of 
points with the copy number of each oligomer designeated.

=back

=head1 SEE ALSO

The seq_oligocount.pl program is part of the DAWG-PAWS package of genome
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

STARTED: 10/11/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 10/11/07 and 10/12/07
# -bulk of program written
# -inital seq_kmer_count subfunction added
#
# 12/12/2007
# - Added SVN tracking of Rev to propset 
# - Updated POD documentation
# - Added the print_help subfunction that extracts
#   usage and help message from the POD documentation
#
# 11/12/2008
# - Removed old print_help subfunction
# - Working to add summary over window
#   and output to wiggle format
# 
# 01/27/2009
# - Updating POD documentation
# - Adding print_help subfunction to extract help from
#   the POD documentation
