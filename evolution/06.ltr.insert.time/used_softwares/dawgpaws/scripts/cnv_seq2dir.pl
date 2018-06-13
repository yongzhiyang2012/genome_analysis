#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_seq2dir.pl - Split sequence file to dir of seq files  | 
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/22/2004                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Takes an input file in any valid bioperl format and      |
#  converts it to a directory containing one fasta file for |
#  each sequence. Files are placed in a subdirectory under  |
#  the current working directory that is given the same     |
#  name as the outfile prefix.                              |
#                                                           |
# VERSION: Release 1.0                                            |
#                                                           |
# Usage:                                                    |
#                                                           |
#  cnv_seq2dir.pl -i infile.fasta -f fasta -o outdir/       |
#                                                           |
#-----------------------------------------------------------+

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Bio::SeqIO;               # Read and write seq files
use Cwd;                      # Get the current working directory
use Getopt::Long;
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path
use Cwd;                       # Get the current working directory
use File::Copy;                # Copy files

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# VARS With Default Values
my $infileformat = "fasta";
my $infile;
my $outdir;
my $root_name;             # Can also set a root name
my $num_pad;               # The number of digits to pad the outfile names

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"   => \$infile,
		    "f|format=s"   => \$infileformat,
                    "o|outdir=s"   => \$outdir,
		    # ADDITIONAL OPTIONS
		    "p|pad=i"      => \$num_pad,
		    "r|rootname=s" => \$root_name,
		    "q|quiet"      => \$quiet,
		    "verbose"      => \$verbose,
		    # ADDITIONAL INFORMATION
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
    print "\ncnv_seq2dir.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

# Get the current working directory
unless ($outdir) {
    my $cwd = cwd();
    $outdir = $cwd."ind_seq_out";
}

#-----------------------------+
# INFILE                      |
#-----------------------------+
my $seq_in;
if ($infile) {
    $seq_in = Bio::SeqIO->new('-file'   => "<$infile",
			      '-format' => $infileformat);
}
else {
    $seq_in = Bio::SeqIO->new('-fh'     => \*STDIN,
			      '-format' => $infileformat);
}

#-----------------------------+
# CREATE OUTPUT DIRECTORY     |
#-----------------------------+
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

unless (-e $outdir) {
    mkdir $outdir;
}

$seq_num = 1;

#-----------------------------+
# WRITE EACH SEQ TO OUTDIR    | 
#-----------------------------+
while (my $inseq = $seq_in->next_seq) {

    my $seq_id = $inseq->primary_id();
    
    # Pad number with zeros so that the total length of the 
    # string is 7 characaters (ie. 0000012)

    # GET NUM
#    $num = sprintf("%7d", $seq_num);
    my $num;
    if ($num_pad) {
	$num = sprintf("%".$num_pad."d", $seq_num);
	$num=~ tr/ /0/;
    }
    else {
	$num = $seq_num;
    }

    my $out_file_path;
    if ($root_name) {
	$out_file_path = $outdir.$root_name."_".$num.".fasta";
    }
    else {
	$out_file_path = $outdir.$seq_id.".fasta";
    }

    my $seq_out = Bio::SeqIO->new('-file' => ">$out_file_path",
				  '-format' => 'fasta');

    print STDERR "Processing: $seq_id\n" if $verbose;
    
    # Write the output file
    $seq_out->write_seq($inseq);  # Write the individual fasta file  
    $seq_num++;                   # Increment SeqNumber
}

print STDERR "The files have all been placed in the directory:\n$outdir\n" 
    if $verbose;

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

=head1 NAME

cnv_seq2dir.pl - Split sequence file to dir of seq files

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_seq2dir.pl -i infile.fasta -f fasta -o dir [-r name]

=head2 Required Arguments

    --infile        # Path to the input file
                    # Expects STDIN if not specified
    --format        # Format of the sequence file
                    # (fasta|genbank)
    --outdir        # Path to the output directory
    --rootname      # Base name for output files
                    # Uses seq id if not specified

=head1 DESCRIPTION

Takes an input file in any valid bioperl format and
converts it to a directory containing one fasta file for
each sequence. Files are placed in a subdirectory under
the current working directory that is given the same
name as the outfile prefix.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

The input file containing multiple sequence records to split into single
sequence files.

=item -o, --outdir

The directory to place all of files with the single sequence records in. If
this option is not specified. The sequences will be placed in a subdirectory
of the current working directory named ind_seq_out.

=back

=head1 OPTIONS

=over 2

=item -r,--rootname

A root name to serve as the individual names of the output files. These will
have a running number suffix appended to the outfile name. If not specified,
the primary id on the input record will be used to name the output files.

=item -p, --pad

An integer representing the number of digits to pad file names with. 
This will be used to append to
file names when the --rootname option is used. For example, this will
transfrom root_1 to root_0001. The default is not to pad the numbers.

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

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=back

=head1 DIAGNOSTICS

=over 2

=item ERROR: Could not create root name dir

The program could not create the directory needed to hold the output for
the sequence file. It is possible that you not have the proper access to
create a directory in the path that you specified. It is also possible
that the directory that you are trying to create this subdirectory in 
does not exist.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuration file or options set
in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

This program does not required any additional software.

=head2 Required Perl Modules

=over

=item * Bio::SeqIO

This module is required to read the sequence files in and to
write the sequence files out.

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

=item * No known limitaitons.

If you discover a limitation to your use of this software please
email the author or file a bug report at
the DAWG-PAWS sourceforge website:
http://sourceforge.net/tracker/?group_id=204962

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

STARTED: 07/22/2004

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 07/11/2007
#  - Updated code for jperl
# 01/26/2008
#  - Changed default outfile name to the header name for
#    the sequence file from the fasta record
# 01/29/2009
#  - Added command line switches
#  - Added print_help subfunction to get help from POD code
#  - Added POD code
#  - Added default input from STDIN when not specified
#  - Added option for rootname of padded output files
#  - Added option for the number of digits to pad with
