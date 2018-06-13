#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_merge.pl - Merge all fasta files in a directory     |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 08/20/2006                                       |
# UPDATED: 04/06/2009                                       |
#                                                           |
# DESCRIPTION:                                              | 
#  Convert TE Nest output to gff file format.               |
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
#
my $indir;
my $outfile;

# BOOLEANS                                                                                  
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode

#-----------------------------+
# GET OPTIONS FROM THE        |
# COMMAND LINE                |
#-----------------------------+
my $ok = GetOptions(
                    # Required Arguments
                    "i|indir=s"    => \$indir,
                    "o|outfile=s"  => \$outfile,
                    # Booleans
                    "verbose"      => \$verbose,
                    "test"         => \$test,
                    "usage"        => \$show_usage,
                    "version"      => \$show_version,
                    "man"          => \$show_man,
                    "h|help"       => \$show_help,
                    "q|quiet"      => \$quiet,);

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
    print "\ncnv_tenest2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

&cat_fasta ($indir, $outfile);

exit;

sub cat_fasta { #  Begin of CatBlast Subfunction
    

    my ($indir, $fasta_file_path) = @_;
    
    my $file_num = 0;
    
    # First warn the user if the output file already exists
    my $ArrayRefI = 0;                # Start the arrray ref at 0
    
    # Set the current working directory to the directory that was selected
    #chdir $Directory;
    
    #-----------------------------+
    # Get the FASTA files from the|
    # directory provided by the   |
    # var $indir                  |
    #-----------------------------+
    opendir( DIR, $indir ) ||
	die "Can't open directory:\n$indir";
    my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
    closedir( DIR );
    
    my $total_files = @fasta_files;
    
    # APPEND ANY EXISTING FILE OR WRITE TO STDOUT
    if ($fasta_file_path) {
	open(FASTAOUT, ">> $fasta_file_path") ||
	    die "Can not open fasta file:\n$fasta_file_path\n";
    }
    else {
	open(FASTAOUT, ">&STDOUT") ||
		die "Can not print to STDOUT\n";
    }

    for my $ind_file (@fasta_files) {
	
	$file_num++;
	my $in_file_path = $indir.$ind_file;

	print STDERR "Processing ".$file_num." of ".$total_files."\n" 
	    if $verbose;

	open(FASTAFILE, $in_file_path) || 
	    die "Can not open input file:\n$in_file_path";
	
	while (<FASTAFILE>) {  # For every line in the FASTA file
	    print FASTAOUT ;   # Print the output to FASTAOUT
	}
	close(FASTAFILE);
	
    } # End of foreach @BlastList loop
    
    close(FASTAOUT);
    
} # End of cat_fasta Subfunction


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

__END__

=head1 NAME

fasta_merge.pl -  Merge all fasta files in a directory

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

   fasta_merge.pl -i in_dir/ -o out_file.fasta

=head2 Required Arguments

    -i,--indir    # Path to the directory containing the fasta files
    -o,--outfile  # Path to the fasta format outfile
                  # Writes output to standard output otherwise

=head1 DESCRIPTION

Given a directory, merge all of the fasta files in
the directory into a single fasta file.

=head1 ARGUMENTS

=over 2

=item -i,--indir

Path to the input directory containing the fasta files.

=item -o,--outfile

Path to the fasta file that will be produced. If an output file is
not specified at the command line, the output will be writtent to STDOUT.

=back

=head1 OPTIONS

=over

=item --verbose

Run the program in verbose mode, this will write the maximum amount of 
output to STDERR.

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

The error message that can be generated will be listed here.

=head1 DEPENDENCIES

This program does not depend on external software or perl modules
outside of those included in a typcial installation of perl.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Please report limitations to using this software

If you find a limitation that makes it difficult to use this program, please
file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

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

GNU LESSER GENERAL PUBLIC LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/20/2006

UPDATED: 04/06/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 08/20/2006
#  - Modified program from a BLAST concatenation program I 
#    wrote a while ago.
# 04/04/2007
#  - Changed from hard coded variables to command line
# 04/30/2007
#  - Added POD style documentation
# 04/06/2009
#  - Updated POD docs
#  - Updated command line options
#  - Updated to print_help subfunction that extracts help
#    statements from POD docs
