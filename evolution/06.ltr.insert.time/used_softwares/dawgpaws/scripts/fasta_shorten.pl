#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# fasta_short.pl - Give fasta file shorter names            |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/17/2007                                       |
# UPDATED: 04/26/2008                                       |
#                                                           |
# DESCRIPTION:                                              |
#   Change headers in a fasta file to give shorter names    |
#   Can shorten header to a set length, or by a delimiting  |
#   character. The output file may also be rename to the    |
#   shorter name, and an lowercase residues can be          |
#   converted to uppercase.                                 |
#                                                           |
# USAGE:                                                    |
#  fasta_shorten -i indir/ -o outdir/                       |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
#
# TO DO:
# Write a general read fasta file from dir for DAWGPAWS

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;                    # Play by the rules
use Getopt::Long;              # Get options from the command line
use File::Copy;                # Copy files in OS independent manner
# The following needed for printing help
use Pod::Select;               # Print subsections of POD documentation
use Pod::Text;                 # Print POD doc as formatted text file
use IO::Scalar;                # For print_help subfunction
use IO::Pipe;                  # Pipe for STDIN, STDOUT for POD docs
use File::Spec;                # Convert a relative path to an abosolute path

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
# Extract program version form the SVN Revision number
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Input dir of fasta files to shorten
my $outdir;                    # Output dir for shortened fasta files
my $new_len = 20;              # New length of the header
my $head_new;                  # New, shorter header
my $test_len;                  # Test length of the header
my $cur_len;                   # Current length of the header
my $outfile;                   # Set scope for outfile
my $delim_char;                # Character used to delimit header
my $delim_pos=1;               # By default the delim position is first postion

# Booleans
my $verbose=0;
my $rename=0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $test = 0;
my $uppercase = 0;             # Convert sequence strings to uppercase

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required options
		    "i|indir=s"      => \$indir,
                    "o|outdir=s"     => \$outdir,
		    # Additional options
		    "l|length=s"     => \$new_len,
		    "d|delim-char=s" => \$delim_char,
		    "p|delim-pos=s"  => \$delim_pos,
		    # Booleans
		    "rename"         => \$rename,
		    "uppercase"      => \$uppercase,
		    "verbose"        => \$verbose,
		    "usage"          => \$show_usage,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,
		    "test"           => \$test,
		    "q|quiet"        => \$quiet,);

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
    print "\nfasta_shorten.pl:\nVersion: $VERSION\n\n";
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
if ( (!$indir) || (!$outdir) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print_help ("usage", $0 );
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

# Convert position as passed at command line to array position
$delim_pos = $delim_pos - 1;

# The test length is the length of the header plus one
#$test_len = $new_len + 1;
$test_len = int($new_len + 1);

#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
# If the indir does not end in a slash then append one
# TO DO: Allow for backslash
unless ($indir =~ /\/$/ ) {
    $indir = $indir."/";
}

unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
my @fasta_files;
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
# Old fasta file aray
# my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
@fasta_files = grep /\.fasta$|\.fa$|\.tfa$|\.fsa$|\.fna$|\.faa$|\.frn$|\.ffn$/, readdir DIR ;

# If none of exptected fasta file extensions found in the input dir, load all 
# files in the input directory and assume that they are fasta files
my $num_fasta_files;
$num_fasta_files = @fasta_files;
if ($num_fasta_files == 0) {
    @fasta_files = grep !/^\.\.?$/, readdir DIR ;
}

closedir( DIR );

# Throw error if no files are found in the directory
$num_fasta_files = @fasta_files;
if ($num_fasta_files == 0) {
    die "ERROR: No fasta files were found in the input directory.\n";
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print STDERR "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# PROCESS EACH INDIVIDUAL     |
# FILE IN THE DIRECTORY       |
#-----------------------------+
for my $ind_file (@fasta_files) {

    # The fasta file name is not changed, it is just moved to a 
    # new directory
    my $infile = $indir.$ind_file;

    open (IN, $infile) ||
	die "Can not open infile:\n$infile\n";
        
    unless ($rename) {
	my $outfile = $outdir.$ind_file;    
	open (OUT, ">".$outfile) ||
	    die "Can not open outfile:\n$outfile\n";
    }

    while (<IN>) {
	
	chomp;

	# If the first characeter is a > we are in the fasta header
	if (/^\>/) {

	    #-----------------------------------------------------------+
	    # FASTA HEADER PROCESSING                                   |
	    #-----------------------------------------------------------+
	    
	    print "Processing $_\n" if $verbose;
	    
	    if ( ( $delim_char)  ) {
		#-----------------------------+
		# DELIMITED BY CHARACTER      |
		#-----------------------------+
		# This allows for a few contingencies if things in the
		# header line are strange. ie. The character is not present

		# Remove the leading >
		s/^\>//;

		my @hp; # Header parts

		#-----------------------------+
		# SPLIT HEADER BY DELIMITER   |
		#-----------------------------+
		if ( $delim_char =~ "pipe") {
		    @hp = split ( /\|/ , $_ );
		}
		elsif ( $delim_char =~ "comma" ) {
		    @hp = split ( /\,/ , $_ );
		}
		elsif ( $delim_char =~ "space" ) {
		    @hp = split ( /s+/ , $_ );
		}
		elsif ( $delim_char =~ "tab" ) {
		    @hp = splilt ( /\t/ , $_ );
		} 
		else {
		    my $die_msg = "An appropriate delimit character".
			"must be passed to use the delimit option.\n";
		    die "$die_msg\n";
		}

		# PRINT OUTPUT
		# This is where the opening of a new outfile for
		# renmae would need to go
		$head_new = $hp[$delim_pos];

		if ($rename) {
		    my $outfile = $outdir.$head_new.".fasta";    

		    if (-e $outfile) {
			print STDERR "\aWARNING: Outfile already exists".
			    "\n\t$outfile\n".
			    "\texisting file will be overwritten.\n";
		    }
		    
		    open (OUT, ">".$outfile) ||
			die "Can not open outfile:\n$outfile\n";
		}
		
		print OUT ">$head_new\n" || 
		    die "Can not print to output file handle\n";

	    }
	    else {

		#-----------------------------+
		# LENGTH THRESHOLD            |
		#-----------------------------+
		# Using a length criteria is the default when a delim
		# character is not defined at the command line
		$cur_len = int(length($_));
		
		print STDERR "\tLEN: $cur_len\n" if $verbose;
		if ( $cur_len > $test_len ) {
		    # If the fasta header is longer then the desired length 
                    # then create a shortened header print to the outfile
		    # We will start at 1 instead of zero to ignore
		    # the > character
		    $head_new = substr ( $_, 1, $new_len );
		    print STDERR "\tNEW: >$head_new\n" if $verbose;

		    # OPEN RENAMED OUTPUT FILE HANDLE
		    if ($rename) {
			my $outfile = $outdir.$head_new.".fasta"; 

			if (-e $outfile) {
			    print STDERR "\aWARNING: Outfile already exists".
				"\n\t$outfile\n".
				"\texisting file will be overwritten.\n";
			}
   
			open (OUT, ">".$outfile) ||
			    die "Can not open outfile:\n$outfile\n";
		    }

		    print OUT ">$head_new\n";

		} 
		else {
		    
		    # If not longer, then use existing string
		    $head_new = substr ( $_, 1);
		    if ($rename) {
			my $outfile = $outdir.$head_new.".fasta";    

			if (-e $outfile) {
			    print STDERR "\aWARNING: Outfile already exists".
				"\n\t$outfile\n".
				"\texisting file will be overwritten.\n";
			}

			open (OUT, ">".$outfile) ||
			    die "Can not open outfile:\n$outfile\n";
		    }
		    
		    # print the header to the outfile unchanged		    
		    print OUT "$_\n";
		    
		} # End of if header is too long

	    } # End of chose between delim and length
	}
	
	else {
	    #-----------------------------+
	    # SEQUENCE STRING PROCESSING  |
	    #-----------------------------+
	    # print the string to the outfile unchanged if this is
	    # not a fasta header line
	    # T
	    #if ($)
	    if ($uppercase) {
		# Convert the sequence string to uppercase
		# I will try this as a bare regexp since 
		# we are workign with $_
		tr/a-z/A-Z/;
		print OUT "$_\n";

	    }
	    else {
		# If not converting the sequence string to uppercase
		# print the output file unchanged.
		print OUT "$_\n";
	    }
	}


    } # End of while IN

    close (IN);
    close (OUT);
 
    # This rename by move could really wreak havoc
    # There are definitely better ways to do this, but this
    # should work for now. This will rewrite any existing file
    # that has the same name at the new location. Therefore
    # if fasta headers are not unique you will loose data.
    # REMOVED THE FOLLOWING 04/06/2009
#    if ( ($rename) ) {
#	my $new_outfile = $outdir.$head_new.".fasta";
#	move ($outfile,$new_outfile) ||
#	    die "Count not move:\n$outfile to\n$new_outfile";
#    }

} # End of for each file in the directory

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

1;
__END__

=head1 NAME

fasta_shorten.pl - Change headers in a fasta file to give shorter names.

=head1 VERSION

This documentation refers to fasta_shorten version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    fasta_shorten.pl -i InDir -o OutDir

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

This program will take all of the fasta files in a input directory 
and will shorten the name in the fasta headers. 
Creating shorter names is often required for programs that have a 
maximum length that they can use for fasta headers.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item -l,--length

New length. The fasta header will be shortened to this length.
Default length is 20.

=item -d,--delim-char

The delimiting character to use. Note that since PERL regular expressions
can not use variables, I must limit this to a list of valid characters.
Valid choices are:

=over

=item pipe [-d pipe]

The pipe character '|' as used by NCBI.

=item comma [-d comma]

A comma ',' is used to delimit the header line.

=item space [-d space]

To delimit by any whitespace type in the word space.

=item tab [-d tab]

To delimit by the tab character, type tab in the command line.

=back

=item -p,--delim-pos

The position in the delimited set to use. Be default this will
be the fist position in the split array. If delim-pos is greater 
then the number of split characters, the first position will be used
and a message sent to STDERR.

=item --upercase

Convert all sequence residues to uppercase.

=item --rename

Rename the output fasta files to the new FASTA header name.

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

=item ERROR: No fasta files were found in the input directory

The input directory does not contain fasta files in the expected format.
This could happen because you gave an incorrect path or because your sequence 
files do not have the expected *.fasta extension in the file name.

=item ERROR: Could not create the output directory

The output directory could not be created at the path you specified. 
This could be do to the fact that the directory that you are trying
to place your base directory in does not exist, or because you do not
have write permission to the directory you want to place your file in.

=item WARNING: Outfile already exits

This will occur when you are choosing to rename the output fasta file, and the
new file name is not unique. The existing file will be overwritten by
the new file. You can avoid this problem by keeping the original
sequence name, by selecting a longer string (-l) to generate unique
new strings, or by choosing a different delimit character to generate
unique names.

This may also occur when you are writing new fasta files to a directory
that already contaings a fasta file with the same name.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not currently make use of configuration files
or settings in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

No external software is currently required to use this program

=head2 Required Perl Modules

=over

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

=item * Renaming outfile limited to single record fasta files

Multiple record fasta files may be reanamed to multiple files or 
throw an error.

=item * No additional known limitations

If you find limitations to your use of this software please email the
author with information regarding your operating system and what
limitations you are experiencing with your use of the software.

=back

=head1 SEE ALSO

The fasta_shorten.pl program is part of the DAWG-PAWS package of genome
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

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 07/17/2007

UPDATED: 04/26/2008

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/17/2007
# - Program started.
#
# 12/10/2007
# - Moved POD docs to end of script
# - Changed program version reference to SVN revision
# - Added new help subfunction that extracts help and usage
#   text from POD documentation
# - Added check for required arguments
#
# 02/18/2008
# - Adding varialbles
#    - delim-char - Delimiting character    
#    - delim-pos  - Char position to select from delimited set
#                   This starts at zero
#    - rename     - Rename the output fasta file to the new
#                   shortened name. Otherwise the original
#                   file name is kept.
# - The delim option has only been fully tested on pipe delim data
# - Added the dependency on the File::Copy module to rename output file
#
# 04/26/2008
# - Added comments to code
# - Minor updates to POD documentation
# - Added code to implement the rename option. 
#    - The new output file will be placed in the output directory
#      and will have the same name as the new header.
#    - Default option is not to rename the file
#    - This option will only be stable with single record fasta
#      files, otherwise it will attempt to open a new output file
#      for each fasta record in the file
# - Adding options for fasta file extensions that are 
#   recognized. These are used to load the @fasta_files array:
#    - Added *.tfa, *.fsa as commonly used extensions
#    - Added NCBI extensions:
#       -fna - whole genomic DNA sequences
#       -faa - protein coding sequences (CDS)
#       -ffn - untranslated nucleotide sequences
#       -frn - nucleotide sequences of RNA related features
# - Adding default mode to use all the files in the directory
#
# 04/06/2009
# -Minor code cleanup
# -Fixed OUT file handle for length delim output with rename
