#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# clust_write_shell.pl - Write shell scripts for r cluster  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/26/2007                                       |
# UPDATED: 04/06/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given information from the command line, write the shell |
#  scripts required to run jobs on cluster environments     |
#  using the LSF queuing system.                            |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
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

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode
my $submit_job = 0;            # Submit the shell script that is created
# PACKAGE LEVEL SCOPE
my $num_dir;                   # Number of dirs to split into
my $file_to_move;              # Path to the file to mask
my $file_new_loc;              # New location of the file to moave
#my $hard_mask_out;             # Path to the hardmasked file output
my $logfile;                   # Path to a logfile to log error info
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $msg;                       # Message printed to the log file
my $search_name;               # Name searched for in grep command
my $bac_out_dir;               # Dir for each sequnce being masked
my $base_name;                 # Base name to be used for output etc
my $prog_name;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required
		    "p|program=s"   => \$prog_name,
                    "o|outdir=s"    => \$outdir,
		    "n|num-dir=s"   => \$num_dir,
		    "b|base-name=s" => \$base_name,
		    # Optional strings
		    "logfile=s"     => \$logfile,
		    # Booleans
		    "verbose"       => \$verbose,
		    "test"          => \$test,
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,
		    "q|quiet"       => \$quiet,);

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
    print "\nclust_write_shell.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
# are not present
if ( (!$outdir) || (!$prog_name) || (!$num_dir) || (!$base_name) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: The base name of the directories was not".
	" specified at the command line" if (!$base_name);
    print STDERR "ERROR: The number of directories was not specified".
	" at the command line" if (!$num_dir);
    print STDERR "ERROR: A program name was not specified at the".
	" command line\n" if (!$prog_name);
    print STDERR "ERROR: An output directory was not specified at the".
	" command line\n" if (!$outdir);
    print_help ("usage", $0 );
}


#-----------------------------+
# OPEN THE LOG FILE           |
#-----------------------------+
if ($logfile) {
    # Open file for appending
    open ( LOG, ">>$logfile" ) ||
	die "Can not open logfile:\n$logfile\n";
    my $time_now = time;
    print LOG "==================================\n";
    print LOG "  batch_hardmask.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
}


#-----------------------------+
# CHECK FOR SLASH IN DIR      |
# VARIABLES                   |
#-----------------------------+
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# WRITE SHELL SCRIPTS         |
#-----------------------------+

my $max_num = $num_dir + 1;
for (my $i=1; $i<$max_num; $i++) {

    # open the shell script
    my $shell_path = $outdir.$base_name."_shell".$i.".sh";
    
    print "Writing to $shell_path\n";

    if ($prog_name =~ "batch_blast") {
	open (SHOUT, ">".$shell_path);
	print SHOUT "batch_blast.pl".
	    " -i /scratch/jestill/wheat_in/$base_name$i/".
	    " -o /scratch/jestill/wheat_out/".
	    " -c /home/jlblab/jestill/scripts/dawg-paws/batch_blast_full.jcfg".
	    " -d /db/jlblab/paws/".
	    " --logfile /home/jlblab/jestill/$base_name$i.log";
	close SHOUT;

    }
    elsif ($prog_name =~ "batch_repmask") {
	open (SHOUT, ">".$shell_path);
	close SHOUT;

    }
    else {
	print STDERR "The program name is not recognized: $prog_name\n";
	exit;
    }
    
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

1;
__END__

=head1 NAME

clust_write_shell.pl - Write shell scripts for the r cluster

=head1 VERSION

This documentation refers to fasta_dirsplit version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    clust_write_shell.pl -p program -o outdir -b job -n 16

=head2 Required Variables

    -p   # program to write the shell script for
    -o   # output dirctory to 
    -b   # base name for the shell scripts to write
    -n   # number of shell scripts to write
         # this will also be used to reference dirs

=head1 DESCRIPTION

Given information from the command line, write the shell
scripts required to run jobs on cluster environments
using the LSF queuing system.

=head1 REQUIRED ARGUMENTS

=over 2

=item -p,--program

Program to write shell scripts for

=item -o,--outdir

Path of the directory to place the program output.

=item -n,--num-dir

Number of dirs to split the parent dir into.

=item -b, --base-name

Name to use a base name for creating the subdirectory

=item 

=back

=head1 OPTIONS

=over 2

=item --logfile

Path to a file that will be used to log program status.
If the file already exists, additional information will be concatenated
to the existing file.

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

=item --test

Run the program without doing the system commands.

=item --verbose

Run the program in verbose mode.

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

This program does not currently rely on an external configuration file
or variables set in the user's environment.

=head1 DEPENDENCIES

This program does not require external software, but the following
Perl modules are required.

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

=item * Limited to cluster scripting for the LSF queueing system.

=back

=head1 SEE ALSO

The clust_write_shell.pl program is part of the DAWG-PAWS package of genome
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

STARTED: 07/26/2007

UPDATED: 04/06/2009

VERSION: Release 1.0

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/24/2007
# - Program started
# - Basic outline with help, usage and man working
# - Reads dir of fasta files and moves names to array
#
# 12/17/2007
# - Moved POD documentation to the end
# - Added SVN tracking of $Rev
# - Added print_help subfunction that extracts help
#   from the POD documentation
# - Changed program version to SVN $Rev ID
# - Changed $ver to $VERSION
# - Added more informative ERROR messages when the 
#   required arguments are not present
