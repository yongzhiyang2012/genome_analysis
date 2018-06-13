#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dir_merge.pl - Merge directories                          |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 07/27/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Merge a two directories containing subdirs with the same |
#  name into a single dir. One dir can serve as the template|
#  or an entirely new directory set could be created.       |
#                                                           |
# VERSION: Release 1.0                                                   |
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
use File::Copy;
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
my $indir_list;                # List of dirs passed from the cmd line
my $outdir;                    # Directory to serve as parent for output
my @indirs;                    # List of directories
my $num_indirs;                # Total number of indirs

# Booleans
my $do_overwrite = 0;          # Overwrite existing files
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $err = 0;                  # Boolean to indicate program encountered error
my $verbose = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions("i|indirs=s"  => \$indir_list,
                    "o|outdir=s"  => \$outdir,
		    # Booleans
		    "overwrite"   => \$do_overwrite,
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
		    "verbose"     => \$verbose,
		    "q|quiet"     => \$quiet,);


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
    print "\ndir_merge.pl:\nVersion: $VERSION\n\n";
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
# Exit if indir_list and outdir not specified at cmdline
if ( (!$indir_list) || (!$outdir) ) {
    print "\a";
    print STDERR "ERROR: An input list of dirs must be ".
	" specified at the command line.\n" if (!$indir_list);
    print STDERR "ERROR: An base output directory must be specified".
	" at the command line.\n" if (!$outdir);
    print_help ("usage", $0 );
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

#-----------------------------+
# GET THE INPUT DIRS TO MERGE |
#-----------------------------+
@indirs = split(/\,/, $indir_list);

# List dirs to split
$num_indirs = @indirs;

# Got rid of this
# Throw error if multiple dirs not indicated
#if ($num_indirs < 2) {
#    print "\a";
#    print "ERROR: Multiple dirs shold be indicated by --indir\n";
#    print ""; 
#}

print "\nDirs to merge ($num_indirs):\n";
for (my $i=0; $i < $num_indirs; $i++) {

    # ADD SLASH TO END OF DIR PATH
    unless ($indirs[$i] =~ /\/$/ ) {
	$indirs[$i] = $indirs[$i]."/";
    }

    # CHECK FOR EXISTENCE OF INPUT DIRECTORIES
    if (-e $indirs[$i]) {
	print "\t-".$indirs[$i]."\n";
    }
    else {
	$err = 1;
	print "\a";
	print "ERROR: Directory does not exist\n\t".$indirs[$i]."\n";
    }
}

#-----------------------------+
# GET THE OUTPUT DIR          |
#-----------------------------+
unless ($outdir =~ /\/$/ ) {
    $outdir = $outdir."/";
}
print "New dir: $outdir\n\n";

# CREATE OUTDIR IF IT DOES NOT EXIST
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}


#-----------------------------+
# FOR EACH DIR TO COPY VIA THE
# merge SUBFUNCTION
#-----------------------------+
for my $ind_dir (@indirs) {
    merge ($ind_dir, $outdir);
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

sub merge {
# SUBFUCNTION SRC:
# http://software.hixie.ch/utilities/unix/merge/merge.pl

    my($a, $b) = @_;
    #   SOURCE: not -e    -f    -d   empty -d else
    # TARGET:
    #  not -e    err      mv    mv    mv      err
    #
    #      -f    err      err   err   err     err
    #
    #      -d    err      err   loop  rmdir   err
    #
    #    else    err      err   err   err     err
    print "merge: merging $a and $b\n";
    if (not -e $a) {
	#-----------------------------+
	# $a DOES NOT EXIST           |
	#-----------------------------+
        print STDERR "merge:$a: doesn't exist\n";

    } 
    elsif (not (-f $a or -d $a)) {
	#-----------------------------+
	# $a IS NOT A FILE OR DIR     |
	#-----------------------------+
        print STDERR "merge:$a: not a normal file\n";

    } 
    elsif (not -e $b) {
	#-----------------------------+
	# $b DOES NOT EXIST           |
	#-----------------------------+
        print "merge: moving $a to $b\n" if $verbose;
        rename($a, $b) || 
	    print STDERR "merge:$a: could not rename to $b, $!\n";;
	# JAMIE MOD BELOW -- curently not working
        #print "dir_merge: moving $a to $b\n" if $verbose;
	#copy ($a, $b) ||
	#    print STDERR "merge:$a: could not rename to $b, $!\n";
    } 
    elsif (-d $b) {
	#-----------------------------+
	# $b IS A DIRECTORY           |
	#-----------------------------+
        if (-d $a) {
            my @entries = getdir($a);
            if (@entries) {
                # not empty
                # recurse through it to give us a chance to make it empty
                print "merge: going through contents of $a\n" if $verbose;
                foreach my $entry (@entries) {
                    my $c = "$a/$entry";
                    $c =~ s|//|/|gos;
                    my $d = "$b/$entry";
                    $d =~ s|//|/|gos;
                    &merge($c, $d);
                }
            }
            # empty now?
            @entries = getdir($a);
            if (not @entries) {
                print "merge: deleting empty directory $a\n" if $verbose;
                rmdir($a) or print STDERR "merge:$a: could not delete ".
		    "directory, $!\n";
            } 
	    else {
                print STDERR "merge:$a: could not delete directory,".
		    " directory is not empty\n";
            }
        } 
	else {
            print STDERR "merge:$a: conflicts with directory $b\n";
        }
    } 
    else {
	#-----------------------------+
	# THE FILE ALREADY EXIST WITH |
	# A FILE IN B                 |
	#-----------------------------+
	print  "merge:$a: conflicts with non-directory $b\n";
	
	if ($do_overwrite) {
	    # DELETE EXISTING FILE AND REPLACE WITH NEW FILE
	    #print STDERR "File $b exists and will be overwritten\n";
	    #my $rmcmd = "rm $b";
	    #system ($rmcmd);
	    #rename($a, $b) || 
	    #	print STDERR "merge:$a: could not rename to $b, $!\n";
	    move($a, $b) || 
		print STDERR "\n\nmerge:$a: could not rename to $b, $!\n";
	}
	else {
	    print STDERR "merge:$a: conflicts with non-directory $b\n";
	}
	
    }
}


sub getdir {
# SUBFUCNTION SRC:
# http://software.hixie.ch/utilities/unix/merge/merge.pl
    my($a) = @_;
    local *DIR;

    unless (opendir(DIR, $a)) {
        print STDERR "dir_merge:$a: can't open directory\n";
        return;
    }

    my @entries;
    while (my $entry = readdir(DIR)) {
        if ($entry !~ m/^\.\.?$/o) {
            push(@entries, $entry);
        }
    }
    
    closedir(DIR) || 
	print STDERR "dir_merge:$a: could not close directory,$!\n";
    
    return @entries;
}

1;
__END__

# Deprecated print_help subfunction
sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"dir_merge.pl -i DirsToMerge -o OutDir";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Path to the directories to merge\n".
	"                 # Dirs separated by comma\n".
	"  --outdir       # Path to the output directory\n".
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

dir_merge.pl - Merge directories

=head1 VERSION

This documentation refers to program version 0.1

=head1 SYNOPSIS

=head2 Usage

  dir_merge.pl -i 'dir01,dir02' -o new_dir

=head2 Required Arguments

    -i, --indirs   # List of directories to merge
    -o, --outdir   # Path to the base output directory

=head1 DESCRIPTION

Given a set of directories, this will merge them into a single output
directory. The output directory could be a new location, or merge into
the parent directory.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indirs

The directories to merge. Multiple directories are separated by commas.

=item -o,--outdir

The directory the merged dirs will be moved to.

=back

=head1 OPTIONS

=over 2

=item --overwrite

Overwrite any existing file. This option should scare the living
hell out of you. Only use it if you are certain you will not be overwriting
anything that you do not want to lose. Use this option at your own risk.

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

This program does not currently require an external configuration
file or make use of variables set in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

The dir_merge.pl program does not rely on external software.

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

=over 2

=item * Limited testing

This program has only been tested using the Linux operating system
I have no idea what this would do in windows. If you try this in widows 
and it works, drop me an email E<lt>JamesEstill at gmail.comE<gt>.

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

GNU GENERAL PUBLIC LICENSE, VERSION 3

http://www.gnu.org/licenses/gpl.html   

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 07/27/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 07/27/2007
# - Program started
# - Added getdir from external source
# 12/11/2007
# - Added SVN tracking of Id and Rev
# - Moved POD documentation to the end of the program
# - Update POD documentation
# - Added print_help subfunction that extracts help and
#   usage information from the POD documentation
#
# 11/06/2008
# - Added some info to the POD documentation
