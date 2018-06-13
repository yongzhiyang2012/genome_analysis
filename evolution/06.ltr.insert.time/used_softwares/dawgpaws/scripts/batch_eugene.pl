#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_eugene.pl - Run EuGene gene prediction in batch mode|
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 11/07/2008                                       |
# UPDATED: 01/20/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the EuGene gene prediction progrm in batch mode      |
#  on a directory of fasta files.                           |
#                                                           |
# USAGE:                                                    |
#  batch_eugene.pl -i indir/ -o outdir/                     |
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
use File::Copy;                # Copy file from one location to another

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;
my $organism = "rice";           # The organism

my $eugene_bin = "eugene";       # Path to the eugene binary
my $get_sites =  "egn_getsites4eugene.pl";
my $param_file = "hex_eugene.par";
# The path to the  egn_getsites4eugene.pl program
#my $get_sites = "/home/jestill/projects/wheat_annotation/20081104_ann/".
#    "eug/egn_getsites4eugene.pl";
# The following links to the parameters file to use with eugene
#my $param_file = "/home/jestill/projects/wheat_annotation/20081104_ann/".
#    "hex_eugene.par";

# BOOLEANS
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_test = 0;                  # Run the program in test mode

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"       => \$indir,
                    "o|outdir=s"      => \$outdir,
		    # ADDITIONAL OPTIONS
		    "eugene-binary=s" => \$eugene_bin,
		    "get-sites=s"     => \$get_sites,
		    "organism=s"      => \$organism,
		    "p|param=s"       => \$param_file,
		    "q|quiet"         => \$quiet,
		    "verbose"         => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"           => \$show_usage,
		    "test"            => \$do_test,
		    "version"         => \$show_version,
		    "man"             => \$show_man,
		    "h|help"          => \$show_help,);

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
    print "\nbatch_eugene.pl:\n".
	"Version: $VERSION\n\n";
    exit;
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
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );

my $count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print STDERR "NUMBER OF FILES TO PROCESS: $count_files\n" if $verbose;

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

for my $ind_file (@fasta_files) {

    my $name_root;
    
    #-----------------------------+
    # Get the root name of the    |
    # file to mask                |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.hard\.fasta$/) {
	# file ends in .hard.fasta
	$name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	# file ends in .fasta
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	# file ends in .fa
	$name_root = "$1";
    } 
    else {
	$name_root = $ind_file;
    }

   
    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n";
    }
    
    #-----------------------------+
    # CREATE THE EUGENE DIR       |
    #-----------------------------+
    my $eugene_dir = $name_root_dir."eugene/";
    unless (-e $eugene_dir) {
	mkdir $eugene_dir ||
	    die "Could not create eugene out dir:\n$eugene_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # This will hold the gff file modified from the
    # original gff fine
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create GFF out dir:\n$gff_dir\n";
    }

    #-----------------------------+
    # COPY THE FASTA FILE TO THE  |
    # EUGENE DIRECTORY            |
    #-----------------------------+
    # Use the File:Copy
    my $src_fasta = $indir.$ind_file;
    my $eug_fasta = $eugene_dir.$name_root.".fasta";
#    my $eug_fasta = $eugene_dir.$ind_file;
    copy ($src_fasta, $eug_fasta) ||
	die "Can not copy\n$src_fasta TO \n$eug_fasta\n";
    
    
    #-----------------------------+
    # RUN THE GET SITES PROGRAM   |
    #-----------------------------+
    # This gets splices and starts using online resources
    my $get_sites_cmd = "$get_sites $eug_fasta";
    
    print STDERR "$get_sites_cmd \n" if $verbose;

    system ($get_sites_cmd);
    
    #-----------------------------+
    # RUN EUGENE                  |
    #-----------------------------+
    # 
#    my $eugene_cmd = $eugene_bin." -p g -DEugene.organism=rice ".
#    my $eugene_cmd = $eugene_bin." -p g -DEugene.organism=rice ".
#	$eug_fasta;

    # -p g will just create gff format output
    # -g gh will do gff, and html with png images
    # gff is non standard gff

    # Use the parameters file
    my $eugene_cmd = $eugene_bin." -p g -A $param_file ".
	" -O $eugene_dir".
	" $eug_fasta";


    print STDERR "$eugene_cmd" if $verbose;
    
    system ($eugene_cmd);

    #-----------------------------+
    # CONVERT EUGENE GFF TO       |
    # APOLLO GFF                  |
    #-----------------------------+
    
#    my $eug_gff = $eug_fasta.".gff";               # The EuGene gff file
    my $eug_gff = $eugene_dir.$name_root.".gff";    # The EuGene gff file
    my $ap_gff = $gff_dir.$name_root."_eugene.gff"; # Apollo gff file
    open (EUGIN, "<$eug_gff") ||
	die "Can not open EuGene gff file $eug_gff";
    open (GFFOUT, ">$ap_gff") ||
	die "Can not open gff output\n";
    
    while (<EUGIN>) {

	# Can choose to ignore UTRs if necessary
	my @gff_parts = split;
	my $num_gff_parts = @gff_parts;
	
	# The middle part of this is gene model number
	# Name as HEX20.9.0
	# Prefix is first five characters of the fasta name
	# Middle is gene model number
	# Last part is 0 for UTRS, exon order number otherwise
	my @name_parts = split ( /\./, $gff_parts[0] );
	my $model_num = $name_parts[1];

	print "MODEL NUMBER: $model_num\n" if $verbose;

	# Positions are listed as
	# UTR5
	# UTR3
	# E.Init
	# E.Term
	# E.Intr
	# E.Sngl
	# How to translate this to Apollo
	# For now ignore 5 and 3 UTR label all as exon
	# otherwise add 5 and three to either end
	unless ($gff_parts[2] =~ "UTR5" ||
		$gff_parts[2] =~ "UTR3") {
	    
	    my $start;
	    my $end;

	    if ($gff_parts[3] < $gff_parts[4]) {
		$start = $gff_parts[3];
		$end = $gff_parts[4];
	    } else {
		$end = $gff_parts[4];
		$start = $gff_parts[3];
	    }

	    # Need to replace the 
	    # Addiing simp to eugene model name
	    # to indicate that external information
	    # was not used besides the rice matrix
	    print GFFOUT "$name_root\t".   # 1
		"eugene\t".                # 2
		"exon\t".                  # 3
		"$start\t".                # 4
		"$end\t".                  # 5
		$gff_parts[5]."\t".        # 6
		$gff_parts[6]."\t".        # 7
		$gff_parts[7]."\t".        # 8
		"eugene_simp_".$model_num.      # 9
		"\n";
	    
	    print "\t$gff_parts[2]\n" if $verbose;
	}

    }
    
    close (EUGIN);
    close (GFFOUT);
    
}


exit 0;

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

batch_eugene.pl - Run EuGene gene prediction in batch mode

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    batch_eugene.pl -i InDir -o OutDir

=head2 Required Arguments

    --indir         # Path to the input directory
    --outdir        # Path to the output directory

=head1 DESCRIPTION

This script runs the EuGene gene annotation program in batch mode on a
directory of FASTA files.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=back

=head1 OPTIONS

=over 2

=item --eugene-binary

Path to the binary for the eugene program. This file should be named
eugene. If this path is not supplied at the command line, the program
assumes that the directory containing this program is listed in the
user's PATH environment variable.

=item --get-sites

Path to the program egn_getsites4eugene.pl. If this path is not supplied at 
the command line, the program
assumes that the directory containing this program is listed in the
user's PATH environment variable.

=item -p,--param

Path to the parameter file to use. If this path is not supplied at the 
command line, the program will use the hex_eugene.par file and will assume
that the directory cotaining this file is listed in the user's
PATH evnironment variable.

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

=item --test

Run the program without doing the system commands. This will
test for the existence of input files.

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

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuartion file. Instead a
eugene specific parameter file must be set with the --param option.

=head2 Environment

This program does not make use of variables in the user environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * EuGene

The batch_eugene.pl program requires that you have a local installation of
the EuGene gene annotation program. This program is available from:
http://www.inra.fr/internet/Departements/MIA/T/EuGene/

=back

=head2 Required Perl Modules

=over

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

=item * Limited testing on version of the EuGene program

This program has been tested and known to work with 
EuGene rel 3.3.

=back

=head1 SEE ALSO

The program is part of the DAWG-PAWS package of genome
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

STARTED: 11/07/2008

UPDATED: 01/20/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 11/07/2008
# -Main body of program designed and written
#
# 01/20/2009
# - Updated POD documentation
