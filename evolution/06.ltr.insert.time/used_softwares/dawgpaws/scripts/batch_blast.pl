#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_blast.pl - Run blast on a set of fasta files.       |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 07/23/2007                                       |
# UPDATED: 04/23/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Given a directory of softmasked fasta files, this will   |
#  BLAST the files against a standard set of BLAST          | 
#  databases used for wheat genome annotation.              |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+

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
use Bio::SearchIO;             # Parse BLAST output

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
my $do_gff = 1;                # Try to convert the blast results to gff

# PACKAGE LEVEL SCOPE

#Counters
my $i;                         # Number of blast files
my $count_db;                  # Number of databases to process
my $count_files;               # Number of files to process
my $count_proc;                # Number of processes
my $proc_num;                  # Number of the current process
my $file_num = 0;              # Number of the current file

# Files
my $file_config;               # Path to the Batch blast config file
my $file_to_blast;             # Path to the file to BLAST, the qry file
my $file_blast_out;            # Path to BLAST output file
my $logfile;                   # Path to a logfile to log error info
# Path to the NCBI blastall binary
my $blast_path = $ENV{DP_BLAST_BIN} || "blastall";
my $dir_blast_db = $ENV{DP_BLAST_DIR} || 0;

#my $blast_path = "blastall";   # Path to the NCBI blastall binary
#                               # by default will assume blastall works
#                               # without providing full path

# Directories
my $dir_gff_out;               # Dir to hold the GFF output
my $dir_blast_out;             # Dir to hold the BLAST output for qry file 
my $indir;                     # Directory containing the seq files to process
my $outdir;                    # Directory to hold the output
my $bac_out_dir;               # Dir for each sequnce being masked

my $msg;                       # Message printed to the log file
my $search_name;               # Name searched for in grep command
my $name_root;                 # Root name to be used for output etc

# ARRAYS
my @dbs;                       # Information for databases
                               # 2d Array filled from the config file

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(
		    # Required Arguments
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    "d|db-dir=s"   => \$dir_blast_db,
		    "c|config=s"   => \$file_config,
		    # Optional strings
		    "blast-path=s" => \$blast_path,
		    "logfile=s"    => \$logfile,
		    # Booleans
		    "verbose"      => \$verbose,
		    "test"         => \$test,
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,
		    "q|quiet"      => \$quiet,);


my $bac_parent_dir = $outdir;  

my ( $ind_lib , $RepMaskCmd, $MakeGffDbCmd, $MakeGffElCmd );
my ( @RepLibs );

#//////////////////////
my $proc_num_max = 5;
#\\\\\\\\\\\\\\\\\\\\\\

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
    print "\nbatch_mask.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# CHECK REQUIRED ARGS         |
#-----------------------------+
if ( (!$indir) || (!$outdir) || (!$file_config) || (!$dir_blast_db) ) {
    print "\a";
    print STDERR "\n";
    print STDERR "ERROR: An input directory was not specified at the".
	" command line\n" if (!$indir);
    print STDERR "ERROR: An output directory was specified at the".
	" command line\n" if (!$outdir);
    print STDERR "ERROR: A configuration file was not specified at the".
	" command line\n" if (!$file_config);
    print STDERR "ERROR: The Blast DB directory must be specified".
	" at the command line to by using the DP_BLAST_DIR variable".
	" in the user environment." if (!$dir_blast_db);
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
    print LOG "  batch_blast.pl\n";
    print LOG "  JOB: $time_now\n";
    print LOG "==================================\n";
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

unless ($dir_blast_db =~ /\/$/ ) {
    $dir_blast_db = $dir_blast_db."/";
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

$count_files = @fasta_files;

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($count_files == 0) {
    print "\a";
    print "\nERROR: No fasta files were found in the input directory\n".
	"$indir\n".
	"Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

print "NUMBER OF FILES TO PROCESS: $count_files\n";


#-----------------------------+
# GET INFO FROM CONFIG FILE   |
#-----------------------------+
# Load to the 2d dbs array
open (CONFIG, $file_config) ||
    die "Can't open the config file:\n$file_config";
$i = 0;
my $line_num = 0;
while (<CONFIG>) {
    $line_num ++;

    unless (m/^\#/) {
	chomp;
	my @tmpary = split (/\t/);
	my $count_tmp = @tmpary;
	#print "COUNT: $count_tmp\n";
	#print "$_\n";
	#print "\t".$tmpary[0]."\n";
	if ($count_tmp == 6) {

	    #-----------------------------+
	    # TEST THE ENTRY              |
	    #-----------------------------+
	    #test_blast_db  ($prog, $db, $line, $dbdir)
	    test_blast_db( $tmpary[0], 
			   $tmpary[4], 
			   $line_num, 
			   $dir_blast_db );

	    $dbs[$i][0] = $tmpary[0];  # Blast Prog
	    $dbs[$i][1] = $tmpary[1];  # Outfile extension
	    $dbs[$i][2] = $tmpary[2];  # Alignment output
	    $dbs[$i][3] = $tmpary[3];  # E Value threshold
	    $dbs[$i][4] = $tmpary[4];  # DBNAME
	    $dbs[$i][5] = $tmpary[5];  # CMD SUFFIX
	    $i++;

	} # End of if count_tmp = 6
	else {
	    print "ERROR: Config file line number $line_num\n";
	    print "       Only $line_num variables were found\n" 
	}

    } # End of unless this is a comment line

} # End of while CONFIG


#-----------------------------+
# CHECK SANITY OF CONFIG FILE |
#-----------------------------+
# This will check for existence of blast database etc
print "Checking config file ...\n" unless $quiet;
for my $ind_db (@dbs) {
#    print "PROG:".$ind_db->[0]."\n";
#    print " EXT:".$ind_db->[1]."\n";
#    print "ALIG:".$ind_db->[2]."\n";
#    print "EVAL:".$ind_db->[3]."\n";
#    print "  DB:".$ind_db->[4]."\n";
#    print "SUFF:".$ind_db->[5]."\n";
#    print "\n\n";
}

$count_db = @dbs;
print "NUMBER OF DATABASES: $count_db\n";
$count_proc = $count_db * $count_files;
print "NUMBER OF PROCESSES: $count_proc\n";

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
	die "ERROR: Could not create the output directory:\n$outdir";
}


#-----------------------------+
# RUN BLAST FOR EACH FILE     |
#-----------------------------+
for my $ind_file (@fasta_files)
{
    
    $file_num++;

    #-----------------------------+
    # Get the root name of the    |
    # file to mask                |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.masked\.fasta$/) {
	# file ends in .masked.fasta
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
    # Create parent dir if it     |
    # does not already exist      |
    #-----------------------------+
    my $dir_parent = $outdir.$name_root."/";
    unless (-e $dir_parent) {
	print "creating dir: $dir_parent\n";
	mkdir $dir_parent ||
	    die "Could not creat the output dir:\n$dir_parent\n";
    }


    #-----------------------------+
    # Create the dir to hold the  |
    # BLAST output                |
    #-----------------------------+
    $dir_blast_out = $outdir.$name_root."/blast/";
    unless (-e $dir_blast_out ) {
	print STDERR "Creating output dir\n: $dir_blast_out\n" 
	    if $verbose;
	mkdir $dir_blast_out ||
	    die "Could not create the output directory:\n$dir_blast_out";
    }
    
    $file_to_blast = $indir.$ind_file;

    #-----------------------------+
    # CREATE THE DIR TO HOLD      |
    # THE GFF OUTPUT              |
    #-----------------------------+
    if ($do_gff) {
	$dir_gff_out = $outdir.$name_root."/gff/";
	unless (-e $dir_gff_out) {
	    print STDERR "Creating output dir\n:$dir_gff_out"
		if $verbose;
	    mkdir $dir_gff_out ||
		die "Could not create the gff output directory:\n$dir_gff_out";
	}
    }
    
    #-----------------------------+
    # FOR EACH DB IN THE CONFIG   |
    # FILE                        |
    #-----------------------------+
    my $gff_count=0;               # Count of the number of gff files
    for my $ind_db (@dbs) {
	
	$proc_num++;

	# Temp exit for debug
	#if ($proc_num > $proc_num_max) { exit; }

	print "\n" unless $quiet;
	print "+-----------------------------------------------------------+\n"
	    if $verbose;
	print "| BLAST PROCESS: $proc_num of $count_proc \n" unless $quiet;
	print "+-----------------------------------------------------------+\n"
	    if $verbose;

	my $file_blast_out = $dir_blast_out.$name_root."_".$ind_db->[4]."."
	    .$ind_db->[1];

	my $blast_cmd = "$blast_path".
	    " -p ".$ind_db->[0].
	    " -i ".$file_to_blast.
	    " -d ".$dir_blast_db.$ind_db->[4].
	    " -e ".$ind_db->[3].
	    " -o ".$file_blast_out.
	    " -m ".$ind_db->[2].
	    " ".$ind_db->[5];
	
	# PRINT VERBOSE OUTPUT
	if ($verbose) {
	    print "\tFILE: $name_root\n";
	    print "\t  DB: ".$ind_db->[4]."\n";
	    print "\t CMD: $blast_cmd\n";
	}

	unless ($test) {
	    system ($blast_cmd);
	}

	# Convert the result to gff
	my $gff_out_file = $dir_gff_out.$name_root."_blast.gff";
	if ($do_gff) {
	    $gff_count++;

	    my $append_gff;
	    if ($gff_count == 1) {
		$append_gff = 0;
	    }
	    else {
		$append_gff = 1;
	    }

	    blast2gff ( $file_blast_out, 
			$gff_out_file,
			$append_gff,
			$name_root,
			$ind_db->[2]);
	}


    } # End of for each ind_db



} # End of for each file in the input folder

close LOG if $logfile;

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

sub test_blast_db {
# Check for the existence of a BLAST database

    # prog ----- blast program to use
    # $db ------ blast database name
    # $line ---- line number from the config file
    # $dbdir --- path to the blast database directory
    # $db_path - path to the expected blast database

    my ($prog, $db, $line, $dbdir) = @_;
    my $db_path;
   
    #-----------------------------+
    # Nucleotide BLAST            |
    #-----------------------------+
    if ( ($prog =~ "blastn") || ($prog =~ "megblast") ||
	 ($prog =~ "tblastx") || ($prog =~ "tblastn")  ) {
	$db_path = $dbdir.$db.".nin";
    }
    #-----------------------------+
    # Protein BLAST               |
    #-----------------------------+
    elsif ( ($prog =~ "blastx") || ($prog =~ "blastp") ||
	    ($prog =~ "psi-blast") || ($prog =~ "phi-blast") ) {
	# Protein BLAST
	$db_path = $dbdir.$db.".pin";
	
    }
    #-----------------------------+
    # Unrecognized BLAST          |
    #-----------------------------+
    else {
	print "\a";
	print "ERROR: Config file line $line\n";
	print "The blast program $prog is not recognized by batch_blast.pl\n";
	exit;
    }

    #-----------------------------+
    # CHECK FOR DB FILE EXISTNECE |
    #-----------------------------+
    unless (-e $db_path) {
	print "\a";
	print "ERROR: Config file line $line\n";
	print "The database file for $db could not be found at:\n".
	    "$db_path\n";
    } else {
	print "The db file is okay at: $db_path\n" if $verbose;
    }

}

sub blast2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # blastin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    # bopt    - The type of blast report to parse (0,8,9)
    my ($blastin, $gffout, $append, $seqname, $bopt ) = @_;
    my $blastprog;        # Name of the blast program used (blastn, blastx)
    my $dbname;           # Name of the database blasted
    my $hitname;          # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit
    my $blast_report;     # The handle for the blast report
    my $blast_score;      # The score for the blast hit

    my $seq_name_len = length($seqname);

    #-----------------------------+
    # GET BLAST REPORT            |
    #-----------------------------+
    if ($bopt == 8 || $bopt == 9) {

	# PARSE m8 or m9 ALIGNMENTS
	$blast_report = new Bio::SearchIO ( '-format' => 'blasttable',
					    '-file'   => $blastin)
	    || die "Could not open BLAST input file:\n$blastin.\n";
	
    }
    else {

	# PARSE DEFAULT ALIGNMNETS (m0)
	$blast_report = new Bio::SearchIO ( '-format' => 'blast',
					    '-file'   => $blastin) 
	    || die "Could not open BLAST input file:\n$blastin.\n";
	
    }


    # Open file handle to the gff outfile    
    if ($append) {
	open (GFFOUT, ">>$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    else {
	open (GFFOUT, ">$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    
    while (my $blast_result = $blast_report->next_result())
    {
	$blastprog = $blast_result->algorithm;

	# FOR m8 output the database name does not work
	# so I will use the file name as given by blastin

	if ($bopt == 8 || $bopt == 9) { 
	    $dbname = $blastin;

	    # Since DAWG-PAWS uses a specific naming stragety, I will
	    # try to extract the database name from the file name
	    # Basic pattern is 'path_root/file name'_'database name'
	    # can get the 
	    if ($dbname =~ m/(.*)\/(.*)\.bln$/ ) {
		$dbname = $2;
		$dbname = substr($dbname, $seq_name_len+1 );

	    }
	    elsif ($dbname =~ m/(.*)\/(.*)\.blx$/ ) {
		$dbname = $2;
		$dbname = substr($dbname, $seq_name_len+1);
	    }

	}
	else {
	    $dbname = $blast_result->database_name();
	}

    	while (my $blast_hit = $blast_result->next_hit())
	{

	    while (my $blast_hsp = $blast_hit->next_hsp())
	    {

		my $hitname = $blast_hit->name();

		$strand = $blast_hsp->strand('query');
		
		if ($strand =~ "-1") {
		    $strand = "-";
		}
		elsif ($strand =~ "1") {
		    $strand = "+";
		}
		else {
		    die "Error parsing strand\n";
		}
		
		# Make certain that start coordinate is
		# less then the end coordinate
		if ( $blast_hsp->start() < $blast_hsp->end() ) {
		    $start = $blast_hsp->start();         # Start
		    $end = $blast_hsp->end();
		    #$strand = "+";
		    
		}
		else {
		    #
		    $start = $blast_hsp->end();         # Start
		    $end = $blast_hsp->start();
		    #$strand = "-";
		}


		#-----------------------------+
		# GET BLAST SCORE             |
		#-----------------------------+
		if ($bopt == 8 || $bopt == 9) {
		    #$blast_score = ".";
		    # trying bits
		    $blast_score = $blast_hsp->bits();
		}
		else {
		    $blast_score = $blast_hsp->score();
		}

		#-----------------------------+
		# PRINT OUTPUT TO GFF FILE    |
		# WITH BAC DATA               |
		#-----------------------------+
		# Changing BLASTN to the Bac Name appears to allow 
		# these to be drawn on different levels.
		print GFFOUT 
		    "$seqname\t".                     # Seqname
		    "$blastprog:$dbname\t".           # Source
		    "exon\t".                         # Feature type name
		    "$start\t".                       # Start
		    "$end\t".                         # End
		    $blast_score."\t".                # Score
		    "$strand\t".                      # Strand
		    ".\t".                            # Frame
		    "$hitname\n";                     # Feature name


		# PRINT GFF TO STDERR IF IN VERBOSE MODE
		if ($verbose) {
		    print STDERR "\t   SEQ:\t$seqname\n";
		    print STDERR "\t SOURC:\t$blastprog:$dbname\n";
		    print STDERR "\t START:\t$start\n";
		    print STDERR "\t   END:\t$end\n";
		    print STDERR "\t SCORE:\t".$blast_score."\n";
		    print STDERR "\tSTRAND:\t$strand\n";
		    print STDERR "\t   HIT:\t$hitname\n";
		}

	    } # End of while next hsp
	} # End of while next hit
    } # End of while next result
    
    close GFFOUT;
    
}

1;
__END__

=head1 NAME

batch_blast.pl - Do NCBI-BLAST searches for a set of fasta files

=head1 VERSION

This documentation refers to batch_blast version $Rev: 747 $

=head1 SYNOPSIS

=head2 Usage

    batch_blast.pl -i DirToProcess -o OutDir -d DbDir -c ConfigFile

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -d, --db-dir   # Directory to hold program output
    -o, --outdir   # Path to the base output directory
    -c, --config   # Path to the config file

=head1 DESCRIPTION

Given a directory of softmasked fasta files, this will
BLAST the files against a the set of BLAST formatted
databases specified in the configuration file.

All of the BLAST output files will be stored in a directory
name for the intput file in a subdirectory named blast.
(ie /home/myhome/infile/blast/infile_blastdb.blo).

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the directory containing the sequences to process.

=item -o,--outdir

Path of the directory to place the program output.

=item -d,--db-dir

Path of the directory containing the blast formatted databases.

=item -c, --config

Path to the batch_blast config file. This is a tab delimited text file
indicating the required information for each of the databases to blast
against. Lines beginning with # are ignored.

=back

=head1 OPTIONS

=over 2

=item --blast-path

Full path to the NCBI blastall program. Default is blastall.

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

=item --verbose

Run the program with maximum output.

=item -q,--quiet

Run the program with minimal output.

=item --test

Run the program without doing the system commands.

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

=head2 Configuration File

The location of the configuration file is indicated by the --config option
at the command line.
This is a tab delimited text file
indicating required information for each of the databases to blast
against. Lines beginning with # are ignored, and data are in six 
columns as shown below:

=over 2

=item Col 1. Blast program to use [ tblastx | blastn | blastx ]

The blastall program to use. DAWG-PAWS will support blastn,
tblastx, and blastx format.

=item Col 2. Extension to add to blast output file. (ie. bln )

This is the suffix which will be added to the end of your blast
output file. You can use this option to set different extensions
for different types of blast. For example *.bln for blastn
output and *.blx for blastx output.

=item Col 3. Alignment output options (-m options from blast)

DAWG-PAWS supports output in the default pairise format (0),
the tabular format (8), and the tabular format with
comments (9). See the blastall documentation for a full list
of all possible alignment output options.

=item Col 4. Evalue threshold

The maximum e-value to include in the blast report

=item Col 5. Database name

The name of the database that is being blasted against.

=item Col 6. Additional blast command line options

This is the place to indicate additional options in your BLAST command such 
as multiple processors (-a 2) or use lowercase filtering (-U). Options
should be space separated.
For a list of all options available in blast, type blastall --help
at the comand line.

=back

An example config file:

 #-----------------------------+
 # BLASTN: TIGR GIs            |
 #-----------------------------+
 blastn	bln	8	1e-5	TaGI_10	-a 2 -U
 blastn	bln	8	1e-5	AtGI_13	-a 2 -U
 blastn	bln	8	1e-5	ZmGI_17	-a 2 -U
 #-----------------------------+
 # TBLASTX: TIGR GIs           |
 #-----------------------------+
 tblastx	blx	8	1e-5	TaGI_10	-a 2 -U
 tblastx	blx	8	1e-5	AtGI_13	-a 2 -U
 tblastx	blx	8	1e-5	ZmGI_17	-a 2 -U

=head2 Environment

The following variables in the user environment can be set for
batch_blast.pl:

=over 2

=item DP_BLAST_BIN

This is the location of the blastall binary that you would like to use 
for the batch_blast.pl program. If this is not specified at the command
line the program will assume that the 'blastall' program is located
in the user's PATH. Alternatively, you can specify the location of the
blastall program using the --blast-path option at the command line. The
command line use of --blast-path will override the path set by 
DP_BLAST_BIN.

=item DP_BLAST_DIR

This is the location of the directory that the compiled blast databases
are stored in. Alternatively, you can specify the location of this
directory using the command line option --db-dir. The use of the --db-dir
option in the command line will override the directory selection set
by the DP_BLAST_DIR option.

=back

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * NCBI blastall

The latest version of the NCBI blastall program can be downloaded from:
ftp://ftp.ncbi.nih.gov/blast/executables/LATEST

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

=item * Limited to NCBI BLAST

The current version is limited to using the NCBI version of BLAST.

=item * Config file must use UNIX format line endings

The config file must have UNIX formatted line endings. Because of
this any config files that have been edited in programs such as
MS Word must be converted to a UNIX compatible text format before
being used with batch_blast.

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

STARTED: 07/23/2007

UPDATED: 04/23/2009

VERSION: Release 1.0

=cut

#-------------------------------------------------------+
# HISTORY                                               |
#-------------------------------------------------------+
#
# 07/23/2007
# - Program started
# - Imported functions from batch_mask.pl as well as
#   Blast_PAWS.pl
# - Added use of a config file to define attributes of
#   BLAST jobs for each database to BLAST against
#   where lines starting with # are ignored
#
# 07/24/2007
# - Adding test_blast_db to test the existence of blast
#   databases
#
# 12/03/2007
# - Updating program POD documentation 
# - Changed version to SVN Revision Number
#
# 12/05/2007
# - Added new print_help subfunction
# - Moved POD documentation to the end
#
# 04/23/2009
# - Added code to convert to gff format
