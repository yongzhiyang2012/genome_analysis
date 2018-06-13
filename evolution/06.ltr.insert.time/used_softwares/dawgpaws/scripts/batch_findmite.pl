#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_findmite.pl - Batch run the findmite program        |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 08/30/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run the FINDMITE program in batch mode using a config    |
#  file to allow the program to run in a high-throughput    |
#  batch mode.                                              |
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


# Array to hold the findmite parameters
my @fm_params = ();

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;                     # Input directory
my $outdir;                    # Output directory
my $parfile;                   # Parameters file
my $name_root;                 # Root name of the input file
my $gff_dir;                    # Dir to hold the gff output file

# Counters/Index Vals
my $i = 0;                     # Array index val
my $file_num = 0;              # File number
my $proc_num = 0;              # Process number

# Booleans
my $verbose = 0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_gff_convert = 1;
my $do_fasta_out = 0;          # Create a fasta file with all predicted MITEs

my $findmite_bin = $ENV{FINDMITE_BIN} || "FINDMITE1New.bin";

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# Required Arguments
		    "i|indir=s"       => \$indir,
                    "o|outdir=s"      => \$outdir,
		    "c|config=s"      => \$parfile,
		    # Additional options
		    "f|fasta"         => \$do_fasta_out,
		    "q|quiet"         => \$quiet,
		    "verbose"         => \$verbose,
		    "gff"             => \$do_gff_convert,
		    # Additional information
		    "usage"           => \$show_usage,
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

if ($show_version) {
    print "\nbatch_findmite.pl:\nVersion: $VERSION\n\n";
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
if ( (!$indir) || (!$outdir) || (!$parfile) ) {
    print "ERROR: Input directory path required\n" if (!$indir);
    print "ERROR: Output directory path required\n" if (!$outdir);
    print "ERROR: Config file required\n" if (!$parfile);
    print_help ("usage", $0 );
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

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
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print "Creating output dir ...\n" if $verbose;
    mkdir $outdir ||
	die "Could not create the output directory:\n$outdir";
}

#-----------------------------+
# LOAD FINDMITE PARAMATERS    |
#-----------------------------+
open (PARFILE, "<$parfile")
    || die "Can not open the paramter file:\n$parfile\n";

while (<PARFILE>) {
    chomp;

#    # Use the following if an error crops up later
#    print "IN: $_\n";
    
    # Ignores comment lines starting with #
    unless (m/^\#/) {
	#my @in_line = split(/\t/, $_);
	my @in_line = split;           # Implicit split of $_ by whitespace
	my $num_in_line = @in_line; 
	$fm_params[$i][0] = $in_line[0] || "NULL"; 
	$fm_params[$i][1] = $in_line[1] || "NULL";
	$fm_params[$i][2] = $in_line[2] || "NULL";
	$fm_params[$i][3] = $in_line[3] || "NULL";
	$fm_params[$i][4] = $in_line[4] || "NULL";
	$fm_params[$i][5] = $in_line[5] || "NULL";
	$fm_params[$i][6] = $in_line[6] || "NULL";
	$fm_params[$i][7] = $in_line[7] || "NULL";
	$fm_params[$i][8] = $in_line[8] || "NULL";
	$fm_params[$i][9] = $in_line[9] || "NULL";
	
#        # Use the following if a error crops up later
#	if ($verbose) {
#	    print STDERR "\tProcessed as:\n";
#	    print STDERR "\t".$fm_params[$i][0]."\n";
#	    print STDERR "\t".$fm_params[$i][1]."\n";
#	    print STDERR "\t".$fm_params[$i][2]."\n";
#	    print STDERR "\t".$fm_params[$i][3]."\n";
#	    print STDERR "\t".$fm_params[$i][4]."\n";
#	    print STDERR "\t".$fm_params[$i][5]."\n";
#	    print STDERR "\t".$fm_params[$i][6]."\n";
#	    print STDERR "\t".$fm_params[$i][7]."\n";
#	    print STDERR "\t".$fm_params[$i][8]."\n";
#	    print STDERR "\t".$fm_params[$i][9]."\n";
#	}
	
	$i++;
    }

} # End of while INFILE

my $num_par_sets = $i;
my $max_i = $i-1;

close PARFILE;


# TO DO: CHECK PARAM FILE VALIDITY


#-----------------------------+
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) || 
    die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $num_files = @fasta_files;

my $num_proc_total = $num_files * $num_par_sets;

print STDERR "$num_proc_total findmite runs to process\n";

for my $ind_file (@fasta_files) {

    $file_num++;
    
    # Get root file name
    if ($ind_file =~ m/(.*)\.masked\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fa$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    #-----------------------------+
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root."/";
    unless (-e $name_root_dir) {
	mkdir $name_root_dir ||
	    die "Could not create dir:\n$name_root_dir\n"
    }

    #-----------------------------+
    # CREATE FINDMITE OUTDIR      |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $findmite_dir = $name_root_dir."findmite/";
    unless (-e $findmite_dir) {
	mkdir $findmite_dir ||
	    die "Could not create genscan out dir:\n$findmite_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    if ($do_gff_convert) {
	$gff_dir = $name_root_dir."gff/";
	unless (-e $gff_dir) {
	    mkdir $gff_dir ||
		die "Could not create genscan out dir:\n$gff_dir\n";
	}
    }

    #-----------------------------------------------------------+
    # RUN FINDMITE FOR EACH SET OF PARAM VALUES                 |
    #-----------------------------------------------------------+
    for ($i=0; $i<$num_par_sets; $i++) {
	$proc_num++;
	
	# Load array vals to usefull short names
	my $mite_name    = $fm_params[$i][0];
	my $rep          = $fm_params[$i][1];
	my $tir_len      = $fm_params[$i][2];
	my $num_mismatch = $fm_params[$i][3];
	my $filt_at      = $fm_params[$i][4];
	my $filt_cg      = $fm_params[$i][5];
	my $filt_atta    = $fm_params[$i][6];
	my $filt_2_base  = $fm_params[$i][7];
	my $min_dist     = $fm_params[$i][8];
	my $max_dist     = $fm_params[$i][9];

	if ($verbose) {
	    print STDERR "\n\n========================================\n";
	    print STDERR " FindMITE Process $proc_num of $num_proc_total\n";
	    print STDERR " ".$name_root."_".$mite_name."\n";
	    print STDERR "========================================\n";
	}

	#-----------------------------+
	# RUN THE FINDMITE PROGRAM    |
	#-----------------------------+
	# Then should be made a subfunction that can return true
	my $fm_in = $indir.$ind_file;
	my $fm_out = $findmite_dir.$name_root."_".$mite_name.".mite.txt";
#	my $fm_cmd = "FINDMITE1New.bin $fm_in $fm_out";
	my $fm_cmd = "$findmite_bin $fm_in $fm_out";
	print "CMD: $fm_cmd\n" if $verbose;


	# It seems like FINDMTE is closing
	# Some j
	open(FINDMITE,"|$fm_cmd") ||
	    die "Could not open findmite\n";
	#sleep 1;
	print FINDMITE "150000\n" ||
	    die "Could not set length";
	print FINDMITE "$rep\n";
	#sleep 1;
	print FINDMITE "$tir_len\n";
	#sleep 1;
	print FINDMITE "$num_mismatch\n";
	#sleep 1;
	print FINDMITE "$filt_at\n";
	#sleep 1;
	print FINDMITE "$filt_cg\n";
	#sleep 1;
	print FINDMITE "$filt_atta\n";
	#sleep 1;
	print FINDMITE "$filt_2_base\n";
	#sleep 1;
	print FINDMITE "$min_dist\n";
	#sleep 1;
	print FINDMITE "$max_dist\n";
	#sleep 10;
	close FINDMITE;


	#-----------------------------+
	# CONVERT OUTPUT TO GFF       |
	#-----------------------------+
	if ($do_gff_convert) {
	    my $findmite_gff = $findmite_dir.$name_root.
		"_".$mite_name."mite.gff";
	    my $findmite_fasta = $outdir."FINDMITE_".
		$mite_name."mite.fasta";

	    if (-e $fm_out) {
		# Currently will do_fasta_file creation if requested
		# This last boolean set to 1 will append to any existing
		# fasta file that is present
		findmite2gff ($fm_out, $findmite_gff, 0, $name_root, 
			      $do_fasta_out, 1, $findmite_fasta);
		
		# Copy the gff file file
		my $gff_copy_path = $gff_dir.$name_root.
		    "_".$mite_name."mite.gff";
		my $cp_err_msg = "Can not copy from:\n".
		    " $findmite_gff TO $gff_copy_path\n";
		copy ($findmite_gff, $gff_copy_path) ||
		    die "$cp_err_msg\n";
		

	    }
	    else {
		print STDERR "\nERROR: Could not find the findmite out file\n".
		    "$fm_out\nto convert to gff\n";
	    }
	}

	# print return to clean up the command line
	print "\n";

    } # End of run find mite for each set of param vals
    
    

} # End of for each $ind_file


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub findmite2gff {
    
    #-----------------------------+
    # SUBFUNCTION VARS            |
    #-----------------------------+
    my ($findmite_in, $gff_out, $append_gff, $seqname, 
	$fasta_out, $append_fasta, $fasta_out_path) = @_;

    # GFF OUT VARS
    my $gff_seq_id;                 # Seq id for use in gff
    my $gff_start;                  # Start of the feature
    my $gff_end;                    # End of the feature
    my $gff_source;                 # Source ie findmite_at_11
    my $gff_name;                   # Name of the feature 

    # FINDMITE HEADER VARS
    my $fm_dir_repeats;             # Direct repeats
    my $fm_tir_len;                 # Lenght of the TIR
    my $fm_mismatch_max;            # Max number of mistmatches
    my $fm_filter_at;               # Were A/T strings filtered 
    my $fm_filter_cg;               # Were C/G strings filtered
    my $fm_filter_atta;             # Were AT/TA strings filtered
    my $fm_filter_2_base;           # Percent (0 to 100)
    my $fm_min_dist;                # Minimum distance used
    my $fm_max_dist;                # Maximum distance used
    my $fm_file_an;                 # Name of the file analyzed

    # FINDMATE INDIVIDUAL MITE VARS
    my $fm_pattern;                 # Pattern
    my $fm_seq_id;                  # Id of the query sequence
    my $fm_num_mismatch;            # Num of mismatches
    my $fm_seq = "";                # Sequence string as parsed from findmite
    my $mite_seq_string;            # Sequence string of the putatitve mite
    my $mite_context_5;             # 5' Context of the mite (40bp)
    my $mite_context_3;             # 3' Context of the mite (40bp)
    my $mite_start;                 # Start of the mite as parsed from FM
    my $mite_end;                   # End of the mite as parsed from FM

    # Counter
    my $mite_num = 0;               # Incremented ID number for the mite

    # BOOLEANS
    my $in_seq = 0;                 # Boolean, in seq data (past header info)

    #-----------------------------+
    # OPEN FILES                  |
    #-----------------------------+

    open (INFILE, "<$findmite_in") ||
	die "Can not open input file:\n$findmite_in\n";

    
    if ($append_gff) {
	open (GFFOUT, ">>$gff_out") ||
	    die "Can not open output file:\n$gff_out\n";	
    }
    else {
	open (GFFOUT, ">$gff_out") ||
	    die "Can not open gff output file:\n$gff_out\n";
    }


    if ($fasta_out) {
#	my $fasta_out_path = $outdir."FINDMITE_".
#	    $fm_dir_repeats."_".$fm_tir_len.".fasta";

	if ($append_fasta) {
	    open (FASTAOUT, ">>$fasta_out_path") ||
		die "Can not FASTA output file:\n$fasta_out_path\n";
	} 
	else {
	    open (FASTAOUT, ">$fasta_out_path") ||
		die "Can not FASTA output file:\n$fasta_out_path\n";
	}
    } # End of if $fasta_out


    #-----------------------------+
    # PROCESS FINDMITE FILE       |
    #-----------------------------+
    while (<INFILE>) {
	chomp;
 	#print $_."\n";

	if(m/^Pattern  : (.*)/) {

	    # If we have previously loaded seq data then
	    # futher parse the sequence string and and
	    # print the results to the gff output
	    if ($in_seq) {
		                   
		if ($fm_seq =~ m/(.*)\((.*)\)\-{5}(.*)\-{5}\((.*)\)(.*)/) {
		    $mite_context_5 = $1;
		    $mite_start = $2;
		    $mite_seq_string = $3;
		    $mite_end = $4;
		    $mite_context_3 = $5;

		    if ($verbose) {
			print STDERR "\t5CON: $mite_context_5\n";
			print STDERR "\tSTAR: $mite_start\n";
			print STDERR "\tMITE: $mite_seq_string\n";
			print STDERR "\tEND : $mite_end\n";
			print STDERR "\t3CON: $mite_context_3\n";
			print STDERR "\n\n";
		    }

		    #-----------------------------+
		    # PRINT TO GFF FILE           |
		    #-----------------------------+
		    # Parse seq id for shorter name
		    if ($fm_seq_id =~ m/^>(\S*)\s./) {
			$gff_seq_id = $1;
		    }
		    elsif ($fm_seq_id =~ m/^>(.*)/) {
			$gff_seq_id = $1;
		    }
		    else {
			$gff_seq_id = $fm_seq_id;
		    }

		    $gff_source = "findmite:".
			$fm_dir_repeats."_".$fm_tir_len;
		    $gff_name = $gff_seq_id."_".$fm_dir_repeats.
			"_".$fm_tir_len."_".$mite_num;
		    $gff_start = $mite_start;
		    $gff_end = $mite_end;

		    print GFFOUT 
			"$gff_seq_id\t".            # Seq name
			"$gff_source\t".            # Source
			"mite\t".                   # Feature type
			"$gff_start\t".             # Start
			"$gff_end\t".               # End
			".\t".                      # Score
			".\t".                      # Strand
			".\t".                      # Frame
			"$gff_name\n";              # Feature name

		    #-----------------------------+
		    # PRINT MITE TO FASTA FILE    |
		    #-----------------------------+
		    if ($fasta_out) {
			print FASTAOUT ">$gff_name\n";
			print FASTAOUT "$mite_seq_string\n";
		    }

		}

		# Reset vals to null
		$fm_seq = "";
	    }

	    $in_seq = 1;
	    $mite_num++;
	    $fm_pattern = $1;
	    print STDERR "$fm_pattern\n" if $verbose;
	}
	elsif(m/^Sequence : (.*)/) {
	    $fm_seq_id = $1;
	    print STDERR "\t$fm_seq_id\n" if $verbose;
	}
	elsif(m/^Mismatch : (.*)/){
	    $fm_num_mismatch = $1;
	    print STDERR "\t$fm_num_mismatch\n" if $verbose;
	}
	elsif($in_seq) {
	    $fm_seq = $fm_seq.$_;
	}
	#-----------------------------+
	# HEADER INFORMATION          | 
	#-----------------------------+
	elsif(m/Direct repeats.{17}(.*)/) {
	    $fm_dir_repeats = $1;
	}
	elsif(m/Length of TIR.{18}(.*)/) {
	    $fm_tir_len = $1;
	}
	elsif(m/Number of mis.{18}(.*)/) {
	    $fm_num_mismatch = $1;
	}
	elsif(m/Filtering A\/T.{18}(.*)/) {
	    $fm_filter_at = $1;
	}
	elsif(m/Filtering C\/G.{18}(.*)/) {
	    $fm_filter_cg = $1;
	}
	elsif(m/Filtering AT\/.{18}(.*)/) {
	    $fm_filter_atta = $1;
	}
	elsif(m/Filtering 2.{20}(.*)/) {
	    $fm_filter_2_base = $1;
	}
	elsif(m/Minimum dist.{19}(.*)/) {
	    $fm_min_dist = $1;
	}
	elsif(m/Maximum dist.{19}(.*)/) {
	    $fm_max_dist = $1;
	}
	elsif(m/The results from the input.{17}(.*)/) {
	    $fm_file_an = $1;
	}

    }

    
    #-----------------------------+
    # PRINT OUT OF DATA IN VARS   |
    #-----------------------------+
    if ($fm_seq =~ m/(.*)\((.*)\)\-{5}(.*)\-{5}\((.*)\)(.*)/) {
    #if ($fm_seq =~ m/(.*)\((.*)\)\-\-\-\-\-(.*)\-\-\-\-\-\((.*)\)(.*)/) {
	$mite_context_5 = $1;
	$mite_start = $2;
	$mite_seq_string = $3;
	$mite_end = $4;
	$mite_context_3 = $5;
	

	if ($verbose) {
	    print STDERR "\t5CON: $mite_context_5\n";
	    print STDERR "\tSTAR: $mite_start\n";
	    print STDERR "\tMITE: $mite_seq_string\n";
	    print STDERR "\tEND : $mite_end\n";
	    print STDERR "\t3CON: $mite_context_3\n";
	    print STDERR "\n\n";
	}	

    }

    #-----------------------------+
    # PRINT TO GFF FILE           |
    #-----------------------------+
    # Parse seq id for shorter name
    if ($fm_seq_id =~ m/^>(\S*)\s./) {
	$gff_seq_id = $1;
    }
    elsif ($fm_seq_id =~ m/^>(.*)/) {
	$gff_seq_id = $1;
    }
    else {
	$gff_seq_id = $fm_seq_id;
    }
    
    $gff_source = "findmite:".
	$fm_dir_repeats."_".$fm_tir_len;
    $gff_name = $gff_seq_id."_".$fm_dir_repeats.
	"_".$fm_tir_len."_".$mite_num;
    $gff_start = $mite_start;
    $gff_end = $mite_end;
    
    print GFFOUT 
	"$gff_seq_id\t".            # Seq name
	"$gff_source\t".            # Source
	"mite\t".                   # Feature type
	"$gff_start\t".             # Start
	"$gff_end\t".               # End
	".\t".                      # Score
	".\t".                      # Strand
	".\t".                      # Frame
	"$gff_name\n";              # Feature name
    
    #-----------------------------+
    # PRINT MITE TO FASTA FILE    |
    #-----------------------------+
    if ($fasta_out) {
	print FASTAOUT ">$gff_name\n";
	print FASTAOUT "$mite_seq_string\n";
    }

    # Print while debu
    if ($verbose) {
	print STDERR "FILE  : $fm_file_an\n";
	print STDERR "DIRREP: $fm_dir_repeats\n";
	print STDERR "TIRLEN: $fm_tir_len\n";
	print STDERR "NUMMIS: $fm_num_mismatch\n";
	print STDERR "F_AT  : $fm_filter_at\n";
	print STDERR "F_CG  : $fm_filter_cg\n";
	print STDERR "F_ATTA: $fm_filter_atta\n";
	print STDERR "F_2BAS: $fm_filter_2_base\n";
	print STDERR "MIN   : $fm_min_dist\n";
	print STDERR "MAX   : $fm_max_dist\n";
    }

    close INFILE;
    close GFFOUT;
    close FASTAOUT;

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

# Deprecated print_help
sub print_help {

    # Print requested help or exit.
    # Options are to just print the full 
    my ($opt) = @_;

    my $usage = "USAGE:\n". 
	"batch_findmite.pl -i InDir -o OutDir -p ParamFile";
    my $args = "REQUIRED ARGUMENTS:\n".
	"  --indir        # Directory of fasta files to process\n".
	"  --outdir       # Path to the output direcotry\n".
	"  --params       # Path to the params file\n".
	"\n".
	"OPTIONS:\n".
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

batch_findmite.pl - Run the findmite program in batch mode.

=head1 VERSION

This documentation refers to batch_findmite version $Rev: 791 $

=head1 SYNOPSIS

=head2 Usage

    batch_findmite.pl -i InDir -o OutDir -c ConfigFile [-gff]

=head2 Required Arguments

    -i, --indir    # Directory of fasta files to process
    -o, --outdir   # Path to the base output directory
    -c, --config   # Path to the config file
    --gff          # Produce output in GFF format

=head1 DESCRIPTION

The batch_findmite program will do a FINDMITE analysis for each parameter
set in your configuration file for each query sequence in your 
input directory. The results from FINDMITE have a VERY high false positive
rate so you will need to further evaluate your results to find the true MITEs
in your query sequence.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the input directory. This is a directory that contains the fasta
files to process.

=item -o,--outdir

Path of the base output directory.

=item -c,--config

Path to the configuration file that includes the parameter sets
for running the FINDMITE program.
These parameters represent the answers to the series of
questions that must be answered when running FINDMITE. Any lines
starting with # are ignored.

B<EXAMPLE>

  #------------------------------------------------------------------
  # Name      	Rep	TIR  Mis AT GC ATTA	2Base	Min	Max
  #------------------------------------------------------------------
  TA_11		TA	11   1	 y  y  y	85	30	700
  TA_12		TA	12   1   y  y  y	85	30	700	

For a detail description, see the Paramters File heading under 
the CONFIGURATION AND ENVIRONMENT section of the full program 
documentation.

=back

=head1 OPTIONS

=over 2

=item -f,--fasta

Create a fasta file of all predicted MITEs. A different fasta file will be 
created for each of the parameter set names. This will currently append to
the existing fasta data file in the output directory.

=item -gff

Produce a GFF output file that indicates where the predicted MITEs are on
the query sequence.

=item -q,--quiet

Run the program with minimal output. Does not require user interaction.

=item --verbose

Run the program with maximal output.

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

The path to the configuration file is indicated at the command line with
-c or --config.

This file is a space delimited text
file that indicates the parameters to use when running the findmite
program. These parameters represent the answers to the series of
questions that must be answered when running FINDMITE.

B<EXAMPLE>

  #------------------------------------------------------------------
  # Name      	Rep	TIR  Mis AT GC ATTA	2Base	Min	Max
  #------------------------------------------------------------------
  TA_11		TA	11   1	 y  y  y	85	30	700
  TA_12		TA	12   1   y  y  y	85	30	700	

The columns above represent the following information:

=over 2

=item Col. 1

Base name to assign to putative mites

=item Col. 2

Direct Repeat

=item Col. 3

Length of the Terminal Inverted Repeat (TIR)

=item Col. 4

Number of mismatches

=item Col. 5

Boolean to fileter the A/T.
This must be set to y or n. 

=item Col. 6

Boolean to Filter C/G
This must be set to y or n.

=item Col. 7

Boolean to filter AT/TA
This must be set to y or n.

=item Col. 8

Proporiton of 2Base to filter.
This must be an integer between 0 and 100.

=item Col. 9

Minimum distance between TIRs
This must be an integer.

=item Col. 10

Maximum distance between TIRs
This must be an integer.

=back

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * FINDMITE

The batch_findmite.pl program is dependent on the FINDMITE program. 
A version of 
FINDMITE compiled for RedHat Linux is available at:
L<http://jaketu.biochem.vt.edu/dl_software.htm>

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the BLAST results.

=item * Getopt::Long

This module is required to accept options at the command line.

=back

=head1 BUGS AND LIMITATIONS

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over

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

STARTED: 08/30/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/30/2007
# - Program started
# 08/31/2007
# - Added the gff conversion subfunction from cnv_findmite2gff.pl
# 09/06/2007
# - Updating POD documentation.
# - Deleted the empty run_findmite subfunction
# 09/28/2007
# - Change variable for the file giving parameter sets from
#   param to config
# - Made gff name SONG complient
# - Added throw error and exit when the required variables not present
# - Adding code to copy the gff output to the gff directory
#   for the contig
# - Moved POD documentation to the end of the program
#
# 12/13/2007
# - Updated POD documentation
# - Added print_help subfunction to extract help and
#   usage messages from the POD documentation
