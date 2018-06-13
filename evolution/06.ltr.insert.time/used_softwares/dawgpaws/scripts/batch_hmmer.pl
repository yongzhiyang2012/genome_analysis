#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_hmmer.pl - Run hmmer searches in batch mode.        |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/17/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run hmmer against repeat hmmer models in batch mode,     |
#  produce results in gff format if requested.              |
#                                                           | 
# VERSION: $Rev: 768 $                                      |
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
use Bio::Tools::HMMER::Results;# Bioperl HMMER results parser
use Text::Wrap;                # Allows word wrapping and hanging indents
                               # for more readable output for long strings.
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
my $indir;
my $outdir;
my $config_file;
my $name_root;                  # Root name of the sequence file
my $i;

# Array of rephmmer parameters
my @bh_params = ();            # Batch hmmer parameters

# Booleans
my $quiet = 0;
my $do_gff = 1; 
my $test = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# Counters
my $file_num = 0;
my $db_parent_dir;


#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"   => \$indir,
                    "o|outdir=s"  => \$outdir,
		    "c|config=s"  => \$config_file,
		    # ADDITIONAL OPTIONS
		    "db-parent=s" => \$db_parent_dir,
		    "gff"         => \$do_gff,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
		    "test"        => \$test,
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
    print "\nbatch_hmmer.pl:\nVersion: $VERSION\n\n";
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
if ( (!$indir) || (!$outdir) || (!$config_file) ) {
    print "\a";
    print STDERR "ERROR: An input directory must be specified\n" 
	if (!$indir);
    print STDERR "ERROR: An output directory must be specified\n" 
	if (!$outdir);
    print STDERR "ERROR: A config file must be specified\n" 
	if (!$config_file);
    print_help ("usage", $0 );
    exit;
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

if ($db_parent_dir) {
    unless ($db_parent_dir =~ /\/$/ ) {
	$db_parent_dir = $db_parent_dir."/";
    }
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
# LOAD THE CONFIG FILE        |
#-----------------------------+
$i=0;
my $config_line_num=0;

open (CONFIG, "<$config_file") ||
    die "ERROR Can not open the config file:\n $config_file";

while (<CONFIG>) {
    chomp;
    $config_line_num++;
    unless (m/^\#/) {
       	my @in_line = split (/\t/);     # Implicit split of $_ by tab
	my $num_in_line = @in_line; 
	
	if ($num_in_line == 3) { 
	    $bh_params[$i][0] = $in_line[0] || "NULL";  # Name
	    $bh_params[$i][1] = $in_line[1] || "NULL";  # HMMER model dir
	    $bh_params[$i][2] = $in_line[2] || "NULL";  # HMMER suffix

	    # Append db_parent dir to dir path if present
	    if ($db_parent_dir) {
		$bh_params[$i][1] = $db_parent_dir.$bh_params[$i][1];
	    }

	    $i++;

	    if ($verbose) {
		print STDERR "NAME:\t".$bh_params[$i][0]."\n";
		print STDERR "DIR:\t".$bh_params[$i][1]."\n";
		print STDERR "".$bh_params[$i][2]."\n";
	    } # End of if verbose
	}
	else {
	    print "\a";
	    print STDERR "WARNING: Unexpected number of lines in config".
		" file line $config_line_num\n$config_file\n";
	}

   } # End of unless comment line
} # End of while CONFIG file
close CONFIG;

# Number of parameter sets specified in the config file
my $num_par_sets = $i;

if ($num_par_sets == 0) {
    print "\a";
    print "No parameter sets were found in the config file:\n";
    print "$config_file\n\n";
    print "batch_hmmer.pl --man \n";
    print "to see how to format configuration files\n";
    exit;
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
my $num_files = @fasta_files;

if ($num_files == 0) {
    print "\a";
    print "No fasta files were found in the input direcotry:\n";
    print "$indir\n";
    print "Fasta file MUST have the fasta or fa extension to be".
	" recognized as fasta files\n";
    exit;
}

my $num_proc_total = $num_files * $num_par_sets;

#-----------------------------+
# RUN hmmsearch FOR EACH OF   |
# THE FASTA FILES IN THE      |
# INPUT DIRECTORY             |
#-----------------------------+

for my $ind_file (@fasta_files) {
    
    $file_num++;
    
    # Get root name
    if ($ind_file =~ m/(.*)\.masked\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }

    #-----------------------------+
    # MAKE OUTPUT DIR             |
    #-----------------------------+
    # The base output dir for the BAC
    my $bac_out_dir = $outdir.$name_root."/";
    mkdir $bac_out_dir, 0777 unless (-e $bac_out_dir); 

    #-----------------------------+
    # MAKE HMMER OUTPUT DIR       |
    #-----------------------------+
    # dir for the repeat maske output
    my  $hmmer_out_dir = "$bac_out_dir"."hmmer/";
    mkdir $hmmer_out_dir, 0777 unless (-e $hmmer_out_dir); 
    
    #-----------------------------+
    # MAKE GFF OUTPUT DIR         |
    #-----------------------------+
    my $gff_out_dir = "$bac_out_dir"."gff/";
    mkdir $gff_out_dir, 0777 unless (-e $gff_out_dir); 


    my $file_to_search = $indir.$ind_file;

    #-----------------------------+
    # FOR EACH PARAM SET IN THE   |
    # bh_params ARRAY             |
    #-----------------------------+
    for ($i=0; $i<$num_par_sets; $i++) {
 	
	my $hmm_par_name = $bh_params[$i][0];
	my $hmm_par_model_dir = $bh_params[$i][1];
	my $hmm_par_suffix = $bh_params[$i][2]; 

	#my $hmm_out_path = $hmmer_out_dir.$name_root."_hmm_".
	#    $hmm_par_name.".hmmout";
	my $hmm_in_path = $indir.$ind_file;
    
	#-----------------------------+
	# CREATE DIR TO HOLD PARAMETER|
	# SET HMMSEARCH OUTPUT        |
	#-----------------------------+
	my $hmmer_par_out_dir = $hmmer_out_dir.$hmm_par_name."/";
	mkdir $hmmer_par_out_dir, 0777 unless (-e $hmmer_par_out_dir); 
	
	run_hmmer( $name_root, $file_to_search, $hmm_par_model_dir,
#		   $hmm_par_suffix, $hmmer_out_dir, $hmm_par_name);
		   $hmm_par_suffix, $hmmer_par_out_dir, $hmm_par_name);
		

	if ($do_gff) {
	    my $gff_out_path = $gff_out_dir.$name_root."_hmm_".
		$hmm_par_name.".gff";
	    
	    # 0 at end indicates to not append to an existing gff
	    hmmer2gff ($hmmer_par_out_dir, $gff_out_path, 0,
		       $name_root, $hmm_par_name);
	    
	} # End of if $do_gff
	

    } # End for each i in the bh_params array
    
    
} # End of for each ind_file in the fasta_files array

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

sub run_hmmer {
# Run the HMMER program for a set of models

    #-----------------------------+
    # VARS PASSED TO THE SUBFUN   |
    #-----------------------------+
    my ($seq_name, $seq_path, $model_dir, $hmm_suffix,
	$out_dir, $par_name) = @_;

    # $seq_name   # name of the seq record id HEX0014K09
    # $seq_path   # Path to the seq 
    # $model_dir  # Path to the dir containing the seq models 
    # $hmm_suffix # Options to add to the hmmsearch command line
    #             # May need to chomp the suffix
    # $out_dir    # Out dir for the location of the hmm out files
    # $par_name   # Name given to the hmm parameter set

    # Append slash to end of model_dir if needed
    unless ($model_dir =~ /\/$/ ) {
	$model_dir = $model_dir."/";
    }

    unless ($out_dir =~ /\/$/ ) {
	$out_dir = $out_dir."/";
    }

    #-----------------------------+
    # VAR SCOPE AND INITIALIZE    |
    #-----------------------------+
    my $sub_proc_num = 0;           # Process number starts at zero

    my $hmm_cmd;                   # hmmsearch command line 
#    my $out_path = $out_dir.$seq_name."_hmm_".$par_name.".hmmout";
    
    # Will need to delete an already existing outdir	

    #-----------------------------+
    # OPEN MODEL DIR              |
    #-----------------------------+
    opendir( MODDIR, $model_dir ) ||
	die "Can not open the hmm model directory:\n $model_dir\n";
    # The following does not search for the hmm extension
    # and assumes that all files in this directory are hmm models

    # The following will only use files with the hmm extension
    my @models = grep /hmm?$/, readdir MODDIR;

    closedir( MODDIR );    

    my $num_models = @models;   # Number of models
    
    for my $ind_model (@models) {
	
	$sub_proc_num++; # Increment the process number

	# Model root name will be used in the naming of files
	# default is just the name of the individual model file
	# will try to trim off the hmm when it is present
	my $model_root_name = $ind_model;

	if ($ind_model =~ m/(.*)\.hmm$/) {
	    $model_root_name = $1;
	}

 	# HMMER MODEL PATH
	my $model_path = $model_dir.$ind_model;
	my $out_path = $out_dir.$seq_name."_hmm_".$par_name."_".
	    $model_root_name.".hmmout";

	#-----------------------------+
	# DOES THE HMM MODEL EXIST    |
	#-----------------------------+
	if ( (-e $model_path) && (-e $seq_path) ) {
	    

	    print "\nMOD PATH: $model_path\n" if $verbose;
	    print "SEQ PATH: $seq_path\n" if $verbose;
	    print "OUT PATH: $out_path\n" if $verbose;

#	    $hmm_cmd = "hmmsearch " . 
#		"--domT 2 $model_path $seq_path >$out_path";
	    
	    # Alternative may be to set a single $out_path
	    # However, the hmmer parser may have a problem with this
	    $hmm_cmd = "hmmsearch " . 
		"--domT 2 $model_path $seq_path >>$out_path";

	    #-----------------------------+
	    # PRINT STATUS OF EACH HMM    |
	    # PROCESS                     |
	    #-----------------------------+	
	    if ($verbose) {
		print "HMM model ".$sub_proc_num." of ".$num_models."\n";
		print wrap("\t", "\t", $hmm_cmd );
		print "\n";
	    }
	 
	    # Run the hmm command if this is not a test
	    system ($hmm_cmd) if (!$test);

	}
	else {
	    
	    # Show error if the model_path or input sequence
	    # could not be located
	   
	    unless (-e $model_path) {
		print STDERR "ERROR: Can not find model:\n$ind_model\n";
	    }
	    
	    unless (-e $seq_path) {
		print STDERR "ERROR: Can not find sequense:\n$seq_path\n";
	    }

	} # End of test for existence of model_path and seq_path
       
	
    } # End of for each database loop


} # End of run hmmer subfunction


sub hmmer2gff {
#-----------------------------+
# PARSE THE OUTPUT FROM A     |
# HMMER RUN AGAINST A SET OF  |
# HMM PROFILES                |
#-----------------------------+ 

    #-----------------------------+
    # VARIABLES                   |
    #-----------------------------+
    # THE FOLLOWING SHOULD BE PASSED TO THE SUBFUNCTION

    my ($hmm_in_dir, $gff_out_path, $append_gff, 
	$seq_name, $gff_suffix) = @_;

    # hmm_in_dir     # The directory containing the hmm files
                     # This will need to be all output for the parameter set
    # gff_out_path   # Out file for the gff file
    # append_gff     # Boolean, should the gff file be appended
                     # By default this will append.
    # seq_name       # Name of the seq file ie HEX0014K09
    # gff_suffix     # Suffix to add to the hmmsearch:suffix
                     # This will be the parameter set name

# Old variables
#    my $seq_name = $_[0];
#    my $WorkDir = $_[1];
#    my $HmmDir = $_[2];
#    my $class = $_[3];

    # Set to subfun scope
    # Added 09/18/2007
    my $domain;
    my $dataset;
    my $tot_res;  # Total results
    my $BestName = "";          # Name of the best hit
    my $FilePath;                # The file path for the inidividual HMM output
    my $TotRes = '0';            # Total number of results for the sequence
    my $BestBit = '0';           # Best BIT Score
    my $BestHit;                 # The name of the Best Hit
    my $BestEval = 'null'; # Set scope of this variable

    #-----------------------------+
    # LOAD FILES TO PARSE INTO    |
    # THE @hmm_files ARRAAY       |
    #-----------------------------+
    # Currently gets all files in the directory
    # May need to create a subdir for each basic model
    # parameter set.
    
    unless ($hmm_in_dir =~ /\/$/ ) {
	$hmm_in_dir = $hmm_in_dir."/";
    }
    # THIS CURRENTLY ASSUMES THAT ALL FILES IN hmm_in_dir ARE 
    # HMMSEARCH OUTPUT FILES
    opendir( HMMDIR, $hmm_in_dir ) ||
	die "ERROR: Can not open hmm input dir:\n$hmm_in_dir\n";
    my @hmm_files = grep !/^\.\.?$/, readdir HMMDIR ;
    closedir( HMMDIR );
    
    # Open up an Output file
    if ($append_gff) {
	open ( GFFOUT, ">>$gff_out_path") ||
	    die "ERROR: Can not open gff output file:\n$gff_out_path\n";
    }
    else {
	open ( GFFOUT, ">$gff_out_path") ||
	    die "ERROR: Can not open gff output file:\n$gff_out_path\n";;
    }

    $tot_res = 0;
    #-----------------------------------------------------------+
    # FOR EACH FILE IN THE ARRAY OF FILES                       |
    #-----------------------------------------------------------+
    for my $ind_file (@hmm_files)
    {
	
	my $hmm_in_path = $hmm_in_dir.$ind_file;
	
        # OPEN THE HMM RESULT AS HMMER RESULTS OBJECT
	my $hmm_res = new Bio::Tools::HMMER::Results ( -file => $hmm_in_path ,
						       -type => 'hmmsearch') 
	    || die "Could not open file\n$hmm_in_path\n";
	
	my $num_res = $hmm_res->number;
	$tot_res = $tot_res + $num_res;
	
	# ONLY PRINT OUTPUT FOR QUERIES WITH MATCHES
	if ($num_res >> 0) #08/07/2006
	{	              #08/07/2006
	    foreach my $seq ( $hmm_res->each_Set ) 
	    {
		foreach $domain ( $seq->each_Domain ) 
		{		
		    #my $CurBit = $domain->bits;
		    my $cur_name = $domain->hmmname;
		    $cur_name =~ m/.*\/(.*)\.hmm/;  # Returns what I want
		    $cur_name = $1;
		    # RECORD THE NAME AND SCORE OF THE
		    # BEST HIT
		    if ($domain->bits > $BestBit)
		    {
			# ASSUMES BIT SCORE AND E VALUE
			# HAVE THE SAME RANK
			$BestBit = $domain->bits;
			$BestName = $cur_name;
			$BestEval = $domain->evalue || "NULL";
		    }

		    # Determine the length of the match of the
		    # element length
		    #my $dom_len = $domain->end - $domain->start;
		    
		    my $gff_source = "hmmer:".$gff_suffix;
		    
		    #-----------------------------+
		    # PRINT RESULTS TO GFF OUT    |
		    #-----------------------------+
		    print GFFOUT "$seq_name\t".    # Name of sequence
			"$gff_source\t".           # Source
			"transposable element\t".  # Feature name
			$domain->start."\t".       # Feature start
			$domain->end."\t".         # Feature end
			$domain->evalue."\t".      # Feature score
			".\t".                     # Feature strand
			".\t".                     # Feature frame
			"$cur_name\n";             # Feature name
		    
		} # End of for each domain in the HMM output file
	    } # End of if greater then zero hits 08/07/2006
	    
	} # End of for each seq in the HMM output file
	
    }
    
#    if ($verbose) {
#	print "\n\t===================================\n";
#	print "\t$class HMMER RESULT\n";
#	print "\t===================================\n";
#	print "\tSEQ:       \t".$seq_name."\n";
#	print "\tDATASET:   \t".$dataset."\n";
#	print "\tCLASS:     \t".$class."\n";
#	print "\tTOTAL:     \t".$tot_res."\n";
#	print "\tBEST HIT:  \t".$BestName."\n";
#	print "\tBEST BIT:  \t".$BestBit."\n";
#	print "\tBEST EVAL: \t".$BestEval."\n";
#	print "\t===================================\n";
#    }

    
    close (GFFOUT);

} # END OF THE CONVERT TO GFF SUBFUNCTION

1;
__END__

# Deprecated print_help function

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

batch_hmmer.pl - Run hmmsearch in batch mode

=head1 VERSION

This documentation refers to program version $Rev: 768 $

=head1 SYNOPSIS

=head2 Usage

    batch_hmmer.pl -i InDir -o OutDir -c Config.cfg [--gff]

=head2 Required Arguments

    --indir         # Path to the input dir containing fasta files
    --outdir        # Path to the base output directory
    --config        # Configuration file for running hmmer
    --gff           # Path to gff output file

=head1 DESCRIPTION

Given a config file this will run the hmmsearch program for each parameter
set in the configuration file for each fasta sequence file in the input
directory. This will also convert the results to GFF format if requested.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the input directory. This is the directory that contains
the 

=item -o,--outdir

Path of the base output directory.

=item -c,--config

Path of the configuration file specifying the parameters used to run
the hmmer program. This config file includes (1) the name of the
parameter set, (2) the directory of the hmm models, and (3)
command line arguments for hmmsearch.

=back

=head1 OPTIONS

=over

=item --db-parent

The parent directory that the directories listed in the configuration
file are located in.

=item --gff

Produce a GFF output file for each of the query sequence for each
of the parameter sets in the input file.

=item -q,--quiet

Run the program with minimal output.

=item --verbose

Run the program in verbose mode.

=item --version

Show program version.

=item --usage

Short overview of how to use program from command line.

=item --help

Show program usage with summary of options.

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

The location of the configuration file is indicated by the --config
option at the command line. This is a tab delimited text file. 
The columns of data represent  the following information

=over 2

=item col 1.

Name of the parameter set

=item col 2. Dir of the hmmer models

All hmm profile models in this directory will be searched using the
hmmer arguments in col 3.

=item col 3.Arguments for hmmer

The arguments for hmmer should be separate by spaces and can include the 
following arguments:

=over 2

=item -A <n>

sets alignment output limit to <n> best domain alignments

=item -E <x>

sets E value cutoff (globE) to <= x

=item -T <x>

sets T bit threshold (globT) to >= x

=item -Z <n>

sets Z (# seqs) for E-value calculation

=item --compat

make best effort to use last version's output style

=item --cpu <n>

run <n> threads in parallel (if threaded)

=item --cut_ga

use Pfam GA gathering threshold cutoffs

=item --cut_nc

use Pfam NC noise threshold cutoffs

=item --cut_tc 

use Pfam TC trusted threshold cutoffs

=item --domE <x>

sets domain Eval cutoff (2nd threshold) to <= x

=item --domT <x>

sets domain T bit thresh (2nd threshold) to >= x

=item --forward

use the full Forward() algorithm instead of Viterbi

=item --informat <s>

sequence file is in format <s>

=item --null2

turn OFF the post hoc second null model

=item --pvm

run on a Parallel Virtual Machine (PVM)

=item --xnu

turn ON XNU filtering of target protein sequences

=back

=back

An example that will do hmmsearch for four directories of hmm models
using an evalue threshold of 0.00001 follows:

    rice_mite  /$HOME/HMMData/db/hmm/rice_mite_models/	-E 0.00001
    rice_mule  /$HOME/HMMData/db/hmm/rice_mule_models/	-E 0.00001
    tpase      /$HOME/HMMData/db/hmm/tpase_models/	-E 0.00001
    pfam       /$HOME/HMMData/db/hmm/pfam/              -E 0.00001

=head1 DEPENDENCIES

=head2 Required Software

=over 2

=item * HMMER

The hmmer program is required: http://hmmer.janelia.org/. Specifically
this program requires the hmmsearch program.

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the BLAST results.

=item * Getopt::Long

This module is required to accept options at the command line.

=item * Text::Wrap

This module allows for word wrapping and hanging indents of printed
text. This allows for a more readable output for long strings

=item * Bio::Tools::HMMER::Results

This module is part of the BioPerl package. This allows for
parsing of the the output from the hmmsearch program.

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

=item * Must created profile hmms external

The DAWG-PAWS package currently does not include programs for 
creating hmm profile models. These must be created externally using
the hmmer package of programs.

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

STARTED: 09/17/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/17/2007
# - Program started
# - Basic boilerplate stuff laid out
# - Main body of program written with subfunctions included
#   from existing code.
# 
# 12/13/2007
# - Updated POD documentation
# - Added print_help subfunction that extracts help and usage
#   messages from the POD documentation
#
# 04/23/2009
# - Modified to always produce gff results
