#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# vennseq.pl - Venn diagram overlap of sequence features    |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 03/06/2007                                       |
# UPDATED: 04/06/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Creates a Venn Diagram comparing the overlap of sequence |
#  features along a query sequence. The purpose is to allow |
#  me to visualize the overlap of Repeat Databases or gene  |
#  models along a BAC. This is also used to compare gene    |
#  annotation curation results among curators               |
#                                                           |
# DEPENDENCIES:                                             |
#  -BioPERL                                                 |
#  -NCBI BLAST                                              |
#  -Venn Master                                             |
#                                                           |
# USAGE:                                                    |
#   vennseq.pl -i HEX3045G05.fasta -o test_out_txt.txt      |
#              -s test_out_svg.svg -d work_gff -f gff       |
#                                                           |
# NOTES:                                                    |
#  -All array indices start at one.                         |
#   This will facilitate parsing information from sequence  |
#   features that are all indexed from 1 to length of seq   |
#  -q is quiet, Q is super quiet                            |
#                                                           |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# To do.- If no fasta file then the sequence length may
#         be passed as a variable from the command line.
#       - Take a single gff file as input       
#
# TEST CMD: vennseq.pl -i HEX3045G05_TREP9.masked.fasta -o test_out_txt.txt -s test_out_svg.svg -d work_gff -f gff
#
# This is a test of svn 

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use Getopt::Long;              # Get cmd line options
use Bio::SearchIO;             # Parse BLAST output
use Bio::SeqIO;                # Parse fasta sequence input file
#use GD;                       # Draw using the GD program
#                              # In case I want to draw a coverage map
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
# SET VARIABLE SCOPE          |
#-----------------------------+
my $infile;                    # Path to sequence file in fasta format
my $outfile;                   # Output file path for Venn text file
my $svg_outfile;               # Output file path for Venn SVG image output
my $feature_dir;               # Directory containing sequence features
my $feat_format = "gff";       # Format of the sequence feature files
                               #  - blast or b
                               #  - tab or t
                               #  - gff or g

# VARS WITH EVN OPTIONS
my $vmaster_dir = $ENV{VMASTER_DIR} || 
    "/Applications/VennMaster-0.37.3/";
my $java_bin = $ENV{VMASTER_JAVA_BIN} || 'java';
my $java_mem = $ENV{VMASTER_JAVA_MEM} || 512;

# BOOLEANS
my $show_help = 0;             # Show program help
my $show_version = 0;          # Show program version
my $show_man = 0;              # Show program manual page using peldoc
my $show_usage = 0;            # Show program usage command             
my $quiet = 0;                 # Boolean for reduced output to STOUT
my $test = 0;                  # Run the program in test mode
my $verbose = 0;               # Run the program in verbose mode
my $do_audit_blast = 0;        # Audit the blast output
my $do_audit_rm = 0;           # Audit the repeatmasker output
my $do_audit_ta = 0;           # Audit the TriAnnotation output
my $do_full_audit = 0;         # Audit everything
my $debug = 0;
my $no_launch_vm = 0;          # Don't launch the venn master program

my $seq_len;                   # Full length of the sequence that the
                               # features are being mapped onto
my @feat_files;                # Array to hold the path of the 
                               # feature files.
my @feature_files;             # The original list of potential feature
                               # files. The zero length files will
                               # be ignored leaving only files with data
my @empty_files;                # Files with no hits (0 length files)
my $num_feat_files;              # The number of feature files
my @FeatMatrix;                # Count of feature occurrences
my @BinFeatMatrix;             # Presence/abscence binary form of the
                               # feature matrix
my @FeatLabels;                # Category labels for the features
                               # this will take the 
my @CrossTab;                  # Cross table matrix
my $FileTestNum = 0;           # Initialize counter of file test number
my $EmptyFileNum = 0;          # Initialize counter of empty files
my $seq_id;

#-----------------------------------------------------------+
# COMMAND LINE VARIABLES                                    |
#-----------------------------------------------------------+


my $ok = GetOptions(
                    # Required Variables
		    "i|infile=s"     => \$infile,
                    "o|outfile=s"    => \$outfile,
                    "d|dir=s"        => \$feature_dir,
    		    # Options
                    "f|format=s"     => \$feat_format,
                    "s|svg-out=s"    => \$svg_outfile,  
                    "venn-master=s"  => \$vmaster_dir,  
                    "no-vm"          => \$no_launch_vm,
                    "java-path"      => \$java_bin,
                    "java-mem"       => \$java_mem,
                    "seq-id"         => \$seq_id,
		    "seq-len"        => \$seq_len,
		    # Booleans
                    "debug"          => \$debug,
		    "verbose"        => \$verbose,
		    "test"           => \$test,
		    "usage"          => \$show_usage,
		    "version"        => \$show_version,
		    "man"            => \$show_man,
		    "h|help"         => \$show_help,
		    "q|quiet"        => \$quiet,
                    );


my $Usage = "\nVennSeq.pl -i InputFastaFile -o OutputFile\n".
	"-d FeatureFileDirecotry \n";

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
    print "\nvennseq.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------+
# REQUIRED VARIABLE CHECK     |
#-----------------------------+
if ( (!$outfile) || (!$feature_dir) ) {
    print "\a";
#    print "ERROR: The fasta format infile  must be specified at the".
#	" command line\n" if !$infile;
    print "ERROR: The outfile must be specified at the command line\n"
	if !$outfile;
    print "ERROR: The feature directory must be specified at the command line\n"
	if !$feature_dir;
 #   print "ERROR: The feature file format must be specified at the command".
	#" line\n" if !$feat_format;
    exit;
}

#-----------------------------+
# CHECK FOR DIR SLASH         |
#-----------------------------+
unless ($feature_dir =~ /\/$/ ) {
    $feature_dir = $feature_dir."/";
}

#-----------------------------+
# CHECK FILE EXISTENCE        |
#-----------------------------+
# This is the fasta sequence file
unless (-e $infile) {
    print STDERR "ERROR:The input sequence file could not be found:\n$infile\n";
}

#-----------------------------+
# OPEN THE OUTPUT FILE        |
#-----------------------------+
open (OUT, ">$outfile") ||
    die "Could not open output file for writing:\n$outfile\n";

#-----------------------------+
# GET LENGTH OF THE SEQUENCE  |
#-----------------------------+
unless ($seq_len) {
    # This will be used to determine the length of the array
    # This includes an internal check for only one sequence
    # in the input file.
    my $seq_in = Bio::SeqIO->new('-file' => "<$infile",
				 '-format' => 'fasta') ||
				     die "Can not open input file:n\$infile\n";
    
    my $num_seq = 0;
    while (my $inseq = $seq_in->next_seq) {
	$num_seq++;
	$seq_len = $inseq->length();
	if ($num_seq >> 1) {
	    print "ERROR. Input sequence file contains more then one sequence.\n";
	    exit;
	}
    }
    
}


#-----------------------------------------------------------+
# GET SEQUENCE FEATURES                                     |
#-----------------------------------------------------------+

#-----------------------------+
# TAB DELIM BLAST OUTPUT      |
#-----------------------------+
# Currently only supporting the blo extension
if ($feat_format =~ "tab") {

    opendir(INDIR, $feature_dir) ||
	die "Can not open input directory:\n$feature_dir: $!";
    @feature_files = grep /blo$/, readdir INDIR;

    # TO DO, CHECK NUMBER FEATURE FILES, IF NONE FOUND
    # ATTEMPT TO USE ALL FILES IN THE INPUT DIRECTORY

    #-----------------------------+
    # TEST FOR FILE LENGTH        |
    #-----------------------------+

    # Check for files with no hits here by checking for zero length
    # 
    # Can load indexes to delete ??
    foreach my $FileTest (@feature_files) {
	$FileTestNum++;
	
	#LOAD ZERO LENGTH FILE TO EmptyFiles array
	# Otherise push them to the FeatFiles array

        # If the file is zero length
	if (-z $feature_dir.$FileTest) {
	    $EmptyFileNum++;
#	    $EmptyFiles[$EmptyFileNum][1] = $EmptyFileNum;
#	    $EmptyFiles[$EmptyFileNum][2] = $FileTestNum;
#	    $EmptyFiles[$EmptyFileNum] = &GetShortName($FileTest);
#	    $EmptyFiles[$EmptyFileNum][3] = $FileTest;
	    # Load the short name of the database to the 
	    # empty files array, if the GetShort Name does
	    # not return a database then I will use the File Name
	    push (@empty_files, &GetShortName($FileTest) || $FileTest);
	}
	else{
	    push (@feat_files, $FileTest); 
	}# End of if FileSize = 0
    } # End of for each file in File test
    
    #-----------------------------+
    # PRINT OUT THE EMPTY FILES   |
    #-----------------------------+
    # This is a place to drop empty files from the array
    if ($EmptyFileNum > 0) {
	print STDERR "The following files had no BLAST hits:\n";
	foreach my $IndEmptyFile (@empty_files) {
	    print STDERR "\t".$IndEmptyFile."\n";
	}
    }
    
    $num_feat_files = @feat_files;
    print "\n$num_feat_files files have hits\n";

    if ($num_feat_files == 0) {
	print "No feature files of type \"$feat_format\" were found in:\n".
	    "$feature_dir\n";
	exit;
    }
    else {
	print "$num_feat_files of type $feat_format were found in\n".
	    "$feature_dir\n";
    }

    #-----------------------------+
    # FILL THE FEATURE MATRIX     |
    # WITH ZEROS                  |
    #-----------------------------+
    print STDERR "Initializing feature matrix and binary matrix with zeros\n" 
	if $verbose;
    for (my $row=1; $row<=$num_feat_files; $row++) {
	my $NumCol=0;
	
	print STDERR "\tInitializing row $row..." if $verbose;
	for (my $col=1; $col<=$seq_len; $col++) {
	    $FeatMatrix[$row][$col] = 0;
	    $BinFeatMatrix[$row][$col] = 0;
	    $NumCol++;
	}
	print STDERR "$NumCol cols.\n" if $verbose;
    }

    #-----------------------------+
    # FOR EACH OF THE FEATURE     |
    # FILES INCREMENT THE PROPER  |
    # ROW IN @FeatMatrix          |
    #-----------------------------+
    # Use uppercase Row and Col here to avoid
    # any problems with var names
    my $Row = 0;                   # Reset row count to zero
                                   # Each feature has a unique row
    print STDERR "Loading data from feature files.\n" if $verbose;
    foreach my $IndFile (@feat_files) {
	$Row++;                    # Increment row for each file

	$FeatLabels[$Row] = &GetShortName($IndFile);

	print "\tLoading data from feature file $FeatLabels[$Row]\n";

	my $FeatPath = $feature_dir.$IndFile;
	my $NumLines = 0;                  # Number of lines in the feat file
	my $BinFeatCount = 0;              # Counter for the "binary" features

	open (FEATIN, "<$FeatPath") ||
	    die "Can not open the feature file:\n$FeatPath\n";
	while (<FEATIN>) {
	    chomp;
	    $NumLines++;

	    # Use unless to ignore comment lines from -m 9 output
	    unless (m/^\#.*/) {
		# TEST THAT THE LENGTH OF THE ARRAY IS 
		# WHAT IS EXPECT
		my @TabSplit = split(/\t/);
		my $TestLen = @TabSplit;
		if ($TestLen > '12') {
		    print STDERR "ERROR: The BLAST file has an unexpected".
			" number of tab delimited columns.\n";
		}

		my ($QryId, $SubId, $PID, $Len,
		    $MisMatch, $GapOpen,
		    $QStart,$QEnd, $SStart, $SEnd,
		    $EVal, $BitScore) = split(/\t/);
		
		#-----------------------------+
		# PRINT SOME FEATURE VALUES   |
		# FOR DEBUG PURPOSES          |
		#-----------------------------+
		if ($debug) {
		    if ($NumLines < '3') {
			print "\n\t\tLINE: $NumLines\n";
			print "\t\t   QRY: $QryId\n";
			print "\t\t   SUB: $SubId\n";
			print "\t\tQSTART: $QStart\n";
			print "\t\t  QEND: $QEnd\n\n";
		    } # End of If NumLines less then three
		} # End of unless quiet
		
		for (my $Col=$QStart; $Col<=$QEnd; $Col++) {
		    $BinFeatCount++;
		    $FeatMatrix[$Row][$Col] = $FeatMatrix[$Row][$Col] + 1;
		} # End of for my $col from Query Start to End
	    } # End of unless # (unless comment line)
	} # End of while FEATIN
	close FEATIN;

	#-----------------------------+
	# REPORT NUMBER OF FEATURES   |
	# THAT WERE PROCESSED         |
	#-----------------------------+
	# For tab delimited blast the number of features
	# is roughly equal to the number of lines == number of HSPS
	#print "\t\tFILE: $FeatPath\n";
	unless ($quiet) {
	    print STDERR "\t\tFILE: $IndFile\n";
	    print STDERR "\t\tFEAT: $NumLines\n";
	    print STDERR "\t\t BIN: $BinFeatCount\n";
	} # End of unless quiet
    } # End of For each $IndFile in @feat_files

}

#-----------------------------------------------------------+
# GFF FEATURE FILE FORMAT                                   |
#-----------------------------------------------------------+
# Much of this is redundant with the code above

elsif ($feat_format =~ "gff") {

    opendir(INDIR, $feature_dir) ||
	die "Can not open input directory:\n$feature_dir: $!";
    @feature_files = grep /gff$/, readdir INDIR;

    #-----------------------------+
    # TEST FOR FILE LENGTH        |
    #-----------------------------+
    # Check for files with no hits here by checking for zero length
    # 
    # Can load indexes to delete ??
    foreach my $FileTest (@feature_files) {
	$FileTestNum++;
	
	#LOAD ZERO LENGTH FILE TO EmptyFiles array
	# Otherise push them to the FeatFiles array
	if (-z $feature_dir.$FileTest) {
	    $EmptyFileNum++;
	    push (@empty_files, &GetShortName($FileTest) || $FileTest);
	}
	else {
	    push (@feat_files, $FileTest); 
	} # End of if FileSize = 0

    } # End of for each file in File test

    #-----------------------------+
    # CHECK FOR EMPTY FILES       |
    #-----------------------------+
    # This is a place to drop empty files from the array
    if ($EmptyFileNum > 0) {
	if ($verbose) {
	    print STDERR "The following files had no data:\n";
	    foreach my $IndEmptyFile (@empty_files) {
		print STDERR "\t".$IndEmptyFile."\n";
	    }
	}
    }
    
    $num_feat_files = @feat_files;
    print STDERR "\n$num_feat_files files have data\n" if $verbose;

    if ($num_feat_files == 0) {
	print STDERR "No feature files of type \"$feat_format\"".
	    " were found in:\n$feature_dir\n" if $verbose;
	exit;
    }
    else {
	print STDERR "$num_feat_files files of type $feat_format were found in\n".
	    "$feature_dir\n" if $verbose;
    }


    #-----------------------------+
    # FILL THE FEATURE MATRIX     |
    # WITH ZEROS                  |
    #-----------------------------+
    print STDERR "Initializing feature matrix and binary matrix with zeros\n" 
	if $verbose;
    for (my $row=1; $row<=$num_feat_files; $row++) {
	my $NumCol=0;              # Varialbe to count number of cols
	                           # for debug purpsoses
	print STDERR "\tInitializing row $row..." if $verbose;
	for (my $col=1; $col<=$seq_len; $col++) {
	    $FeatMatrix[$row][$col] = 0;
	    $BinFeatMatrix[$row][$col] = 0;
	    $NumCol++;
	}
	print STDERR "$NumCol cols.\n" if $verbose;
    }

    #-----------------------------+
    # LOAD DATA FROM FEAT FILE    |
    #-----------------------------+
    my $Row = 0;                   # Reset row count to zero
                                   # Each feature has a unique row
    print STDERR "Loading data from feature files.\n" if $verbose;
    foreach my $IndFile (@feat_files) {
	$Row++;                    # Increment row for each file

	$FeatLabels[$Row] = &GetShortName($IndFile);

	print STDERR "\tLoading data from feature file $FeatLabels[$Row]\n"
	    if $verbose; 

	my $FeatPath = $feature_dir.$IndFile;
	my $NumLines = 0;                  # Number of lines in the feat file
	my $BinFeatCount = 0;              # Counter for the "binary" features

	open (FEATIN, "<$FeatPath") ||
	    die "Can not open the feature file:\n$FeatPath\n";

	while (<FEATIN>) {
	    chomp;
	    $NumLines++;
	    
	    unless (m/^\#.*/) {
	
		my ($seqname, $source, $feature, $start, $end,
		    $score, $strand, $frame, $attributes) = split(/\t/);
		
		#-----------------------------+
		# PRINT SOME FEATURE VALUES   |
		# FOR DEBUG PURPOSES          |
		#-----------------------------+
		# Just do this for the first few hits
		if ($verbose) {
		    if ($NumLines < '3') {
			print STDERR "\n\t\t  LINE: $NumLines\n";
			print STDERR "\t\t   seq: $seqname\n";
			print STDERR "\t\t   src: $source\n";
			print STDERR "\t\t  feat: $feature\n";
			print STDERR "\t\t start: $start\n";
			print STDERR "\t\t   end: $end\n";
		    } # End of If NumLines less then three
		} # End of unless quiet
		
		for (my $Col=$start; $Col<=$end; $Col++) {
		    $BinFeatCount++;
		    $FeatMatrix[$Row][$Col] = $FeatMatrix[$Row][$Col] + 1;
		} # End of for my $col from Query Start to End
	    } # End of unless # (unless comment line)
	} # End of while FEATIN
	close FEATIN;
	
	#-----------------------------+
	# REPORT NUMBER OF FEATURES   |
	# THAT WERE PROCESSED         |
	#-----------------------------+
	# For tab delimited blast the number of features
	# is roughly equal to the number of lines == number of HSPS
	#print "\t\tFILE: $FeatPath\n";
	print STDERR "\t\tFILE: $IndFile\n" if $verbose;
	print STDERR "\t\tFEAT: $NumLines\n" if $verbose;
	print STDERR "\t\t BIN: $BinFeatCount\n" if $verbose;
    } # End of For each $IndFile in @feat_files

}

#-----------------------------------------------------------+
# UNRECOGNIZED FEATURE FILE FORMAT                          |
#-----------------------------------------------------------+
else {
    print "A proper sequence feature file format was not selected.";
    exit;
} # END OF IF $feat_format

#-----------------------------------------------------------+
# FILL THE BINARY MATRIX AND PRINT OUTPUT FOR VENN MASTER   |
#-----------------------------------------------------------+

print "Filling the binary matrix.\n" if $verbose;
for ($row=1; $row<=$num_feat_files; $row++) {
    print STDERR "\tLoading row $row\n" if $verbose;

    my $BinFeatLen = 0;                    # Binary feature length
    my $BinFeatOccur = 0;                  # Occurrence of features

    for ($col=1; $col<=$seq_len; $col++) {
	$BinFeatLen++;
	# If the feature matrix contains data
	if ($FeatMatrix[$row][$col] >> 0) {
	    $BinFeatOccur++;
	    $BinFeatMatrix[$row][$col] = 1;
	    
	    if ($seq_id) {
		print OUT "$seq_id$col\t$FeatLabels[$row]\n";
	    }
	    else {
		print OUT "$col\t$FeatLabels[$row]\n";
	    }
	} 
    } # End of for $col

    #-----------------------------+
    # REPORT THE PERCENT COVERAGE |
    # OF THE FEATURE AFTER        |
    # TESTING FOR DIVIDE BY ZERO  |
    #-----------------------------+
    if ($BinFeatOccur == 0) {
	print STDERR "\t\tERROR: The feature has zero length.\n" if $verbose;
    } 
    else {
	# Calculate feature coverage
	my $BinFeatCov = $BinFeatOccur/$BinFeatLen;
	# Convert feature coverage to percent
	$BinFeatCov = sprintf("%.4f",$BinFeatCov);
	print STDERR "\t\tLEN: $BinFeatLen\tCOV: $BinFeatCov\n" if $verbose;
    }

} # End of for row

# close the VENN Master text file to make it available
close OUT;

#-----------------------------------------------------------+
# GENERATE COMPARISON MATRIX @CrossTab                      |
#-----------------------------------------------------------+
# This should be made and option, ie. don't do under quiet
# Generates the NxN CrossTab matrix 
# Where N is the number of sequence features being comparedh zeros
print STDERR "Generating the Cross tab\n" if $verbose;

print STDOUT "\n";

#-----------------------------+
# PRINT TABLE HEADER          |
#-----------------------------+


#-----------------------------+
# CREATE SEPERATOR STRING AND |
# PRINT TOP OF CROSS TABLE    |
#-----------------------------+
my $Sep = "+";
for (my $i=1; $i<=$num_feat_files ; $i++) {
    $Sep = $Sep."--------+";
}
print STDOUT $Sep."\n";

#-----------------------------+
# PRINT BODY OF CROSS TABLE   |
#-----------------------------+
for (my $i=1; $i<=$num_feat_files ; $i++) {
    print STDOUT "|"; # Print the left hand side of the cross tab

    for (my $j=1; $j<=$num_feat_files ; $j++) {
	my $SharedCount=0;    # Initialize the count of shared
	my $iCount=0;         # Initialize the i feature count
	my $CrossRatio; # Set scope of the cross ratio

	if ($i != $j) {

	    #-----------------------------+
	    # CALCULATE CROSS TABLE RATIO |
	    #-----------------------------+
	    # k is the sequence position
	    for (my $k=1; $k<=$seq_len; $k++) {
		if ($BinFeatMatrix[$i][$k] == 1) {
		    # Increment the i feature count
		    $iCount++;
		    if ($BinFeatMatrix[$j][$k] == 1) {
			# Increment the shared feature count
			$SharedCount++;
		    } # End of if j = 1 in binary matrix
		} # End of if i = 1 in binary matrix
	    } # End of for k
	    
	    #-----------------------------+
	    # LOAD RESULT TO CROSS TABLE  |
	    # AFTER DIVIDE BY ZERO TEST   |
	    #-----------------------------+
	    if ($iCount >0) {
		$CrossRatio = $SharedCount/$iCount;
		$CrossRatio = sprintf("%.4f",$CrossRatio);
	    }
	    else {
		# If if = 0 this is the NULL case
		$CrossRatio = "  NULL  ";
	    } # End of if $iCount > 0
	    
	    # Load to the Cross Tab Array
	    $CrossTab[$i][$j]=$CrossRatio;	    
	    # To print in table format
	    print STDOUT " ".$CrossRatio." |";
	}
	else {
	    # If i = j then print the data source
	    # label in the table, this will need to
	    # be formatted to fit in the box
	    # - sign below indicates left justification
	    my $DatSrcLab = sprintf ('%-*s', 8, $FeatLabels[$i]);
	    print STDOUT $DatSrcLab."|";
	}
	
    } # End of for j (Seq feature set two)

    print STDOUT "\n"; # Print new line for end of row in cross tab
    print STDOUT $Sep."\n";
} # End of for i (Seq feature set one)

#print $Sep."\n";
print STDOUT "\n"; # Print final new row at the very end of cross tab

#-----------------------------------------------------------+
# OPEN VENN MASTER TO DISPLAY VENN DIAGRAMS                 |
#-----------------------------------------------------------+
print STDERR "Running the VennMaster program.\n" if $verbose;


# NOTE: VENN MASTER OPTIONS
#-help,-?
#    displays help information
#--version, -v
#    displays the version number of VennMaster
#--cfg file.xml
#    configuration file
#--ocfg file.xml
#    configuration file output
#--list file.list
#    list file input
#--gce input.gce
#    GoMiner gene-category export file
#--se input.se
#    GoMiner summary export file
#--htgce ht-input.gce
#    High-Throughput GoMiner gce file
#--filter file.xml
#    filter file (for GoMiner)
#--ofilter file.xml
#    filter file output
#--svg file.svg
#    SVG file output
#--sim file.txt
#    simulation errors output
#--prof file.txt
#    error profile output for the last simulation 


#-----------------------------+
# RUN VENN MASTER             |
#-----------------------------+
my $vm_cmd; # Command to run in vmatch
if (!$svg_outfile) {
    $vm_cmd = "$java_bin -Xms256m -Xmx256m -jar $vmaster_dir".
	"venn.jar --list $outfile $@";
}
else {
    $vm_cmd = "$java_bin -Xms256m -Xmx256m -jar $vmaster_dir".
	"venn.jar --list $outfile --svg $svg_outfile";
}

print STDERR "VennMaster Command:\n$vm_cmd\n" if $verbose;

# LAUNCH THE VENN MASTER PROGRAM UNLESS SET OTHERWISE
unless ($no_launch_vm) {
    system ($vm_cmd) ||
	die "Could not run the VennMaster program with cmd:\n$vm_cmd\n";
}

#-----------------------------------------------------------+
# TO DO : CREATE COVERAGE DIAGRAM AS A HEAT MAP             |
#-----------------------------------------------------------+

# This can me made to be an option
# It would also be useful to do this as an exported file
# that could be viewed in Apollo or Gbrowse
# 

exit;

#-----------------------------------------------------------+
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+

sub ParseConfigFile {

#-----------------------------------------------------------+
# This will parase a configuration text file to set         |
# variables that would normally be set at the command line  |
# This is not the fastest way to parse the config file      |
# but this makes the subfuction reusable.                   |
# If the variable name occurs more then once in the text    |
# file, the last occurrence will be used.                   |
#-----------------------------------------------------------+
    my $ConfigFilePath = $_[0];
    my $VarName = $_[1];
    my $VarValue;


    open (CONFILE, $ConfigFilePath) ||
        die "Could not open config file:\n\t$ConfigFilePath";

    while (<CONFILE>) {
        chomp;                 # Remove newline character
        unless (m/\#.*/)       # Ignore comment lines
        {
            my @SplitLine = split;
            if ($SplitLine[0] =~ $VarName){
                $VarValue = $SplitLine[1];}
        }
    }
    close CONFILE;
    return $VarValue;

}

sub GetShortName {
#-----------------------------------------------------------+
# CONVERT THE  NAME TO A SHORTENED FORM                     |
#-----------------------------------------------------------+ 
# This will be used to write short names for the CrossTab
# table or any other use of the name in the program. The short
# form of the name will be limited to 8 characters
# NOTE:
# This is a kluge that will only work for the BLAST databases
# that I am working with in wheat (03/12/2007 - JCE)
# I should definited convert this to as hash table
# 
# THIS SHOULD BE DONE AS A HASH AND NOT IF ELSE STATEMETNS
# 10/22/2007 -JCE
# Al alternative is to use a db config file that cotains
# this information.
# Best alternative is probably a text file that can be parsed.

    my $InName = $_[0]; #  Input name string
    my $OutName;        #  Set scope for the name to return

    #-----------------------------------------------------------+
    # REPEAT DATABASES                                          |
    #-----------------------------------------------------------+
    if ($InName =~ m/TREP9_nr/) {
	$OutName = "TREP9nr";
    }
    elsif ($InName =~ m/TREP9_total/) {
	$OutName = "TREP9tot";
    }
    elsif ($InName =~ m/mips_REdat/) {
	$OutName = "MIPS";
    }
    elsif ($InName =~ m/TIGR_Gram/) {
	$OutName = "TIGRGram";
    }
    elsif ($InName =~ m/RB_pln/) {
	$OutName = "RB_pln";
    }
    elsif ($InName =~ m/TATrans/) {
	$OutName = "TATran";
    }
    elsif ($InName =~ m/Wessler/) {
	$OutName = "Wessler ";
    }
    elsif ($InName =~ m/SanMiguel/) {
	$OutName = "SanMig";
    #-----------------------------------------------------------+
    # EST DATABASES                                             |
    #-----------------------------------------------------------+
    # TIGR GENE INDICES	
    }
    elsif ($InName =~ m/TaGI_10/) {
	$OutName = "TaGI_10";    	
    }
    elsif ($InName =~ m/AtGI_13/) {
	$OutName = "AtGI_13";    	
    }
    elsif ($InName =~ m/ZmGI_17/) {
	$OutName = "ZmGI_17";
    }
    elsif ($InName =~ m/SbGI_8/) {
	$OutName = "SbGI_8";
    }
    elsif ($InName =~ m/OsGI_17/) {
	$OutName = "OsGI_17";
    }
    elsif ($InName =~ m/HvGI_9/) {
	$OutName = "HvGI_17";
    #-----------------------------+ 
    # EST: TIGR TAs               |
    #-----------------------------+
    }
    elsif ($InName =~ m/AcTA_1/) {
	$OutName = "AcTA_1";
    }
    elsif ($InName =~ m/AsTA_1/) {
	$OutName = "AsTA_1";
    }
    elsif ($InName =~ m/AsTA_1/) {
	$OutName = "AsTA_1";
    }
    elsif ($InName =~ m/AvenSatTA_2/) {
	$OutName = "AvenSaTa";
    }
    elsif ($InName =~ m/CdTA_2/) {
	$OutName = "CeTA_2";
    }
    elsif ($InName =~ m/FaTA_2/) {
	$OutName = "FaTA_2";
    }
    elsif ($InName =~ m/HvTA_2/) {
	$OutName = "HvTA_2";
    }
    elsif ($InName =~ m/OsTA_2/) {
	$OutName = "OsTA_2";
    }
    elsif ($InName =~ m/PgTA_2/) {
	$OutName = "PgTA_2";
    }
    elsif ($InName =~ m/SbTA_2/) {
	$OutName = "SbTA_2";
    }
    elsif ($InName =~ m/ShTA_2/) {
	$OutName = "ShTA_2";
    }
    elsif ($InName =~ m/SoTA_2/) {
	$OutName = "SoTA_2";
    }
    elsif ($InName =~ m/SpTA_2/) {
	$OutName = "SpTA_2";
    }
    elsif ($InName =~ m/TaTA_2/) {
	$OutName = "TaTA_2";
    }
    elsif ($InName =~ m/TmTA_2/) {
	$OutName = "TmTA_2";
    }
    elsif ($InName =~ m/TtTA_1/) {
	$OutName = "TtTA_1";
    }
    elsif ($InName =~ m/ZmB73TA_1/) {
	$OutName = "ZmB71TA_1";
    #-----------------------------+	
    # EST: NCBI UniGenes	  |
    #-----------------------------+
    }
    elsif ($InName =~ m/AtUniGene_56/) {
	$OutName = "AtUn_56";
    }
    elsif ($InName =~ m/HvUniGene_47/) {
	$OutName = "HvUn_47";
    }
    elsif ($InName =~ m/OsUniGene_65/) {
	$OutName = "OsUn_65";
    }
    elsif ($InName =~ m/SbUniGene_21/) {
	$OutName = "SbUn_21";
    }
    elsif ($InName =~ m/SoUniGene_8/) {
	$OutName = "SoUn_8";
    }
    elsif ($InName =~ m/TaUniGene_46/) {
	$OutName = "TaUn_46";
    }
    elsif ($InName =~ m/ZmUniGene_61/) {
	$OutName = "ZmUn_61";
    #-----------------------------+
    # EST: PLANT GDB PUTS         |
    #-----------------------------+
    }
    elsif ($InName =~ m/AsPUT_157/) {
	$OutName = "AsPUT157";
    }
    elsif ($InName =~ m/AtPUT_157/) {
	$OutName = "AtPUT157";
    }
    elsif ($InName =~ m/HvPUT_157/) {
	$OutName = "HvPUT157";
    }
    elsif ($InName =~ m/OsIndPUT_157/) {
	$OutName = "OsIndPUT";
    }
    elsif ($InName =~ m/OsJap_157/) {
	$OutName = "OsJapPUT";
    }
    elsif ($InName =~ m/OsPut_157/) {
	$OutName = "OsPUT157";
    }
    elsif ($InName =~ m/SbPUT_157/) {
	$OutName = "SbPUT157";
    }
    elsif ($InName =~ m/SpPUT_157/) {
	$OutName = "SpPUT157";
    }
    elsif ($InName =~ m/TaPUT_157/) {
	$OutName = "TaPUT157";
    }
    elsif ($InName =~ m/TmPUT_157/) {
	$OutName = "TmPUT157";
    }
    elsif ($InName =~ m/ZmPUT_157/) {
	$OutName = "ZmPUT157";
    #-----------------------------------------------------------+
    # PLANT PROTEIN DATABASES                                   |
    #-----------------------------------------------------------+
    #-----------------------------------------------------------+
    # FINISHED GENOMES PROTEINS BLASTX                          |
    #-----------------------------------------------------------+
    }
    elsif ($InName =~ m/TAIR6/) {
	$OutName = "TAIR6";
    }
    elsif ($InName =~ m/Os_5_cds_pep/) {
	$OutName = "OS_5_cds_pep";
    }
    elsif ($InName =~ m/Os_5_RAP1_pep/) {
	$OutName = "OS_5_cds_pep";
    #-----------------------------------------------------------+
    # UNIPROT                                                   |
    #-----------------------------------------------------------+
    }
    elsif ($InName =~ m/Os_5_RAP1_pep/) {
	$OutName = "OS_5_cds_pep";
    #-----------------------------------------------------------+
    # COMBINATION DATABASES                                     |
    #-----------------------------------------------------------+
    }
    elsif ($InName =~ m/EST_DB/) {
	$OutName = "EST_DB";
    }
    elsif ($InName =~ m/TE_DB/) {
	$OutName = "TE_DB";
    #-----------------------------------------------------------+
    # TRIANNOTATION FILES                                       |
    #-----------------------------------------------------------+
    }
    elsif ($InName =~ m/1trf/) {
	$OutName = "trf";
    }
    elsif ($InName =~ m/2gmHv/) {
	$OutName = "gmHv";
    }
    elsif ($InName =~ m/2gmOs/) {
	$OutName = "gmOs";
    }
    elsif ($InName =~ m/2gmTa/) {
	$OutName = "gmTa";
    }
    elsif ($InName =~ m/2gmZm/) {
	$OutName = "gmZm";
    }
    elsif ($InName =~ m/2gID/) {
	$OutName = "gID";
    }
    elsif ($InName =~ m/2eugOs/) {
	$OutName = "eugOs";
    }
    elsif ($InName =~ m/2fGh/) {
	$OutName = "2fGh";
    }
    elsif ($InName =~ m/genscan/) {
	$OutName = "gensc";
    #-----------------------------------------------------------+
    # UNKNOWN DATABASES                                         |
    #-----------------------------------------------------------+
    # Use the input name if no match found
    # another alternative is to flag this as unkonwn
    }
    else {

	#$OutName = "UNK";
	$OutName = $InName;
    }

    return $OutName;

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

=head1 NAME

vennseq.pl - Venn diagram overlap of sequence features

=head1 SYNOPSIS

=head2 Usage

    vennseq.pl -i InputFastaFile.fasta -o OutputFile -d FeatureFileDir
               -f tab

=head2 Required Options

    -i       # Fasta file being annotated
    -o       # Output file
    -d       # Feature file directory
    -f       # Format of the annotation data [gff]
           
=head1 DESCRIPTION

Creates a Venn Diagram comparing the overlap of sequence
features along a query sequence. The purpose is to allow
me to visualize the overlap of Repeat Databases or gene
models along a BAC.            

=head1 REQUIRED ARGUMENTS

=over 2

=item -i

Path to the sequence file in fasta format.

=item -o

Path to the output file in the Venn text file.

=item -d

Directory containing sequence features.

=back

=head1 OPTIONS

=over 2

=item -f,--format

The feature file format. This must be one of ( tab | blast | gff ).
Currently only tab delim blast hits are supported. This is set to
gff default.

=item -s,--svg-out

Path to a svg format vector file of the Venn diagram. This will produce
a svg file at this path, otherwise to svg file is created.

=item --venn-master

The directory containing the VennMaster program.

=item --no-vm

This will do the anlysis of the suquence features, and write the 
analysis to STDOUT without using VennMaster to visualize the results.

=item --java-path

The path to the java binary. By default, this program will use 'java'
and assume that the java path you want to use is defined in your PATH.

=item --seq-id

The identifier to use for the sequence id

=item --seq-len

This is an alternative to providing a sequence file to determine the sequence
length from. This allows you to define the lenght of the sequence that the
features are being mapped onto.

=item --debug

Run the program in debug mode.

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

=back

=head1 DIAGNOSTICS

Error messages generated by this program and possible solutions are listed
below.

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

This program does not require a configuration file and does not

=head2 Environment Variables

The following variables may be dfined in the user environment.

=over

=item * VMASTER_DIR

This is the directory where the VennMaster program is located
This option can also be defined wit the --venn-master option 
at the command line.

=item * VMASTER_JAVA_BIN

This is the path of the java binary that you want to use to
run the VennMaster program. If this is not specified in the user
environment, this will use 'java' as defined in the PATH variable.
This option can also be defined with the --java-path option
at the command line.

=item * VMASTER_JAVA_MEM

This is the amount of memory to allot to the VennMaster program
If this is not specified, the program will use 512 GB memory. This
option can also be defined with the --java-mem option at the command line.

=back

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * VennMaster

This progra requires the VennMaster program 
L<http://www.informatik.uni-ulm.de/ni/staff/HKestler/vennm/doc.html>.

=back

=head2 Required Perl Modules

This program does not require Perl modules beyond those modules included in
a normal installation of Perl.

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWGPAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over 2

=item * Limited Testing with VennMaster

I have only tested this will VennMaster on the Mac OS 10.5.

=item * Separate sources must be in separate files

For gff files, the separate sources must be in separate gff files. Currently
sources defined in the second colum will be regarded as the same source
if they are in the same gff file.

=back

=head1 SEE ALSO

The vennseq.pl program is part of the DAWGPAWS package of genome
annotation programs. See the DAWGPAWS web page 
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

Started: 03/06/2007

Updated: 04/06/2009

Version: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 03/06/2007
# - Program started
# - Included ParseConfigFile from previous work
# - Main body of program written
#   Reading input fasta file and sequence features and
#   loading to the FeatMatrix 2d array.
# - Printing output to screen for testing purposes
# - Ability to read tab delim blast output
#
# 03/10/2007
# - Working with test tabl delim blast output
# - Adding a coverage calculation and report
# - Adding quiet flag to turn of debug print information
# - Fixed problem parsing the tab delimited BLAST files:
#   I was using The regular expression: (m/\#.*/) to ignore
#   comment lines. However this will ignore any line that 
#   CONTAINS the # symbol. Changing the reg exp to
#   (m/^\#.*/) fixes the problem. This is a particular issue with
#   BLAST databases that have been formatted in the 
#   RepBase format. (ie EnSpm1_TD#EnSpm)
# - Added CrosTab function
#
# 03/12/2007
# - Added GetShortName subfunction
# - Got GetShortName working for Repeat databases and
#   the for all other BLAST databases used for the 
#   Wheat Annotation project 
# - Added the EmptyFiles array. This will prevent
#   wasting time doing counts for zero length files
# - Added print of the Empty files to the end of the array
#
# 03/14/2007
# - Working on the automated startup of Venn Master
#
# 05/21/2007
# - Added Pod Documentation
#
# 11/11/2007
# - Reworking code to accept gff format
# - Changing options to the long option format
#
# 12/17/2007
# - Updated POD documentation slightly
#
# 01/05/2009
# - Minor code clean up
# - Added test files to the svn repository
# - Added debug option and debug variable
#
# 02/18/2009
# - Additional code clean up
# - Updated POD documentation
#
# - TEST CMD IS:
#   ./vennseq.pl -i /home/jestill/projects/wheat_annotation/VennTest/HEX0075G21.fasta.masked -o /home/jestill/projects/wheat_annotation/VennTest/TestOutTwo.txt -d /home/jestill/projects/wheat_annotation/VennTest/ -f tab -s /home/jestill/projects/wheat_annotation/VennTest/HEX0075G21.test.svg


# vennseq.pl -i HEX2493A05.masked.fasta -d features/ -o test_out.venn
