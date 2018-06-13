#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ltrstruc2gff.pl - Convert ltrstruc rpt to gff format  |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/25/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Convert *rpt.txt output files from LTR_STRUC to          |
#  gff formatted annotation coordinates.                    |
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

# Required variables
my $indir;
my $outdir;
my $repdir;

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_copy = 0;
my $do_seq_data = 0;          # Create files in outdir with sequence data

my $program = "LTR_Struc";
my $parameter;
my $name_root;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"     => \$indir,
                    "o|outdir=s"    => \$outdir,
		    "r|results=s"   => \$repdir,
		    # ADDITIONAL OPTIONS
		    "param=s"       => \$parameter,
		    "program=s"     => \$program,
		    "c|copy"        => \$do_copy,
		    "q|quiet"       => \$quiet,
		    "s|seq-data"    => \$do_seq_data,
		    "verbose"       => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"         => \$show_usage,
		    "version"       => \$show_version,
		    "man"           => \$show_man,
		    "h|help"        => \$show_help,);

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
    print "\ncnv_ltrstruc2gff.pl:\nVersion: $VERSION\n\n";
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
if ( (!$indir) || (!$outdir) || (!$repdir) ) {
    print "ERROR: Input directory must be specified\n" if !$indir;
    print "ERROR: Output directory must be specified\n" if !$outdir;
    print "ERROR: Reports directory must be specified\n" if !$repdir;
    print "\n";
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

unless ($repdir =~ /\/$/ ) {
    $repdir = $repdir."/";
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
# Get the FASTA files from the|
# directory provided by the   |
# var $indir                  |
#-----------------------------+
opendir( DIR, $indir ) 
    || die "Can't open directory:\n$indir"; 
my @fasta_files = grep /\.fasta$|\.fa$/, readdir DIR ;
closedir( DIR );
my $num_fasta_files = @fasta_files;

if ($num_fasta_files == 0) {
    print "\a";
    print "No fasta files were found in the input direcotry:\n";
    print "$indir\n";
    print "Fasta file MUST have the fasta or fa extension to be".
	" recognized as fasta files\n";
    exit;
}

#-----------------------------+
# Get the LTR_STRUC report    |
# files                       |
#-----------------------------+
opendir(REPDIR,$repdir) 
    || die "Can not open results direcotry:\n$repdir\n";
my @report_files = grep /rprt\.txt$/, readdir REPDIR;
@report_files = sort(@report_files);
closedir (REPDIR); 
# Sort the array of reports

my $num_report_files = @report_files;

if ($verbose) {
    print STDERR "\n-----------------------------------------------\n";
    print STDERR " Report files to process: $num_report_files\n";
    print STDERR " Fasta files to process: $num_fasta_files\n";
    print STDERR "-----------------------------------------------\n\n";
}

my $fasta_file_num = 0;

for my $ind_fasta_file (@fasta_files) {
    
    my $ind_report_num=0;
    
    $fasta_file_num++;
    
    #-----------------------------+
    # GET ROOT FILE NAME          |
    #-----------------------------+
    if ($ind_fasta_file =~ m/(.*)\.masked\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_fasta_file =~ m/(.*)\.fasta$/ ) {	    
	$name_root = "$1";
    }  
    elsif ($ind_fasta_file =~ m/(.*)\.fa$/ ) {	    
	$name_root = "$1";
    } 
    else {
	$name_root = "UNDEFINED";
    }
    my $name_root_len = length($name_root);
    
    if ($verbose) {
	print STDERR "\n--------------------------------------------------\n";
	print STDERR "Processing $name_root: $fasta_file_num of".
	" $num_fasta_files\n";
	print STDERR "--------------------------------------------------\n";
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
    # CREATE LTR_STRUC OUTDIR     |
    #-----------------------------+
    # Dir to hold copies of the ltr_struc results
    my $ltrstruc_dir = $name_root_dir."ltr_struc/";
    unless (-e $ltrstruc_dir) {
	mkdir $ltrstruc_dir ||
	    die "Could not create ltr_struc out dir:\n$ltrstruc_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    my $gff_dir = $name_root_dir."gff/";
    unless (-e $gff_dir) {
	mkdir $gff_dir ||
	    die "Could not create gff out dir:\n$gff_dir\n";
    }

    my $fasta_file_path = $indir.$ind_fasta_file;
    my $gff_out_path = $gff_dir.$name_root."_ltrstruc.gff";
    
    #-----------------------------+
    # FIND REPORTS THAT MATCH     |
    #-----------------------------+
    for my $ind_report (@report_files) {

	# This is the id of the sequence the report is for
	my $ind_report_id = substr ($ind_report,0,$name_root_len);

	if ($name_root =~ $ind_report_id) {

	    print STDERR "\tReport: $ind_report\n" if $verbose;
	    $ind_report_num++;

	    my $report_file_path = $repdir.$ind_report;
  
	    if ($ind_report_num == 1) {
		# If first record start new gff file
		ltrstruc2gff ( $fasta_file_path, $report_file_path,
			       $gff_out_path, 0, $ind_report_num, 
			       $do_seq_data, $program, $parameter);
	    }
	    else {
		# If not first record append to existing gff file
		ltrstruc2gff ( $fasta_file_path, $report_file_path,
			       $gff_out_path, 1, $ind_report_num,
			       $do_seq_data, $program, $parameter);

	    }

	} # End of if report is for the fasta seq 
	
    } # End for for each individual report

    print STDERR "\tNum Reports: $ind_report_num\n" if $verbose;

    #-----------------------------+
    # COPY FILES TO ltr_struc DIR |
    #-----------------------------+
    if ($do_copy) {
	opendir(REPDIR,$repdir) 
	    || die "Can not open directory:\n$repdir\n";
	my @ls_files = grep /^$name_root/, readdir REPDIR;
	closedir (REPDIR); 
	
	for my $ind_ls_file (@ls_files) {
	    print STDERR "Copying file: \t$ind_ls_file\n" if $verbose;

	    my $ind_ls_file_path = $repdir.$ind_ls_file;
	    my $file_copy_path = $ltrstruc_dir.$ind_ls_file;
	    
	    copy ($ind_ls_file_path, $file_copy_path) ||
		print STDERR "Can not copy file:\n $ind_ls_file_path TO\n".
		" $file_copy_path\n";
	} # End of for each individual ltr_struc file
	
    }

    # Temp exit while working on the code
    #if ($fasta_file_num == 2) {exit;}

}


exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


sub ltrstruc2gff {
    
    # VARS PASSED TO THE SUBFUNCTION
    # $source is program name
    # $param is parameter name
    my ($fasta_in, $report_in, $gff_out, $gff_append, $ls_id_num, 
	$print_seq_data, $source, $param) = @_;

    # print_seq_data - create fasta files containing the sequences
    #                  for the prediced LTRS, PBS etc. [Boolean]

    # Append parameter name to the program name
    if ($param) {
	$source = $source.":".$param;
    }

    # FASTA RELATED VARS
    my $qry_seq;
    
    # Counters
    my $ltr5_blank_count = 0;

    # LTR STRUC VARS
    my $ls_score;        # Score assigned by LTR_STRUC
    my $ls_contig_len;   # Length of the source contig
    my $ls_orientation;  # Orientation 
    my $ls_retro_len;    # Overall length of the retrotransposon
    my $ls_longest_orf;  # Length of the longest ORF
    my $ls_rt_frame1;
    my $ls_rt_frame2;
    my $ls_rt_frame3;
    my $ls_5ltr_len;     # Length of the 5' LTR
    my $ls_3ltr_len;     # Length of the 3' LTR
    my $ls_ltr_homology; # Percent ID of LTRs
    my $ls_5dinuc;       # 5' dinucleotide sequence
    my $ls_3dinuc;       # 3' dinucleotide sequence
    my $ls_5tsr_seq;     # 5' Target site rep sequence
    my $ls_3tsr_seq;     # 3' Target site rep sequence
    my $ls_5flank_seq;   # 5' Flanking sequence
    my $ls_3flank_seq;   # 3' Flanking sequence
    my $ls_pbs_seq;      # Primer Binding Site Sequence
    my $ls_ppt_seq;      # Polypuring Tract sequence
    my $ls_5id_seq;
    my $ls_3id_seq;
    my $ls_5ltr_seq;     # Sequence of the 5' LTR
    my $ls_3ltr_seq;
    my $ls_full_retro_seq; 
    my $par_5ltr_len;      # Length of the 5' LTR as parsed
    my $par_3ltr_len;      # Length of the 3' LTR as parsed

    # BOOLEANS
    my $in_rt_frame_1 = 0;
    my $in_rt_frame_2 = 0;
    my $in_rt_frame_3 = 0;
    my $in_5id_seq = 0;
    my $in_3id_seq = 0;
    my $in_5ltr = 0;
    my $in_3ltr = 0;
    my $in_complete_seq = 0;
    my $in_aligned_ltrs = 0;

    # Coordinate values
    my $full_retro_start;
    my $full_retro_end;
    my $ltr5_start;
    my $ltr5_end;
    my $ltr3_start;
    my $ltr3_end;
    my $pbs_start;
    my $pbs_end;
    my $ppt_start;
    my $ppt_end;

    # GFF Coordinates values
    # These are from absolute start of the query sequence string
    # starting the index value at one
    my $gff_full_retro_start;
    my $gff_full_retro_end;
    my $gff_ltr5_start;
    my $gff_ltr5_end;
    my $gff_ltr3_start;
    my $gff_ltr3_end;
    my $gff_pbs_start;
    my $gff_pbs_end;
    my $gff_ppt_start;
    my $gff_ppt_end;
    my $gff_5tsr_start;
    my $gff_5tsr_end;
    my $gff_3tsr_start;
    my $gff_3tsr_end;

    # Coordinate substring tests
    my $ppt_from_full_retro;
    my $ppt_from_query_seq;
    my $pbs_from_full_retro;
    my $pbs_from_query_seq;
    my $tsd5_from_query_seq;
    my $tsd3_from_query_seq;

    #-----------------------------+
    # OPEN GFF OUTPUT FILE        |
    #-----------------------------+
    if ($gff_append) {
	open (GFFOUT, ">>$gff_out") ||
	    die "ERROR: Could not open gff output file:\b$gff_out\n";
    }
    else 
    {
	open (GFFOUT, ">$gff_out") ||
	    die "ERROR: Could not open gff output file:\n$gff_out\n";
    }
    
    #-----------------------------+
    # LOAD FASTA SEQUENCE         |
    #-----------------------------+
    open (INFASTA, "$fasta_in") ||
	die "Can not open fasta input file:\n$fasta_in\n";
    while (<INFASTA>) {
	chomp;
	unless(m/^\>/) {
	    # Do uppercase to make sure lowercase masked sequences
	    # will properly match
	    $qry_seq = $qry_seq.uc($_);
	}
    }
    close (INFASTA);
    
    #-----------------------------+
    # GET DATA FROM REPORT FILE   |
    #-----------------------------+
    open (REPIN, "<$report_in")
	|| die "ERROR: Can not open report file:\n$report_in\n";

    while (<REPIN>) {
	chomp;
	if ($in_rt_frame_1) {

	}
	elsif (m/COMPLETE SEQUENCE OF PUTATIVE TRANSPOSON/) {
	    $in_3ltr = 0;
	    $in_complete_seq = 1;
	}
	elsif (m/^ALIGNED LTRS:/) {
	    $in_complete_seq = 0;
	    $in_aligned_ltrs = 1;
	}
	elsif ($in_rt_frame_2) {

	}
	elsif ($in_rt_frame_3) {

	}
	elsif ($in_5id_seq) {

	}
	elsif ($in_3id_seq) {

	}
	elsif ($in_complete_seq) {
	    $ls_full_retro_seq = $ls_full_retro_seq.$_;
	}
	elsif(/CUT-OFF SCORE:\s+(\d\.\d+)/){
	    $ls_score = $1;
	}
	elsif(m/TRANSPOSON IS IN (.*) ORIENTATION/) {
	    $ls_orientation = $1;
	    if ($ls_orientation =~ "NEGATIVE") {
		$ls_orientation = "-";
	    }
	    elsif ($ls_orientation =~ "POSITIVE") {
		$ls_orientation = "+";
	    }
	    else {
		# If return can not be parsed just use dot
		# this indicates unknown orientation if gff
		$ls_orientation = "."; 
	    }
	}
	#-----------------------------+
	# SEQUENCE DATA               |
	#-----------------------------+
	elsif ($in_5ltr) {
	    $ls_5ltr_seq = $ls_5ltr_seq.$_;
	    my $len_inline = length ($_);

	    # The following for debug
	    #print "\tLEN: $len_inline\n";

	    if ($len_inline == 0) {
		$ltr5_blank_count++;
		if ($ltr5_blank_count == 2) {
		    # Set in_5ltr boolean to false
		    $in_5ltr = 0;
		    $in_3ltr = 1;
		}
	    }

	}
	elsif ($in_3ltr) {
	    $ls_3ltr_seq = $ls_3ltr_seq.$_;
	}
	elsif (m/DINUCLEOTIDES: (.*)\/(.*)/) {
	    $ls_5dinuc = $1;
	    $ls_3dinuc = $2;
	}
	elsif (m/DIRECT REPEATS: (.*)\/(.*)/) {
	    $ls_5tsr_seq = $1;
	    $ls_3tsr_seq = $2;
	}
	elsif (m/PBS: (.*)/) {
	    $ls_pbs_seq = $1;
	}
	elsif (m/POLYPURINE TRACT: (.*)/){
	    $ls_ppt_seq = $1;
	}
	elsif (m/5\' FLANK: (.*)/) {
	    $ls_5flank_seq = $1;
	}
	elsif (m/3\' FLANK: (.*)/) {
	    $ls_3flank_seq = $1;
	}
	#-----------------------------+
	# OTHER DATA                  |
	#-----------------------------+
	elsif(m/OVERALL LENGTH OF TRANSPOSON:\s+(\d+)/){
	    $ls_retro_len = $1;
	}
	elsif(m/LENGTH OF LONGEST ORF:\s+(\d+)/){
	    $ls_longest_orf = $1;
	}
	elsif(m/LENGTH OF PUTATIVE 3\' LTR:\s+(\d+)/){
	    $ls_3ltr_len = $1;
	}
	elsif(m/LENGTH OF PUTATIVE 5\' LTR:\s+(\d+)/){
	    $ls_5ltr_len = $1;
	}
	elsif(m/LTR PAIR HOMOLOGY:\s+(\S+)%/){
	    $ls_ltr_homology = $1;
	}
	#-----------------------------+
	# SET BOOLEAN FLAGS           |
	#-----------------------------+
	elsif (m/LTRS:/) {
	    $in_5ltr = 1;
	}
    }

    close (REPIN);

    #-----------------------------+
    # GET COORDINATES             |
    #-----------------------------+
    $par_5ltr_len = length($ls_5ltr_seq);
    $par_3ltr_len = length($ls_3ltr_seq);

    $full_retro_start = index($qry_seq,$ls_full_retro_seq) + 1;
    $full_retro_end = $full_retro_start + $ls_retro_len;
    
    # The following will have a problem on 100% identical LTRs
    # however, telling the search to start at the end of the
    # 5' LTR will solve this problem since the index function
    # will accept an offset at the third argument
    # Index returns the first occurrence of the string
    $ltr5_start = index ($ls_full_retro_seq, $ls_5ltr_seq) + 1;
    $ltr5_end = $ltr5_start + $ls_5ltr_len;
    $ltr3_start = index ($ls_full_retro_seq, $ls_3ltr_seq) + 1;
    $ltr3_end = $ltr3_start + $ls_3ltr_len;
    $pbs_start = index ($ls_full_retro_seq, $ls_pbs_seq) + 1 ;
    $pbs_end = $pbs_start + length($ls_pbs_seq);
    $ppt_start = index ($ls_full_retro_seq, $ls_ppt_seq) + 1;
    $ppt_end = $ppt_start + length($ls_ppt_seq);

    #-----------------------------+
    # GET EXTRACTED SEQS          |
    #-----------------------------+
    # This is to test if the coordinates I am returning matches the
    # observations that LTR_STRUC is reporting
    $ppt_from_full_retro = substr ($ls_full_retro_seq, $ppt_start - 1,
				   length($ls_ppt_seq) );
    $pbs_from_full_retro = substr ($ls_full_retro_seq, $pbs_start - 1,
				   length($ls_pbs_seq) );

    # Note that the coordinates to get the correct substring below
    # start the  string index from 0 so to get this in gff I would need
    # to set 
    # ppt_abs_start to $ppt_start - 1 + $full_retro_start
    $ppt_from_query_seq = substr ( $qry_seq, 
				   $ppt_start - 2 + $full_retro_start, 
				   length($ls_ppt_seq) );
    $pbs_from_query_seq = substr ( $qry_seq, 
				   $pbs_start - 2 + $full_retro_start, 
				   length($ls_pbs_seq) );


    $tsd5_from_query_seq = substr ( $qry_seq, 
				    $full_retro_start - 
				    length($ls_5tsr_seq) - 1,
				    length($ls_5tsr_seq) );

    $tsd3_from_query_seq = substr ( $qry_seq, 
				    $full_retro_start + 
				    $ls_retro_len - 1,
				    length($ls_3tsr_seq) );
    
    #-----------------------------+
    # GFF COORDINATES             |
    #-----------------------------+
    # Full retro (not includings tsds)
    $gff_full_retro_start = $full_retro_start;                      # OK
    $gff_full_retro_end = $gff_full_retro_start + $ls_retro_len - 1;# OK
    # Subunits (including the putative tsds)
    $gff_pbs_start = $pbs_start - 1 + $full_retro_start;            # OK
    $gff_pbs_end = $gff_pbs_start + length($ls_pbs_seq) - 1;        # OK
    $gff_ppt_start = $ppt_start - 1 + $full_retro_start;            # OK
    $gff_ppt_end = $gff_ppt_start + length($ls_ppt_seq) - 1;        # OK
    # LTRs
    $gff_ltr5_start = $gff_full_retro_start;                        # OK
    $gff_ltr5_end = $gff_ltr5_start + $ls_5ltr_len - 1;             # OK
    $gff_ltr3_end = $gff_full_retro_end;                            # OK
    $gff_ltr3_start = $gff_ltr3_end - $ls_3ltr_len + 1;             # OK
    # TSRs - Currently returning correct on positive strand
    $gff_5tsr_start = $gff_full_retro_start - length($ls_5tsr_seq); # OK
    $gff_5tsr_end = $gff_5tsr_start + length($ls_5tsr_seq) - 1;     # OK
    $gff_3tsr_start = $gff_full_retro_start + $ls_retro_len;        # OK
    $gff_3tsr_end = $gff_3tsr_start + length($ls_3tsr_seq) - 1;     # OK

    #-----------------------------+
    # PRINT GFF OUTPUT            |
    #-----------------------------+
    # Data type follows SONG
    # http://song.cvs.sourceforge.net/*checkout*/song/ontology/so.obo
    my $gff_result_id = "ltr_struc_".$ls_id_num;
    print GFFOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"LTR_retrotransposon\t".   # Data type, has to be exon for APOLLO
	"$gff_full_retro_start\t". # Start
	"$gff_full_retro_end\t".   # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print GFFOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"primer_binding_site\t".   # SO:0005850
	"$gff_pbs_start\t".        # Start
	"$gff_pbs_end\t".          # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print GFFOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"RR_tract\t".              # SO:0000435 
	"$gff_ppt_start\t".        # Start
	"$gff_ppt_end\t".          # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print GFFOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"five_prime_LTR\t".        # SO:0000425
	"$gff_ltr5_start\t".       # Start
	"$gff_ltr5_end\t".         # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print GFFOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"three_prime_LTR\t".       # SO:0000426  
	"$gff_ltr3_start\t".       # Start
	"$gff_ltr3_end\t".         # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print GFFOUT "$name_root\t".     # Seq ID
	"$source\t".                 # Source
	"target_site_duplication\t". # SO:0000434
	"$gff_5tsr_start\t".         # Start
	"$gff_5tsr_end\t".           # End
	"$ls_score\t".               # Score
	"$ls_orientation\t".         # Strand
	".\t".                       # Frame
	"$gff_result_id\n";          # Retro Id

    print GFFOUT "$name_root\t".     # Seq ID
	"$source\t".                 # Source
	"target_site_duplication\t". # SO:0000434
	"$gff_3tsr_start\t".         # Start
	"$gff_3tsr_end\t".           # End
	"$ls_score\t".               # Score
	"$ls_orientation\t".         # Strand
	".\t".                       # Frame
	"$gff_result_id\n";          # Retro Id


    #-----------------------------+
    # PRINT TO STDOUT             |
    #-----------------------------+
    # This allows for an easy way to send the results to a single file
    # that contains all of the LTR_Struc results

    print STDOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"LTR_retrotransposon\t".   # Data type, has to be exon for APOLLO
	"$gff_full_retro_start\t". # Start
	"$gff_full_retro_end\t".   # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print STDOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"primer_binding_site\t".   # SO:0005850
	"$gff_pbs_start\t".        # Start
	"$gff_pbs_end\t".          # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print STDOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"RR_tract\t".              # SO:0000435 
	"$gff_ppt_start\t".        # Start
	"$gff_ppt_end\t".          # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print STDOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"five_prime_LTR\t".        # SO:0000425
	"$gff_ltr5_start\t".       # Start
	"$gff_ltr5_end\t".         # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print STDOUT "$name_root\t".   # Seq ID
	"$source\t".               # Source
	"three_prime_LTR\t".       # SO:0000426  
	"$gff_ltr3_start\t".       # Start
	"$gff_ltr3_end\t".         # End
	"$ls_score\t".             # Score
	"$ls_orientation\t".       # Strand
	".\t".                     # Frame
	"$gff_result_id\n";        # Retro Id

    print STDOUT "$name_root\t".     # Seq ID
	"$source\t".                 # Source
	"target_site_duplication\t". # SO:0000434
	"$gff_5tsr_start\t".         # Start
	"$gff_5tsr_end\t".           # End
	"$ls_score\t".               # Score
	"$ls_orientation\t".         # Strand
	".\t".                       # Frame
	"$gff_result_id\n";          # Retro Id

    print STDOUT "$name_root\t".     # Seq ID
	"$source\t".               # Source
	"target_site_duplication\t". # SO:0000434
	"$gff_3tsr_start\t".         # Start
	"$gff_3tsr_end\t".           # End
	"$ls_score\t".               # Score
	"$ls_orientation\t".         # Strand
	".\t".                       # Frame
	"$gff_result_id\n";          # Retro Id


    #-----------------------------+
    # PRINT SUMMARY OUTPUT        |
    #-----------------------------+
    if ($verbose) {
	print STDERR "\t\tScore: $ls_score\n";
	print STDERR "\t\tOrientation: $ls_orientation\n";
	print STDERR "\t\tRetro Len: $ls_retro_len\n";
	print STDERR "\t\tLongest Orf: $ls_longest_orf\n";
	print STDERR "\t\tPair Homology: $ls_ltr_homology\n";
	print STDERR "\t\t5LTR Len LST: $ls_5ltr_len\n";
	print STDERR "\t\t5LTR Lem PAR: $par_5ltr_len\n";
	print STDERR "\t\t3LTR Len LST: $ls_3ltr_len\n";
	print STDERR "\t\t3LTR Len PAR: $par_3ltr_len\n";
	
	print STDERR "\n\t\tEXTRACTED COORDINATES:\n";
	print STDERR "\t\t5\'LTR Start: $ltr5_start\n";
	print STDERR "\t\t5\'LTR End: $ltr5_end\n";
	print STDERR "\t\t3\'LTR Start: $ltr3_start\n";
	print STDERR "\t\t3\'LTR End: $ltr3_end\n";
	print STDERR "\t\tPBS Start: $pbs_start\n";
	print STDERR "\t\tPBS end : $pbs_end\n";
	print STDERR "\t\tPPT Start: $ppt_start\n";
	print STDERR "\t\tPPT End: $ppt_end\n";
	
	print STDERR "\t\tRetro Start: $full_retro_start\n";
	print STDERR "\t\tRetro End: $full_retro_end\n";
	
	# SEQUENCE DATA
	print STDERR "\n\t\tSEQUENCE DATA:\n";
	print STDERR "\t\t5Dinuc: $ls_5dinuc\n";
	print STDERR "\t\t3Dinuc: $ls_3dinuc\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\t5\' TSD TEST:\n";
	print STDERR "\t\tLTR_STRUC: $ls_5tsr_seq\n";
	print STDERR "\t\tFULL SEQ:  $tsd5_from_query_seq\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\t3\' TSD TEST:\n";
	print STDERR "\t\tLTR_STRUC: $ls_3tsr_seq\n";
	print STDERR "\t\tFULL SEQ:  $tsd3_from_query_seq\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\tPBS STRING TEST:\n";
	print STDERR "\t\tLTR_STRUC:  $ls_pbs_seq\n";
	print STDERR "\t\tFROM RETRO: $pbs_from_full_retro\n";
	print STDERR "\t\tFULL SEQ:   $pbs_from_query_seq\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n\t\t=====================\n";
	print STDERR "\t\tPPT STRING TEST:\n";
	print STDERR "\t\tLTR_STRUC:  $ls_ppt_seq\n";
	print STDERR "\t\tFROM RETRO: $ppt_from_full_retro\n";
	print STDERR "\t\tFULL SEQ:   $ppt_from_query_seq\n";
	print STDERR "\t\t=====================\n";
	
	print STDERR "\n";
	print STDERR "\t\t5\'Flank: $ls_5flank_seq\n";
	print STDERR "\t\t3\'Flank: $ls_3flank_seq\n";
	#print "\t\t5\'LTR: $ls_5ltr_seq\n";
	#print "\t\t3\'LTR: $ls_3ltr_seq\n";
	#print "$ls_full_retro_seq\n";
    }

    #-----------------------------+
    # PRINT SEQ DATA              |
    #-----------------------------+
    if ($print_seq_data) {
	
	# This uses the global $outdir variable

	#-----------------------------+
	# 5' LTR                      |
	#-----------------------------+
	my $ltr5_seq_path = $outdir."ltr5_ltr_struc.fasta";
	open (LTR5OUT, ">>$ltr5_seq_path") ||
	    die "Can not open 5\'LTR sequence output file:\n$ltr5_seq_path\n";
	print LTR5OUT ">".$name_root."_".$gff_result_id.
	    "|five_prime_ltr\n";
	print LTR5OUT "$ls_5ltr_seq\n";
	close (LTR5OUT);
	
	#-----------------------------+
	# 3' LTR                      |
	#-----------------------------+
	my $ltr3_seq_path = $outdir."ltr3_ltr_struc.fasta";
	open (LTR3OUT, ">>$ltr3_seq_path") ||
	    die "Can not open 3\'LTR sequence output file:\n$ltr3_seq_path\n";
	print LTR3OUT ">".$name_root."_".$gff_result_id.
	    "|three_prime_ltr\n";
	print LTR3OUT "$ls_3ltr_seq\n";
	close (LTR3OUT);

	#-----------------------------+
	# PBS                         |
	#-----------------------------+
	my $pbs_seq_path = $outdir."pbs_ltr_struc.fasta";
	open (PBSOUT, ">>$pbs_seq_path") ||
	    die "Can not open PBS sequence output file:\n$pbs_seq_path\n";
	print PBSOUT ">".$name_root."_".$gff_result_id.
	    "|primer_binding_site\n";
	print PBSOUT "$ls_pbs_seq\n";
	close (PBSOUT);

	#-----------------------------+
	# PPT                         |
	#-----------------------------+
	my $ppt_seq_path = $outdir."ppt_ltr_struc.fasta";
	open (PBSOUT, ">>$ppt_seq_path") ||
	    die "Can not open PPT sequence output file:\n$ppt_seq_path\n";
	print PBSOUT ">".$name_root."_".$gff_result_id.
	    "|RR_tract\n";
	print PBSOUT "$ls_pbs_seq\n";
	close (PBSOUT);
	
	#-----------------------------+
	# FULL RETRO MODEL            |
	#-----------------------------+ 
	# NOTE THIS INCLUDES NESTED LTR RETROS
	# LTR_retrotransposon
	my $ltr_seq_out = $outdir."full_ltr_retro.fasta";
	open (LTROUT, ">>$ltr_seq_out") ||
	    die "Can not open full ltr retro sequence outfile\n";
	print LTROUT ">".$name_root."_".$gff_result_id.
	    "|LTR_retrotransposon\n";
	print LTROUT "$ls_full_retro_seq\n";
	close (LTROUT);
	    

    }

    close (GFFOUT);

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

cnv_ltrstruc2gff.pl - Convert LTR_STRUC report output files to gff

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_ltrstruc2gff.pl -i InDir -o OutDir -r LStrucOut

=head2 Required Arguments

    --indir         # Directory with the fasta files
                    # This is used to find the root seq names
    --outdir        # Directory for the base output dir
    --results       # Directory containing the LTR_STRUC results

=head1 DESCRIPTION

Given a directory containing the output from LTR_STRUC and a directory
containing the fasta files that the structural predictions were made
from,

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path of the intput directory containing the fasta files that were
analyzed by LTR_STRUC. It may seem awkward to need to provide this
directory of fasta files, but it is currently necessary to be able
to find the reports for that sequence set. LTR_Struc does not provide the
sequence name in a format that can be parsed. Furthermore, LTR_struc
does not provide the location of the sequence features on the sequence file,
so the fasta file is needed to map the LTR retrotransposon model back
onto the sequence file.

=item -o,--outdir

Path of the output directory that will serve as the base for the
output from the conversion to gff. Every sequence will have a subdirectory
created in this parent dir.

=item -r,--results

Path of the directory containing the results from LTR_STRUC. It is 
expected that these file names will end with rprt.txt.

=back

=head1 OPTIONS

=over 2

=item --program

Specify the program name to use in the GFF output file. By default, the program
name used is ltr_sturc

=item --param

Specify the parameter set name to used. This will be appended to the 
program name in the source column of the gff output file.

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

=back

=head1 CONFIGURATION AND ENVIRONMENT

The program cnv_ltrstruc2gff.pl does not currently require an 
external configuration file or make use of variables defined in the
user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over

=item * LTR_Struc

This program is designed to work with output files generated by
the LTR_Struc
program. This program is available for download from:
http://www.genetics.uga.edu/retrolab/data/LTR_Struc.html

=back

=head2 Required Perl Modules

=over

=item * File::Copy

This module is required to copy the results.

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

=head1 SEE ALSO

The cnv_ltrstruc2gff.pl program is part of the DAWG-PAWS package of genome
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

STARTED: 09/25/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
# 09/25/2007
# - Basic input output started 
# 09/26/2007
# - ltrstruc2gff subfunction added
# 09/27/2007
# - Fixed coordinates in the gff output
# - Added code to copy ltr_struc related files to the
#   outdir.
# 09/28/2007
# - Working with gff files produced and wheat.tiers file to
#   to display the LTR-Retro subcomponent. Not really working out.
#   May have to define features of the ltr_retrospans
#   and do annotation of the LTR retromodels in external process
# - Added code to generate fasta files for the following components
#    - Full Retro Sequence
#    - 5' LTR
#    - 3' LTR
#    - PBS
#    - PPT
# 12/12/2007
# - Changed print_help subfunction to extract help and 
#   usage messages from the POD documentation
# - Updated POD documentation
#
# 03/30/2009
# - Added program
# - Added parameter
# - Updated POD
# - Added support for masked.fasta in seq name
# - Added support for the use of lowercase masked fasta files
#   by making the input of the Fasta file use uc()
