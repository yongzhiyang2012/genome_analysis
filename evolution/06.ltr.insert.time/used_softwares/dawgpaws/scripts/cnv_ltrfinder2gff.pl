#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_ltrfinder2gff.pl - Converts ltr_finder to gff         |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 09/14/2007                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Converts the LTR_FINDER results to gff format.           |
#                                                           |
# VERSION: $Rev: 761 $                                      |
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
use Cwd;                       # Get the current working directory
use File::Copy;                # Copy files

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $infile;                    # Infile. textfile result from LTR_FINDER
my $outfile;                   # Outfile.

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $append = 0;                # Append gff output to existing file

my $param;                    # Suffix appended to the end of the gff 
my $seqname;          #
my $program = "ltr_finder";   # The source program

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    # ADDITIONAL OPTIONS
		    "p|param=s"   => \$param,
		    "program=s"   => \$program,
		    "s|seqname=s" => \$seqname,
		    "append"      => \$append,
		    "q|quiet"     => \$quiet,
		    "verbose"     => \$verbose,
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

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\ncnv_ltrfinder2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+

ltrfinder2gff ($program, $seqname, $infile, $outfile, $append, $param);

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

sub ltrfinder2gff {
    
    # seq_id could be extracted from the ltrfinder result
    # but passing it to the subfunction directly allows for cases
    # where the id assigned by ltr_struc differs from the way
    # the user is referring to the assembly
    # MAY WANT TO ALLOW FOR USING THE $ls_seq_id
    my ($gff_src, $seq_id, $lf_infile, $gffout, $do_append, $gff_suffix) = @_;
    
    # The gff src id
    #my $gff_src = "ltr_finder";
    if ($gff_suffix) {
	$gff_src = $gff_src.":".$gff_suffix;
    }

    my $print_gff_out = 0;  # Boolean to print out gff data

    my $gff_str_out;        # A single out string line of gff out

    # LTF FINDER COUTNERS/INDICES
    my $lf_id_num = 0;       # Incremented for every predicted model

    # LTR_FINDER BOOLEANS
    my $in_emp_details = 0;  # Exact match pairs details
    my $in_ltr_align = 0;    # Details of the LTR alignment
    my $in_pbs_align = 0;
    my $in_ppt;

    
    #
    my $lf_prog_name;             # LTR Finder program name
    my $lf_seq_id;                # Seq id
    my $lf_seq_len;               # Sequence length
    my $lf_version;               # Version of LTR finder being parsed
    my $lf_trna_db;               # tRNA database used
    
    my $lf_strand;                # Strand of the LTR
    my $lf_span_start;            # Location start
    my $lf_span_end;              # Location end
    my $lf_length;                # Length
    
    my $lf_score;                 # Score
    my $lf_ltr_similarity;        # Similarity of the LTRs
    
    # Status strings
    my $has_5ltr_tg;                 # TG in 5' END of 5' LTR
    my $has_5ltr_ca;                 # CA in 3' END of 5' LTR
    my $has_3ltr_tg;                 # TG in 5' END of 3' LTR
    my $has_3ltr_ca;                 # CA in 3' END of 3' LTR
    my $has_tsr;                     # Has Target Site Replication
    my $has_pbs;                     # Has Primer Binding Site
    my $has_ppt;                     # Has Poly Purine Tract
    my $has_rt;                      # Has Reverse Transcriptase
    my $has_in_core;                 # Has Integrase Core
    my $has_in_cterm;                # Has Integrase C-term
    my $has_rh;                      # Has RNAseH
    
    my $lf_ltr_id;    # Id number assigned to the LTR retrotransposon

    # LTR COORDINATES
    my $lf_5ltr_start;
    my $lf_5ltr_end;
    my $lf_5ltr_len;
    my $lf_3ltr_start;
    my $lf_3ltr_end;
    my $lf_3ltr_len;

    # TSR COORDINATES
    my $lf_5tsr_start;             # Start of the 5' TSR
    my $lf_5tsr_end;               # End of the 5' TSR
    my $lf_3tsr_start;             # Start of the 3' TSR
    my $lf_3tsr_end;               # End of the 3' TSR 
    my $lf_tsr_string;             # String of bases in the TSR

    # BOUNDARY SHARPNESS
    my $lf_sharp_5;                # Sharpness of 5' Boundary
    my $lf_sharp_3;                # Sharpness of 3' Boundary
    
    # PBS
    my $lf_pbs_num_match;          # Number of matched bases in the PBS
    my $lf_pbs_aln_len;            # PBS alignment length
    my $lf_pbs_start;              # Start of the PBS signal
    my $lf_pbs_end;                # End of the PBS signal
    my $lf_pbs_trna;               # PBS tRNA type and anti-codon

    # PPT 
    my $lf_ppt_num_match;          # Number of matched based in the PPT
    my $lf_ppt_aln_len;            # PPT alignment length
    my $lf_ppt_start;              # Start of the PPT
    my $lf_ppt_end;                # End of the PPT
    
    #-----------------------------+
    # DOMAIN DATA                 |
    #-----------------------------+
    
    # GENERAL DOMAIN VARS
    my $lf_domain_dom_start;
    my $lf_domain_dom_end;
    my $lf_domain_orf_start;
    my $lf_domain_orf_end;
    my $lf_domain_name;            # Type of the domain
    
    # INTEGRASE CORE
    #my $has_in_core = 0;           # Boolean for has integrase core
    my $lf_in_core_dom_start;
    my $lf_in_core_dom_end;
    my $lf_in_core_orf_start;
    my $lf_in_core_orf_end;
    
    # INTEGRASE C-TERM
    #my $has_in_cterm = 0;
    my $lf_in_cterm_dom_start;
    my $lf_in_cterm_dom_end;
    my $lf_in_cterm_orf_start;
    my $lf_in_cterm_orf_end;

    # RNASEH
    #my $has_rh = 0;
    my $lf_rh_dom_start;
    my $lf_rh_dom_end;
    my $lf_rh_orf_start;
    my $lf_rh_orf_end;
    
    # RT
    #my $has_rt = 0;
    my $lf_rt_dom_start;
    my $lf_rt_dom_end;
    my $lf_rt_orf_start;
    my $lf_rt_orf_end;
    

    #-----------------------------+
    # OPEN GFF OUTFILE            |
    #-----------------------------+
    # Default to STDOUT if no arguemtn given
    if ($gffout) {
	if ($do_append) {
	    open (GFFOUT, ">>$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
	else {
	    open (GFFOUT,">$gffout") ||
		die "ERROR: Can not open gff outfile:\n $gffout\n";
	}
    } 
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }

    #-----------------------------+
    # OPEN INPUT FILE             |
    #-----------------------------+
    if ($lf_infile) {
	open (INFILE, "<$lf_infile") ||
	    die "ERROR: Can not open LTR_FINDER result file\n $lf_infile\n";
	
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (INFILE, "<&STDIN") ||
	    die "Can not accept input from standard input.\n";
    }

    while (<INFILE>) {
	chomp;
        #    print $_."\n";
	
	
	# CHECK BOOLEANS
	
	
	# 
	if (m/Nothing Header(.*)/) {
	    
	}
	

	# IN NEW REC, GET ID
	elsif (m/^\[(.*)\]/) {

	    #-----------------------------+
	    # PRINT STORED GFF OUTPUT     |
	    # IF VALUES ARE PRESENT       |
	    #-----------------------------+
	    
	    # override seq id if one passed 
	    if ($seqname) {
		$lf_seq_id = $seqname;
	    }
	    
	    # FULL SPAN
	    if ($lf_span_start) {
		$gff_str_out = "$lf_seq_id\t". # Seq ID
		    "$gff_src\t".                # Source
		    "LTR_retrotransposon\t".     # Data type
		    "$lf_span_start\t".          # Start
		    "$lf_span_end\t".            # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
	    }
	    
	    if ($lf_5ltr_start) {
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "five_prime_LTR\t".         # Data type
		    "$lf_5ltr_start\t".         # Start
		    "$lf_5ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
	    }
	    
	    if ($lf_3ltr_start) {
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "three_prime_LTR\t".        # Data type
		    "$lf_3ltr_start\t".         # Start
		    "$lf_3ltr_end\t".           # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
	    }
	    
	    if ($has_pbs) {
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "primer_binding_site\t".    # Data type
		    "$lf_pbs_start\t" .         # Start
		    "$lf_pbs_end\t".            # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
	    }
	    
	    
	    if ($has_ppt) {
		$gff_str_out = "$lf_seq_id\t".  # Seq ID
		    "$gff_src\t".               # Source
		    "RR_tract\t".               # Data type
		    "$lf_ppt_start\t".          # Start
		    "$lf_ppt_end\t".            # End
		    "$lf_score\t".              # Score
		    "$lf_strand\t".             # Strand
		    ".\t".                      # Frame
		    "ltr_finder_$lf_ltr_id\n";  # Retro ID
		print GFFOUT $gff_str_out;
	    }
	    
	    
	    if ($has_tsr) {
		
		$gff_str_out = "$lf_seq_id\t".   # Seq ID
		    "$gff_src\t".                # Source
		    "target_site_duplication\t". # Data type
		    "$lf_5tsr_start\t".          # Start
		    "$lf_5tsr_end\t".            # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;	    
		
		$gff_str_out = "$lf_seq_id\t".   # Seq ID
		    "$gff_src\t".                # Source
		    "target_site_duplication\t". # Data type
		    "$lf_3tsr_start\t".          # Start
		    "$lf_3tsr_end\t".            # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
	    }
	    
	    
	    # Integrase Core
	    if ($has_in_core) {
		#/////////
		# NOT SONG
		#\\\\\\\\\
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".              # Source
		    "integrase_core_domain\t".  # Data type
		    "$lf_in_core_dom_start\t".   # Start
		    "$lf_in_core_dom_end\t".     # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
		#/////////
		# NOT SONG
		#\\\\\\\\\
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".              # Source
		    "integrase_core_orf\t".      # Data type
		    "$lf_in_core_orf_start\t".   # Start
		    "$lf_in_core_orf_end\t".     # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
	    } # End of has in_core
	    
	    
	    if ($has_in_cterm) {
		
		#/////////
		# NOT SONG
		#\\\\\\\\\
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".                # Source
		    "integrase_cterm_domain\t".  # Data type
		    "$lf_in_cterm_dom_start\t".  # Start
		    "$lf_in_cterm_dom_end\t".    # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".                # Source
		    "integrase_cterm_orf\t".     # Data type
		    "$lf_in_cterm_orf_start\t".  # Start
		    "$lf_in_cterm_orf_end\t".    # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		    
	    } # End of has_in_cterm
	    
	    if ($has_rh) {
		
		#/////////
		# NOT SONG
		#\\\\\\\\\
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".                # Source
		    "rnaseh_domain\t".           # Data type
		    "$lf_rh_dom_start\t".        # Start
		    "$lf_rh_dom_end\t".          # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".                # Source
		    "rnaseh_orf\t".              # Data type
		    "$lf_rh_orf_start\t".        # Start
		    "$lf_rh_orf_end\t".          # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
	    }
	    
	    if ($has_rt) {
		
		#/////////
		# NOT SONG
		#\\\\\\\\\
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".                # Source
		    "rt_domain\t".               # Data type
		    "$lf_rt_dom_start\t".        # Start
		    "$lf_rt_dom_end\t".          # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
		$gff_str_out = "$lf_seq_id\t".     # Seq ID
		    "$gff_src\t".                # Source
		    "rt_orf\t".                  # Data type
		    "$lf_rt_orf_start\t".        # Start
		    "$lf_rt_orf_end\t".          # End
		    "$lf_score\t".               # Score
		    "$lf_strand\t".              # Strand
		    ".\t".                       # Frame
		    "ltr_finder_$lf_ltr_id\n";   # Retro ID
		print GFFOUT $gff_str_out;
		
	    } # End of has_reverse_transcriptase
	    
	    #-----------------------------+
	    # LOAD ID VAR                 |
	    #-----------------------------+
	    $lf_ltr_id = $1;

	}
	
	# SEQ ID AND LENGTH
	elsif (m/>Sequence: (.*) Len:(.*)/){
	    $lf_seq_id = $1;
	    $lf_seq_len = $2;
	}
	
	# SPAN LOCATION, LENGTH, AND STRAND
	elsif (m/^Location : (\d*) - (\d*) Len: (\d*) Strand:(.)/){
	    $lf_span_start = $1;
	    $lf_span_end = $2;
	    $lf_length = $3;
	    $lf_strand = $4;
	}
	
	# SCORE SIMILARITY
	elsif (m/^Score    : (.*) \[LTR region similarity:(.*)\]/){
	    $lf_score = $1;
	    $lf_ltr_similarity = $2;
	}
	
	# STATUS SET
	elsif (m/^Status   : (\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)(\d)/){
	    # Since this is a binary string, it can be split as digits
	    # and used to load the $has_* booleans
	    $has_5ltr_tg = $1;
	    $has_5ltr_ca = $2;
	    $has_3ltr_tg = $3;
	    $has_3ltr_ca = $4;
	    $has_tsr = $5;
	    $has_pbs = $6;
	    $has_ppt = $7;
	    $has_rt = $8;
	    $has_in_core = $9;
	    $has_in_cterm = $10;
	    $has_rh = $11;
	}
	
	# 5' LTR
	elsif (m/^5\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	    $lf_5ltr_start = $1;
	    $lf_5ltr_end = $2;
	    $lf_5ltr_len = $3;
	}
	
	# 3' LTR
	elsif (m/^3\'-LTR   : (\d*) - (\d*) Len: (\d*)/){
	    $lf_3ltr_start = $1;
	    $lf_3ltr_end = $2;
	    $lf_3ltr_len = $3;
	}
	
    # TARGET SITE REPLICATION
	elsif (m/TSR      : (\d*) - (\d*) , (\d*) - (\d*) \[(.*)\]/){
	    $lf_5tsr_start = $1;
	    $lf_5tsr_end = $2;
	    $lf_3tsr_start = $3;
	    $lf_3tsr_end = $4;
	    $lf_tsr_string = $5;
	}
	
	# SHARPNESS METRIC
	elsif (m/^Sharpness: (.*),(.*)/){
	    $lf_sharp_5 = $1;
	    $lf_sharp_3 = $2;
	}
	
	# PBS
	elsif (m/PBS   : \[(\d*)\/(\d*)\] (\d*) - (\d*) \((.*)\)/) {
	    $lf_pbs_num_match = $1;
	    $lf_pbs_aln_len = $2;
	    $lf_pbs_start = $3;
	    $lf_pbs_end = $4;
	    $lf_pbs_trna = $5;
	}
	
	# PPT
	elsif (m/PPT   : \[(\d*)\/(\d*)\] (\d*) - (\d*)/) {
	    $lf_ppt_num_match = $1;
	    $lf_ppt_aln_len = $2;
	    $lf_ppt_start = $3;
	    $lf_ppt_end = $4;
	}
	
	# PROTEIN DOMAINS
	# This will need to be modified and checked after another run
	# using ps_scan to get the additional domains
	#
	#Domain: 56796 - 57326 [possible ORF:56259-59144, (IN (core))]
	elsif (m/Domain: (\d*) - (\d*) \[possible ORF:(\d*)-(\d*), \((.*)\)\]/) {
	    
	    $lf_domain_dom_start = $1;
	    $lf_domain_dom_end = $2;
	    $lf_domain_orf_start = $3;
	    $lf_domain_orf_end = $4;
	    $lf_domain_name = $5;
	    
	    # Temp while I work with this data
	    #print "DOMAIN:".$lf_domain_name."\n";
	    
	    if ($lf_domain_name =~ 'IN \(core\)') {
		
		$lf_in_core_dom_start = $lf_domain_dom_start;
		$lf_in_core_dom_end = $lf_domain_dom_end;
		$lf_in_core_orf_start = $lf_domain_orf_start;
		$lf_in_core_orf_end = $lf_domain_orf_end;
		
		# Temp while debug
		# This is to check the regexp vars can be fetched here
		#print "\tDom Start: $lf_in_core_dom_start\n";
		#print "\tDom End:   $lf_in_core_dom_end\n";
		#print "\tOrf Start: $lf_in_core_orf_start\n";
		#print "\tOrf End:   $lf_in_core_orf_end\n";
		
	    }
	    elsif ($lf_domain_name =~ 'IN \(c-term\)') {
		$lf_in_cterm_dom_start = $lf_domain_dom_start;
		$lf_in_cterm_dom_end = $lf_domain_dom_end;
		$lf_in_cterm_orf_start = $lf_domain_orf_start;
		$lf_in_cterm_orf_end = $lf_domain_orf_end;
	    }
	    elsif ($lf_domain_name =~ 'RH') {
		$lf_rh_dom_start = $lf_domain_dom_start;
		$lf_rh_dom_end = $lf_domain_dom_end;
		$lf_rh_orf_start = $lf_domain_orf_start;
		$lf_rh_orf_end = $lf_domain_orf_end;
	    }
	    elsif ($lf_domain_name =~ 'RT') {
		
		$lf_rt_dom_start = $lf_domain_dom_start;
		$lf_rt_dom_end = $lf_domain_dom_end;
		$lf_rt_orf_start = $lf_domain_orf_start;
		$lf_rt_orf_end = $lf_domain_orf_end;  
		
	    }
	    else {
		print "\a";
		print STDERR "Unknown domain type: $lf_domain_name\n";
	    }
	    
	    
	    
	} # End of elsif Domain
	
	#-----------------------------+
	# FILE HEADER INFORMATION     |
	#-----------------------------+
	
	# PROGRAM NAME
	elsif (m/^Program    : (.*)/) {
	    $lf_prog_name = $1;
	}
	
	# PROGRAM VERSION
	elsif (m/^Version    : (.*)/) {
	    $lf_version = $1;
	}



	




	# MOD HERE
	if ($print_gff_out) {

	}












	
    }
    
    close INFILE;
    
    #-----------------------------+
    # PRINT LAST GFFOUT           |
    #-----------------------------+
    
    # FULL SPAN
    $gff_str_out = "$lf_seq_id\t". # Seq ID
	"$gff_src\t".                # Source
	"LTR_retrotransposon\t".     # Data type
	"$lf_span_start\t".          # Start
	"$lf_span_end\t".            # End
	"$lf_score\t".               # Score
	"$lf_strand\t".              # Strand
	".\t".                       # Frame
	"ltr_finder_$lf_ltr_id\n";   # Retro ID
    print GFFOUT $gff_str_out;
    
    
    if ($has_pbs) {
	$gff_str_out = "$lf_seq_id\t".  # Seq ID
	    "$gff_src\t".               # Source
	    "primer_binding_site\t".    # Data type
	    "$lf_pbs_start\t" .         # Start
	    "$lf_pbs_end\t".            # End
	    "$lf_score\t".              # Score
	    "$lf_strand\t".             # Strand
	    ".\t".                      # Frame
	    "ltr_finder_$lf_ltr_id\n";  # Retro ID
	print GFFOUT $gff_str_out;
    }
    
    if ($has_ppt) {
	$gff_str_out = "$lf_seq_id\t".  # Seq ID
	    "$gff_src\t".               # Source
	    "RR_tract\t".               # Data type
	    "$lf_ppt_start\t".          # Start
	    "$lf_ppt_end\t".            # End
	    "$lf_score\t".              # Score
	    "$lf_strand\t".             # Strand
	    ".\t".                      # Frame
	    "ltr_finder_$lf_ltr_id\n";  # Retro ID
	print GFFOUT $gff_str_out;
    }
    
    
    $gff_str_out = "$lf_seq_id\t".  # Seq ID
	"$gff_src\t".               # Source
	"five_prime_LTR\t".         # Data type
	"$lf_5ltr_start\t".         # Start
	"$lf_5ltr_end\t".           # End
	"$lf_score\t".              # Score
	"$lf_strand\t".             # Strand
	".\t".                      # Frame
	"ltr_finder_$lf_ltr_id\n";  # Retro ID
     print GFFOUT $gff_str_out;
    
    $gff_str_out = "$lf_seq_id\t".  # Seq ID
	"$gff_src\t".               # Source
	"three_prime_LTR\t".        # Data type
	"$lf_3ltr_start\t".         # Start
	"$lf_3ltr_end\t".           # End
	"$lf_score\t".              # Score
	"$lf_strand\t".             # Strand
	".\t".                      # Frame
	"ltr_finder_$lf_ltr_id\n";  # Retro ID
    print GFFOUT $gff_str_out;
    
    if ($has_tsr) {
	
	$gff_str_out = "$lf_seq_id\t".  # Seq ID
	    "$gff_src\t".                # Source
	    "target_site_duplication\t". # Data type
	    "$lf_5tsr_start\t".          # Start
	    "$lf_5tsr_end\t".            # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;	    
	
	$gff_str_out = "$lf_seq_id\t".   # Seq ID
	    "$gff_src\t".                # Source
	    "target_site_duplication\t". # Data type
	    "$lf_3tsr_start\t".          # Start
	    "$lf_3tsr_end\t".            # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
    }
    
    
    # Integrase Core
    if ($has_in_core) {
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".              # Source
	    "integrase_core_domain\t".  # Data type
	    "$lf_in_core_dom_start\t".   # Start
	    "$lf_in_core_dom_end\t".     # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".              # Source
	    "integrase_core_orf\t".      # Data type
	    "$lf_in_core_orf_start\t".   # Start
	    "$lf_in_core_orf_end\t".     # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
    } # End of has in_core
    
    
    if ($has_in_cterm) {
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "integrase_cterm_domain\t".  # Data type
	    "$lf_in_cterm_dom_start\t".  # Start
	    "$lf_in_cterm_dom_end\t".    # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "integrase_cterm_orf\t".     # Data type
	    "$lf_in_cterm_orf_start\t".  # Start
	    "$lf_in_cterm_orf_end\t".    # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
    } # End of has_in_cterm
    
    if ($has_rh) {
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rnaseh_domain\t".           # Data type
	    "$lf_rh_dom_start\t".        # Start
	    "$lf_rh_dom_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rnaseh_orf\t".              # Data type
	    "$lf_rh_orf_start\t".        # Start
	    "$lf_rh_orf_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
    }
    
    if ($has_rt) {
	
	#/////////
	# NOT SONG
	#\\\\\\\\\
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rt_domain\t".               # Data type
	    "$lf_rt_dom_start\t".        # Start
	    "$lf_rt_dom_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
	$gff_str_out = "$lf_seq_id\t".     # Seq ID
	    "$gff_src\t".                # Source
	    "rt_orf\t".                  # Data type
	    "$lf_rt_orf_start\t".        # Start
	    "$lf_rt_orf_end\t".          # End
	    "$lf_score\t".               # Score
	    "$lf_strand\t".              # Strand
	    ".\t".                       # Frame
	    "ltr_finder_$lf_ltr_id\n";   # Retro ID
	print GFFOUT $gff_str_out;
	
    } # End of has_reverse_transcriptase
    
    
    #-----------------------------+
    # PRINT ADDITIONAL DATA       |
    # IF REQUESTED                |
    #-----------------------------+
    if ($verbose) {
	print STDERR "\n\n+-----------------------------+\n";
	print STDERR " RESULTS SUMMARY\n";
	print STDERR " $lf_prog_name\n";
	print STDERR " $lf_version\n";
	print STDERR "+-----------------------------+\n";
	
	# ONLY SHOWS LAST PARSED\
	print STDERR "\n";
	print STDERR "SEQ_ID:\t\t$lf_seq_id\n";
	print STDERR "SEQ_LEN:\t$lf_seq_len\n";
	print STDERR "LTR_ID:\t\t$lf_ltr_id\n";
	print STDERR "LTRSTART:\t$lf_span_start\n";
	print STDERR "LTREND:\t\t$lf_span_end\n";
	print STDERR "LTRLEN:\t\t$lf_length\n";
	print STDERR "STRAND:\t\t$lf_strand\n";
	print STDERR "SCORE:\t\t$lf_score\n";
	print STDERR "SIMIL:\t\t$lf_ltr_similarity\n";
	
	print STDERR "\nLTR LOCATION:\n";
	print STDERR "5LTR START:\t$lf_5ltr_start\n";
	print STDERR "5LTR END:\t$lf_5ltr_end\n";
	print STDERR "5LTR LEN:\t$lf_5ltr_len\n";
	print STDERR "3LTR START:\t$lf_3ltr_start\n";
	print STDERR "3LTR END:\t$lf_3ltr_end\n";
	print STDERR "3LTR LEN:\t$lf_3ltr_len\n";
	
	print STDERR "\nTSR COORDINATES:\n";
	print STDERR "5-TSR START:\t$lf_5tsr_start\n";
	print STDERR "5-TSR END:\t$lf_5tsr_end\n";
	print STDERR "3-TSR START:\t$lf_3tsr_start\n";
	print STDERR "3-TSR END:\t$lf_3tsr_end\n";
	print STDERR "TSR-STRING:\t$lf_tsr_string\n";
    
	print STDERR "\nMETRICS:\n";
	print STDERR "5 SHARP: $lf_sharp_5\n";
	print STDERR "3 SHARP: $lf_sharp_3\n";
	
	print STDERR "\nPBS:\n";
	print STDERR "NUM MATCH:\t$lf_pbs_num_match\n" if $lf_pbs_num_match;
	print STDERR "ALN LENGT:\t$lf_pbs_aln_len\n" if $lf_pbs_aln_len;
	print STDERR "PBS START:\t$lf_pbs_start\n" if $lf_pbs_start;
	print STDERR "PBS END:\t$lf_pbs_end\n" if $lf_pbs_end;
	print STDERR "PBS tRNA:\t$lf_pbs_trna\n" if $lf_pbs_trna;
	
	print STDERR "\nPPT:\n";
	print STDERR "NUM MATCH:\t$lf_ppt_num_match\n";
	print STDERR "ALN LENT:\t$lf_ppt_aln_len\n";
	print STDERR "PPT START:\t$lf_ppt_start\n";
	print STDERR "PPT END:\t$lf_ppt_end\n";
	
	print STDERR "\nCOMPONENT STATUS:\n";
	print STDERR "5LTR_TG\n" if $has_5ltr_tg;
	print STDERR "5LTR_CA\n" if $has_5ltr_ca;
	print STDERR "3LTR_TG\n" if $has_3ltr_tg;
	print STDERR "3LTR_CA\n" if $has_3ltr_ca;
	print STDERR "TSR\n" if $has_tsr;
	print STDERR "PBS\n" if $has_pbs;
	print STDERR "PPT\n" if $has_ppt;
	print STDERR "RT\n" if $has_rt;
	print STDERR "IN(core)\n" if $has_in_core;
	print STDERR "IN(c-term)\n" if $has_in_cterm;
	print STDERR "RH\n" if $has_rh;
    }
    
    close GFFOUT;

} # End of ltrfinder2gff subfunction




1;
__END__

=head1 NAME

cnv_ltrfinder2gff.pl - Converts LTR_Finder output to gff format

=head1 VERSION

This documentation refers to program version $Rev: 761 $

=head1 SYNOPSIS

=head2 Usage

    cnv_ltrfinder2gff.pl -i lf_result.txt -o lf_result.gff

=head2 Required Arguments

    --infile        # Path to the input file
                    # Result from a single record fasta file
    --outdir        # Base output dir

=head1 DESCRIPTION

Convert the ltrfinder output to gff format. This assumes that the output
from ltr_finder correspons to a single BAC. A suffix can be passed with the
--suffix option to provide a suffix for the source column. For example
run default ltr_finder results with --suffix def to create 
ltr_finder:def. T

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--infile

Path of the input file. If an input file is not provided, the program
will expect input from STDIN.

=item -o,--outfile

Path of the output file.  If an output path is not provided,
the program will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item -n,--name

The sequence name to use in the GFF output file. Otherwise, this will
just use 'seq' as the sequence name.

=item -p, --param

The name of the paramter set used. This will be appened to the data in the
second column of the GFF file, and can be used to distinguish among 
parameter combinations
for multiple applications of ltrfinder to the same sequence file.

=item --apend

Append the GFF output to the gff file indicate by the --outfile option. This
allows you to append results from multiple programs to a single GFF file, but
it is generally recommended to create separate GFF files and concatenate
them at a later stage.

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

Error messages that you may encounter and possible solutions are listed below:

=over 2

=item Expecting input from STDIN

If a file is not specified by the -i or --infile option, the program will
expect to receive intput from standard input.

=back

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuartion file or variables set 
in the user's environment.

=head1 DEPENDENCIES

=head2 Software

This program requires the following software:

=over

=item * LTR Finder

This program parses output from the LTR finder program. It is possible to
obtain a linux binary by contacting the authors : xuzh <at> fudan.edu.cn.
It is also possible to obtain these results using the LTR_FINDER web page:
http://tlife.fudan.edu.cn/ltr_finder/

=back

=head2 Perl Modules

This program does not make use of perl modules beyond those installed
with the basic Perl package. If you discover a dependency that is not
documented here, please email the author or file a bug report.

=head1 BUGS AND LIMITATIONS

Any known bugs and limitations will be listed here.

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limitations

=over 2

=item * Limited Testing with LTR_FINDER versions

This program is known to parse the results from LTR_FINDER v 1.0.2. This 
program has not been tested with the results from the LTR_FINDER web site.

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

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/14/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/14/2007
# - Program started
# 10/01/2007
# - Moved main body of program to subfunction ltrfinder2gff
#   This subfunction operates under the assumption that all
#   results in the ltrfinder output file are the results for
#   a single assembly. This will facilitate copying the results
#   to an individual directory for the contig being annotated.
# - Making gff output, making this song complient
# 10/02/2007
# - Finishing gff output, now saving to string to write to
#   both gffout and stdout.
# 01/27/2009
# - Modifying program to use the print_help subfunction that
#   extracts the help message from the POD documentation
# - Updated POD documentation
# - Modified to accept intput from STDIN when --infile not
#   specified at the command line
# - Modified to write output to STOUT when --outfile not
#   specified at the command line
# 01/28/2009
# - Finished update of POD documentation
#
# 03/27/2009
# - Renamed --name to --seqname
# - Fixed program to accept seqname to override parsed name
# - Added program
#
# 04/27/2009
# - modified the parser to work with a multiple record fasta
#   file
# - This was done by removing the statement where I did not
#   print output when the ltr_finder id was set to one
# - I added a logical to check for values, and this will output
#   GFF only when values are present
