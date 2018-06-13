#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# batch_tenest.pl - Run TeNest in batch mode                |
#                                                           |
#-----------------------------------------------------------+
#                                                           |
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_@_gmail.com                          |
# STARTED: 04/29/2008                                       |
# UPDATED: 03/24/2009                                       |
#                                                           |
# DESCRIPTION:                                              |
#  Run TeNest in batch mode.                                |
#                                                           |
# USAGE:                                                    |
#  ShortFasta Infile.fasta Outfile.fasta                    |
#                                                           |
# VERSION: $Rev: 777 $                                      |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#-----------------------------------------------------------+
# TODO: -Convert results to gff format
#       -additional command line options for te_nest
#
# NOTE: TENEST Assumes that directories do not end with /
#

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

#-----------------------------+
# PROGRAM VARIABLES           |
#-----------------------------+
my ($VERSION) = "Release 1.0";

#-----------------------------+
# VARIABLE SCOPE              |
#-----------------------------+
my $indir;
my $outdir;

# Booleans
my $quiet = 0;
my $verbose = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;

# TE NEST VARIABLES WITH DEFAULT VALUES
my $organism_db = "wheat";
my $home_dir = $ENV{'HOME'};
#my $te_db_dir = "$home_dir/apps/te_nest";

# LOCATION OF TE NEST
# Look for env variable or assume in PATH
my $wublast_dir = $ENV{'DP_WUBLAST_DIR'} || "/usr/local/genome/wu_blast/";
my $te_nest_bin = $ENV{'TE_NEST_BIN'} || "TEnest.pl";
my $te_db_dir = $ENV{'TE_NEST_DIR'} || $ENV{'HOME'};

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|indir=s"    => \$indir,
                    "o|outdir=s"   => \$outdir,
		    # ADDITIONAL OPTIONS
		    "org=s"        => \$organism_db,
		    "blast-dir=s"  => \$wublast_dir,
		    "tenest-bin=s" => \$te_nest_bin,
		    "tenest-dir=s" => \$te_db_dir,
		    "q|quiet"      => \$quiet,
		    "verbose"      => \$verbose,
		    # ADDITIONAL INFORMATION
		    "usage"        => \$show_usage,
		    "version"      => \$show_version,
		    "man"          => \$show_man,
		    "h|help"       => \$show_help,);

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
    print "\nbatch_tenest.pl:\n".
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
    print_help("full");
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

#unless ($te_db_dir =~ /\/$/ ) {
#    $te_db_dir = $te_db_dir."/";
#}

#-----------------------------+
# CHECK FOR REQUIRED PROGRAMS |
#-----------------------------+
unless (-e $te_db_dir."/TEnest.pl") {
    die "TENest.pl does not exist at:\n$te_db_dir\n"
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

#-----------------------------+
# SHOW ERROR IF NO FILES      |
# WERE FOUND IN THE INPUT DIR |
#-----------------------------+
if ($num_files == 0) {
    print STDERR "\a";
    print STDERR "\nERROR: No fasta files were found in the input directory\n".
        "$indir\n".
        "Fasta files must have the fasta or fa extension.\n\n";
    exit;
}

#-----------------------------+
# CREATE THE OUT DIR          |
# IF IT DOES NOT EXIST        |
#-----------------------------+
unless (-e $outdir) {
    print STDERR "Creating output dir ...\n" unless $quiet;
    mkdir $outdir ||
        die "Could not create the output directory:\n$outdir";
}


#-----------------------------------------------------------+ 
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+ 
my $proc_num = 0;
my $file_num =0;

for my $ind_file (@fasta_files)
{

    my $name_root;

    $proc_num++;
    $file_num++;

    #-----------------------------+
    # GET THE ROOT NAME OF THE    |
    # FASTA FILE                  |
    #-----------------------------+
    if ($ind_file =~ m/(.*)\.hard\.fasta$/) {
        # file ends in .hard.fasta
        # This is hard masked fasta files
        $name_root = "$1";
    }
    elsif ($ind_file =~ m/(.*)\.masked\.fasta$/) {
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
    # CREATE ROOT NAME DIR        |
    #-----------------------------+
    my $name_root_dir = $outdir.$name_root;
    unless (-e $name_root_dir) {
        mkdir $name_root_dir ||
            die "Could not create dir:\n$name_root_dir\n"
	}

    #-----------------------------+
    # CREATE TENEST OUTDIR        |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $tenest_out_dir = $outdir.$name_root."/tenest";
    unless (-e $tenest_out_dir) {
        mkdir $tenest_out_dir ||
            die "Could not create tenest out dir:\n$tenest_out_dir\n";
    }

    #-----------------------------+
    # CREATE GFF OUTDIR           |
    #-----------------------------+
    # Dir to hold gene prediction output from local software
    my $gff_out_dir = $outdir.$name_root."/gff/";
    #print STDERR "$gff_out_dir\n";
    unless (-e $gff_out_dir) {
        mkdir $gff_out_dir ||
            die "Could not create gff out dir:\n$gff_out_dir\n";
    }

    #-----------------------------+
    # RUN TE NEST                 |
    #-----------------------------+

    # The following does work
    # 04/28/2009
    # TE NEST COMMAND
    my $te_nest_cmd = $te_nest_bin." ".$indir.$ind_file.
	" --org ".$organism_db.
	" --output ".$tenest_out_dir;

# TO DO: Get the following working    
# Trying usage as program OPTIONS INFILE
#    my $te_nest_cmd = $te_db_dir."/TEnest.pl ".
#	" --output ".$tenest_out_dir.
#	$indir.$ind_file
#	" --org ".$organism_db;
#	" --blast ".$wublast_dir.
#	" --current ".$te_db_dir;


    system ($te_nest_cmd);

    #-----------------------------+
    # CONVERT THE OUTPUT TO       |
    # GFF FORMAT                  |
    #-----------------------------+
    my $te_nest_result = $tenest_out_dir."/TE_DIR/".$name_root.".LTR";

    # check to see if this file name contains the masked
    # kludge
    unless (-e $te_nest_result) {
	$te_nest_result = $tenest_out_dir."/TE_DIR/".$name_root.
	    ".masked.LTR";
    }
    

    my $gff_out_path =  $gff_out_dir.$name_root.".tenest.gff";
    my $do_gff_append = 0;
    my $gff_seqname = $name_root;
    
    tenest2gff ($te_nest_result, $gff_out_path, 
		$do_gff_append, $gff_seqname);
    
}

exit;

#-----------------------------------------------------------+ 
# SUBFUNCTIONS                                              |
#-----------------------------------------------------------+


sub tenest2gff {
# CONVERT BLAST TO GFF 
    
    # seqname - ID to use in the seqname field
    # tenestin - path to the blast input file
    # gffout  - path to the gff output file
    # append  - boolean append data to existing file at gff out
    my ($tenestin, $gffout, $append, $seqname) = @_;
    my $tename;           # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit

    my $num_te_data;      # Length of the te data array

    my $i;                # Array index val
    my $j;                # Array index val

    # Initialize counters
    my $numfrag = 0;
    my $numpair = 0;
    my $numsolo = 0;
    my $numnltr = 0;
    my $pair_line = 0;

    # BOOLEANS
    my $in_frag = 0;
    my $in_pair = 0;
    my $in_solo = 0;
    my $in_nltr = 0;

    # SOLO VARS
    # These need to be available outside of the
    # initial observation
    my $sol_entry;
    my $sol_number;
    my $sol_type;
    my $sol_dir;
    my $sol_nest_group;
    my $sol_nest_order;
    my $sol_nest_level;
    
    # PAIR VARS
    my $pair_number;
    my $pair_type;
    my $pair_dir;
    my $pair_bsr;
    my $pair_nest_group;
    my $pair_nest_order;
    my $pair_nest_level;

    # FRAG VARS
    my $frag_number;
    my $frag_type;
    my $frag_dir;
    my $frag_nest_group;
    my $frag_nest_order;
    my $frag_nest_level;

    # NLTR VARS
    my $nltr_number;
    my $nltr_type;
    my $nltr_dir;
    my $nltr_nest_group;
    my $nltr_nest_order;
    my $nltr_nest_level;

    # Open the BLAST report object
    open (TENESTIN,"<$tenestin")
	|| die "Could not open TE-NEST input file:\n$tenestin\n";
    
    # Open file handle to the gff outfile    
    if ($append) {
	open (GFFOUT, ">>$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    else {
	open (GFFOUT, ">$gffout") 
	    || die "Can not open file:\n $gffout\n";
    }
    
    
    while (<TENESTIN>) {
	chomp;                   # Remove line endings
	my @te_data = split;   # Split by spaces
	
	# Four annotation types
	#SOLO (solo LTRs), 
	#PAIR (full LTR retrotransposons), 
	#FRAG (fragmented TEs of all types), 
	#NLTR (Full length non-LTR containing TEs)
 
	# Temp show the line being processed
	#print STDERR "$_\n";

	#-----------------------------------------------------------+
	# ADDITIONAL FEATURE DATA                                   |
	#-----------------------------------------------------------+

	#-----------------------------+
	# ADDITIONAL SOLO DATA        |
	#-----------------------------+
	if ($in_solo) {
      
	    $in_solo = 0;
	    $num_te_data = @te_data;

	    if ($verbose) {
		print STDERR $te_data[0]."\n";
		print STDERR "\tNum in array:\t$num_te_data\n";
	    }

	    # i is used to index in the parent array
	    # j is used to index the sol_coords array where
	    # 1,2,3,4 is seq_start, seq_end, te_start, te_end
	    #
	    my @sol_coords = ();
	    $j = 0;
	    for ($i=1; $i<$num_te_data; $i++) {
		$j++;
		$sol_coords[$j] = $te_data[$i];

		if ($verbose) {
		    print STDERR "\t$i:\t".$te_data[$i]."\n";
		}

		if ($j == 4) {

		    # Appending sol_num to sol type will give a
		    # unique name for every occurrence in the gff file
		    my $sol_name = "solo_".$sol_type."_".$sol_number;
		    # Print output to gff file
		    print GFFOUT 
			"$seqname\t".                # Seqname
			"tenest\t".                  # Source
			"exon\t".                    # Feature type name
			$sol_coords[1]."\t".         # Start
			$sol_coords[2]."\t".         # End
			".\t".                       # Score
			# $sol_dir may give strand
			".\t".                 # Strand
			".\t".                       # Frame
			"$sol_name\n";                # Feature name
		    
		    $j = 0;
		    #$sol_coords=();

		} # End of if $j==4
	    } # End of for $i
	} # End of if $in_solo

	#-----------------------------+
	# ADDITIONAL PAIR DATA        |
	#-----------------------------+
	elsif ($in_pair) {

	    # IS BSR Substitution rate ?

	    # Pair has three additional lines of info
	    $pair_line++;

	    # Get additional info
	    
	    #-----------------------------+
	    # LEFT LTR                    |
	    #-----------------------------+
	    if ($te_data[1] =~ 'L') {
		
		$num_te_data = @te_data;

		if ($verbose) {
		    print STDERR $te_data[0]."\n";
		    print STDERR "\tNum in array:\t$num_te_data\n";
		}

		my @l_pair_coords = ();
		$j = 0;
		for ($i=2; $i<$num_te_data; $i++) {
		    $j++;
		    $l_pair_coords[$j] = $te_data[$i];
		    
		    if ($verbose) {
			print STDERR "\t$i:\t".$te_data[$i]."\n";
		    }
		    
		    if ($j == 4) {
			
			my $l_start;
			my $l_end;
			
			# Make start less then end
			if ($l_pair_coords[1] < $l_pair_coords[2]) {
			    $l_start = $l_pair_coords[1];
			    $l_end = $l_pair_coords[2];
			}
			else {
			    $l_start = $l_pair_coords[2];
			    $l_end = $l_pair_coords[1];
			}

			# Appending sol_num to sol type will give a
			# unique name for every occurrence in the gff file
			my $pair_name = "pair_".$pair_type."_".$pair_number;
			# Print output to gff file
			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "tenest\t".                  # Source
			    "exon\t".                    # Feature type name
			    "$l_start\t".                # Start
			    "$l_end\t".                  # End
			    ".\t".                       # Score
			    # $pair_dir may give strand
			    ".\t".                       # Strand
			    ".\t".                       # Frame
			    "$pair_name\n";              # Feature name
			
			$j = 0;
			#$sol_coords=();
			
		    } # End of if $j==4
		} # End of for $i
		
		
	    }

	    #-----------------------------+
	    # RIGHT LTR                   |
	    #-----------------------------+
	    elsif ($te_data[1] =~ 'R') {


		$num_te_data = @te_data;

		if ($verbose) {
		    print STDERR $te_data[0]."\n";
		    print STDERR "\tNum in array:\t$num_te_data\n";
		}		

		my @r_pair_coords = ();
		$j = 0;
		for ($i=2; $i<$num_te_data; $i++) {
		    $j++;
		    $r_pair_coords[$j] = $te_data[$i];
		    
		    if ($verbose) {
			print STDERR "\t$i:\t".$te_data[$i]."\n";
		    }		    

		    if ($j == 4) {
			
			my $r_start;
			my $r_end;
			
			# Make start less then end
			if ($r_pair_coords[1] < $r_pair_coords[2]) {
			    $r_start = $r_pair_coords[1];
			    $r_end = $r_pair_coords[2];
			}
			else {
			    $r_start = $r_pair_coords[2];
			    $r_end = $r_pair_coords[1];
			}

			# Appending sol_num to sol type will give a
			# unique name for every occurrence in the gff file
			my $pair_name = "pair_".$pair_type."_".$pair_number;
			# Print output to gff file
			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "tenest\t".                  # Source
			    "exon\t".                    # Feature type name
			    "$r_start\t".                # Start
			    "$r_end\t".                  # End
			    ".\t".                       # Score
			    # $pair_dir may give strand
			    ".\t".                       # Strand
			    ".\t".                       # Frame
			    "$pair_name\n";              # Feature name
			
			$j = 0;
			#$sol_coords=();
			
		    } # End of if $j==4
		} # End of for $i
	

	    }
	    
	    #-----------------------------+
	    # MIDDLE                      |
	    #-----------------------------+
	    elsif ($te_data[1] =~ 'M') {
			
		$num_te_data = @te_data;
		if ($verbose) {
		    print STDERR $te_data[0]."\n";
		    print STDERR "\tNum in array:\t$num_te_data\n";
		}		

		my @m_pair_coords = ();
		$j = 0;
		for ($i=2; $i<$num_te_data; $i++) {
		    $j++;
		    $m_pair_coords[$j] = $te_data[$i];
		   
		    if ($verbose) {
			print STDERR "\t$i:\t".$te_data[$i]."\n";
		    }

		    if ($j == 4) {
			
			my $m_start;
			my $m_end;
			
			# Make start less then end
			if ($m_pair_coords[1] < $m_pair_coords[2]) {
			    $m_start = $m_pair_coords[1];
			    $m_end = $m_pair_coords[2];
			}
			else {
			    $m_start = $m_pair_coords[2];
			    $m_end = $m_pair_coords[1];
			}

			# Appending sol_num to sol type will give a
			# unique name for every occurrence in the gff file
			my $pair_name = "pair_".$pair_type."_".$pair_number;
			# Print output to gff file
			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "tenest\t".                  # Source
			    "exon\t".                    # Feature type name
			    "$m_start\t".                # Start
			    "$m_end\t".                  # End
			    ".\t".                       # Score
			    # $pair_dir may give strand
			    ".\t".                       # Strand
			    ".\t".                       # Frame
			    "$pair_name\n";              # Feature name
			
			$j = 0;
			
		    } # End of if $j==4
		} # End of for $i


	    }


	    # END OF PAIR DATA
	    if ($pair_line == 3) {
		$in_pair = 0;
		$pair_line = 0;

		# Print output to gff
	    }

	}

	#-----------------------------+
	# ADDITIONAL FRAG DATA        |
	#-----------------------------+
	elsif ($in_frag) {
	    $in_frag = 0;

	    # Get additinal info
      	    $num_te_data = @te_data;

	    if ($verbose) {
		print STDERR $te_data[0]."\n";
		print STDERR "\tNum in array:\t$num_te_data\n";
	    }

	    # i is used to index in the parent array
	    # j is used to index the sol_coords array where
	    # 1,2,3,4 is seq_start, seq_end, te_start, te_end
	    #
	    my @frag_coords = ();
	    $j = 0;
	    for ($i=1; $i<$num_te_data; $i++) {
		$j++;
		$frag_coords[$j] = $te_data[$i];

		if ($verbose) {
		    print STDERR "\t$i:\t".$te_data[$i]."\n";
		}

		if ($j == 4) {

		    # Appending sol_num to sol type will give a
		    # unique name for every occurrence in the gff file
		    my $frag_name = "frag_".$frag_type."_".$frag_number;
		    # Print output to gff file
		    print GFFOUT 
			"$seqname\t".                # Seqname
			"tenest\t".                  # Source
			"exon\t".                    # Feature type name
			$frag_coords[1]."\t".         # Start
			$frag_coords[2]."\t".         # End
			".\t".                       # Score
			# $sol_dir may give strand
			".\t".                 # Strand
			".\t".                       # Frame
			"$frag_name\n";                # Feature name
		    
		    $j = 0;

		} # End of if $j==4
	    } # End of for $i

	} # End of $in_frag

	#-----------------------------+
	# ADDITIONAL NLTR DATA        |
	#-----------------------------+
	elsif ($in_nltr) {
	    $in_nltr = 0;

	    # Get additinal info
      	    $num_te_data = @te_data;

	    if ($verbose) {
		print STDERR $te_data[0]."\n";
		print STDERR "\tNum in array:\t$num_te_data\n";
	    }

	    # i is used to index in the parent array
	    # j is used to index the sol_coords array where
	    # 1,2,3,4 is seq_start, seq_end, te_start, te_end
	    #
	    my @nltr_coords = ();
	    $j = 0;
	    for ($i=1; $i<$num_te_data; $i++) {
		$j++;
		$nltr_coords[$j] = $te_data[$i];

		if ($verbose) {
		    print STDERR "\t$i:\t".$te_data[$i]."\n";
		}		

		if ($j == 4) {

		    # Appending sol_num to sol type will give a
		    # unique name for every occurrence in the gff file
		    my $nltr_name = "frag_".$nltr_type."_".$nltr_number;
		    # Print output to gff file
		    print GFFOUT 
			"$seqname\t".                # Seqname
			"tenest\t".                  # Source
			"exon\t".                    # Feature type name
			$nltr_coords[1]."\t".         # Start
			$nltr_coords[2]."\t".         # End
			".\t".                       # Score
			# $sol_dir may give strand
			".\t".                 # Strand
			".\t".                       # Frame
			"$nltr_name\n";                # Feature name
		    
		    $j = 0;

		} # End of if $j==4
	    } # End of for $i

	} # End of $in_nltr



	#-----------------------------------------------------------+
	# NEW RECORD STARTING                                       |
	#-----------------------------------------------------------+

	#-----------------------------+
	# SOLO                        |
	#-----------------------------+
	if ($te_data[0] =~ 'SOLO') {
	    $numsolo++;
	    $in_solo = 1;              # Flip boolean to true

	    #$sol_entry = $te_data[0]; 
	    $sol_number = $te_data[1];
	    $sol_type = $te_data[2];
	    $sol_dir = $te_data[3];
	    $sol_nest_group = $te_data[4];
	    $sol_nest_order = $te_data[5];
	    $sol_nest_level = $te_data[6];

	}

	#-----------------------------+
	# PAIR                        |
	#-----------------------------+
	elsif ($te_data[0] =~ 'PAIR') {
	    $numpair++;
	    $in_pair = 1;              # Flip boolean to true

	    $pair_number = $te_data[1];
	    $pair_type = $te_data[2];
	    $pair_dir = $te_data[3]; 
	    $pair_bsr = $te_data[4];
	    $pair_nest_group = $te_data[5];
	    $pair_nest_order = $te_data[6]; 
	    $pair_nest_level = $te_data[7];

	}

	#-----------------------------+
	# FRAG                        |
	#-----------------------------+
	elsif ($te_data[0] =~ 'FRAG') {
	    $numfrag++;
	    $in_frag = 1;              # Flip boolean to true

	    $frag_number = $te_data[1]; 
	    $frag_type = $te_data[2]; 
	    $frag_dir = $te_data[3];
	    $frag_nest_group = $te_data[4]; 
	    $frag_nest_order = $te_data[5]; 
	    $frag_nest_level = $te_data[6];

	}
	
	#-----------------------------+
	# NLTR                        |
	#-----------------------------+
	elsif ($te_data[0] =~ 'NLTR') {
	    $numnltr++;
	    $in_nltr = 1;              # Flip boolean to true

	    #Non-LTR entry, 
	    $nltr_number = $te_data[1]; 
	    $nltr_type = $te_data[2]; 
	    $nltr_dir = $te_data[3]; 
	    $nltr_nest_group = $te_data[4]; 
	    $nltr_nest_order = $te_data[5]; 
	    $nltr_nest_level = $te_data[6];

	}

    } # End of while TESTIN
    

    # Show summary of counts if vebose
    if ($verbose) {
	print STDERR "\n";
	print STDERR "NUM SOLO:\t$numsolo\n";
	print STDERR "NUM PAIR:\t$numpair\n";
	print STDERR "NUM FRAG:\t$numfrag\n";
	print STDERR "NUM NLTR:\t$numnltr\n";
    }

    close GFFOUT;
    
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

batch_tenest.pl - Run TENest in batch mode

=head1 VERSION

This documentation refers to program version $Rev: 777 $

=head1 SYNOPSIS

=head2  Usage

    batch_tenest.pl -i indir -o outdir

=head2 Required Arguments

    -i,--indir         # Dir containing fasta files to process 
    -o,--outdir        # Dir to place the output in

=head1 DESCRIPTION

Runs the program TE Nest in batch mode. This will run TE Nest for all fasta
files located in the input directory.

=head1 REQUIRED ARGUMENTS

=over 2

=item -i,--indir

Path to the directory containing the input files to annotate with TEnest.
These files must be in the fasta format.

=item -o,--outdir

Path to the base output directory to place the batch_tenest.pl output.

=back

=head1 OPTIONS

=over 2

=item --tenest-bin

This is the path of the TEnest.pl program. This option can also be specified
with the TE_NEST_BIN variable in the user environment. If this option is
not given in the environment or at the command line, the program
will assume that TEnest.pl is locate in the user PATH.

=item --org

The organism database to use. The default organism database is maize.
This can be an organism directory in the TEnest database dir defined by
TE_NEST_DIR, or this can be the full path to the organism database to use.

=item --blast-dir

The directory for wu blast. If this directory is not specified at the
command line, batch_tenest.pl will assume that this is at
/usr/local/genome/wu_blast/

=item --tenest-dir

The directory that the TE Nest program is located in. If this directory is not
specified, batch_tenet.pl will assume that the TE Nest program is locate
in your home directory at $HOME/apps/te_nest/

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

=over

=item TENest.pl does not exist at: /some/directory/

The TENest.pl program does not exist at the directory you specified or
in the directory that the program assumes it to be in. The best way to
make sure this works is to specify a directory using the --tenest-dir option.

=item Can't open directory: /your/input/directory/

The input directory that you specified by the -i file has either been
typed in incorrectly, the path you entered does not exist, or you
do not have read access to that directory.

=item ERROR: No fasta files were found in the input directory

Currently this program assumes that all fasta files end with the 
fasta or fa extension. It is possible that you have fasta files in the
input directory, but they do not have the *.fasta or *.fa file names.

=back

=head1 CONFIGURATION AND ENVIRONMENT

=head2 Configuration File

This program does not currently make use of a configuration file.

=head2 Environment Options

This program makes use of the following options in the user environment.

=over 2

=item TE_NEST_BIN

This is path of the location of the TEnest.pl program. This option can
also be specified using the tenest-bin option at the command line. If this
is not specified in the environment or 

=item TE_NEST_DIR

This is the base path of the directory containing the subdirectories 
that have the organism specific directories of TE_NEST databases for
TE annotation with the TE nest program. This option can also be specified
with the --tenest-dir at the command line. If this option is not given
at the command line or using the TE_NEST_DIR variable in the environment,
the program will expect these datbase files to be located in the 
home directory.

=item DP_WUBLAST_DIR

The location of the wublast directory on your current machine. This option
can also be specified at the command line with the --blast-dir option.
If this option is not specified at the command line, the program will assume
that the blast directory is located at '/use/local/genome/wu_blast'.

=back

=head1 DEPENDENCIES

=head2 Required Software

=over 2

=item TE Nest

This program requires the installation of TE Nest. This program is 
available from:
http://www.public.iastate.edu/~imagefpc/Subpages/te_nest.html

=back

=head2 Required Databases

This program requires an organism database is the TE Nest format. This
databases are downloadable from:
http://www.public.iastate.edu/~imagefpc/Subpages/te_nest.html

=head1 BUGS AND LIMITATIONS

=head2 Bugs

=over 2

=item * No bugs currently known 

If you find a bug with this software, file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head2 Limiations

=over

=item * Please report any limitations you encounter

If you encounter limitations while working with this program, please 
file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

=back

=head1 SEE ALSO

This program is part of the DAWG-PAWS package of genome
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

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 09/29/2008

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 09/29/2008
# - Program started
# - Basic idea is to run the tenest program from 
#
# 11/06/2008
# - Added the ability to specify organism
# - Added the tenest2gff subfunction from cnv_tenest2gff.pl
# 
# 01/19/2009
# -Updated POD help file
