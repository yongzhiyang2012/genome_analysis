#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# cnv_tenest2gff.pl - Convert TENest Output to GFF          |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: JamesEstill_at_gmail.com                         |
# STARTED: 08/27/2007                                       |
# UPDATED: 04/03/2009                                       |
#                                                           |  
# DESCRIPTION:                                              | 
#  Convert TE Nest output to gff file format.               |
#                                                           |
# LICENSE:                                                  |
#  GNU General Public License, Version 3                    |
#  http://www.gnu.org/licenses/gpl.html                     |  
#                                                           |
#                                                           |
#-----------------------------------------------------------+
# TODO: -follow ontology names for: 
#    - LTR_fragment
#    - solo_LTR
#    - non-LTR
#    - pair-LTR
#  Allow override of this name using --feature

package DAWGPAWS;

#-----------------------------+
# INCLUDES                    |
#-----------------------------+
use strict;
use Getopt::Long;
use Bio::SearchIO;             # Parse BLAST output
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
# LOCAL VARIABLES             |
#-----------------------------+
my $feature_type;
# Set variable scope
my $seqname;                   # Can specify seq name in the command line
my $infile;
my $outfile;
my $param_name;                # Name for the tenest parameter set.
my $program = "tenest";        # The program, used in source col

# Booleans
my $verbose = 0;
my $help = 0;
my $quiet = 0;
my $show_help = 0;
my $show_usage = 0;
my $show_man = 0;
my $show_version = 0;
my $do_append = 0;

#-----------------------------+
# COMMAND LINE OPTIONS        |
#-----------------------------+
my $ok = GetOptions(# REQUIRED OPTIONS
		    "i|infile=s"  => \$infile,
                    "o|outfile=s" => \$outfile,
		    "n|s|seqname|name=s"  => \$seqname,
		    # OPTIONS
		    "feature=s"   => \$feature_type,
		    "program=s"   => \$program,
		    "p|param=s"   => \$param_name,
		    "verbose"     => \$verbose,
		    "append"      => \$do_append,
		    # BOOLEANS
		    "usage"       => \$show_usage,
		    "version"     => \$show_version,
		    "man"         => \$show_man,
		    "h|help"      => \$show_help,
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

if ($show_man) {
    # User perldoc to generate the man documentation.
    system ("perldoc $0");
    exit($ok ? 0 : 2);
}

if ($show_version) {
    print "\ncnv_tenest2gff.pl:\n".
	"Version: $VERSION\n\n";
    exit;
}

#-----------------------------------------------------------+
# MAIN PROGRAM BODY                                         |
#-----------------------------------------------------------+
tenest2gff ($do_append, $seqname, $infile, $outfile, $program, $param_name,
	    $feature_type);

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
    # feature - the type of feature it is
    my ($append, $seqname, $tenestin, $gffout, $source, $param_set, 
	$feature_type) 
	= @_;
    my $feature;          # The sequence feature
    my $tename;           # Name of the hit
    my $start;            # Start of the feature
    my $end;              # End of the feature
    my $strand;           # Strand of the hit

    my $num_te_data;      # Length of the te data array

    my $i;                # Array index val
    my $j;                # Array index val

    unless ($seqname) {
	$seqname = "seq";
    }

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

    # OPEN THE TENEST INTPUT FILE
    # Default is to expect input from standard input
    if ($tenestin) {
	open (TENESTIN,"<$tenestin")
	    || die "Could not open TE-NEST input file:\n$tenestin.\n";
    }
    else {
	print STDERR "Expecting input from STDIN\n";
	open (TENESTIN, "<&STDIN") ||
	    die "Can not accepte input from standard input.\n";
    }

    # OPEN OUTPUT FILE HANDLE
    # Default to STDOUT if no argument given
    if ($gffout) {
	if ($append) {
	    open (GFFOUT, ">>$gffout") 
		|| die "Can not open file:\n $gffout\n";
	}
	else {
	    open (GFFOUT, ">$gffout") 
		|| die "Can not open file:\n $gffout\n";
	}
    }
    else {
	open (GFFOUT, ">&STDOUT") ||
	    die "Can not print to STDOUT\n";
    }
    

    # SET THE SOURCE
    # This allows for the specification of the paramater set in the
    # gff result fijle
    #my $source = "tenest";
    if ($param_set) {
	$source = $source.":".$param_set;
    }

    while (<TENESTIN>) {
	chomp;                   # Remove line endings
	my @te_data = split;   # Split by spaces
	
	# set feature to default type
	$feature = "transposable_element";

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

		    #-----------------------------+
		    # SOLO LTR                    |
		    #-----------------------------+
		    if ($feature_type) {
			$feature = $feature_type;
		    }
		    else {
			$feature = "LTR_retrotransposon";
		    }

		    print GFFOUT 
			"$seqname\t".                # Seqname
			"$source\t".                 # Source
			"$feature\t".                # Feature type name
			$sol_coords[1]."\t".         # Start
			$sol_coords[2]."\t".         # End
			".\t".                       # Score
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

			#-----------------------------+
			# LEFT LTR                    |
			#-----------------------------+
			if ($feature_type) {
			    $feature = $feature_type;
			}
			else {
			    $feature = "five_prime_LTR_component";
			}

			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "$source\t".                 # Source
			    "$feature\t".                # Feature type name
			    "$l_start\t".                # Start
			    "$l_end\t".                  # End
			    ".\t".                       # Score
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

			#-----------------------------+
			# RIGHT LTR                   |
			#-----------------------------+
			if ($feature_type) {
			    $feature = $feature_type;
			}
			else {
			    $feature = "three_prime_LTR_component";
			}

			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "$source\t".                 # Source
			    "$feature\t".                # Feature type name
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

			#-----------------------------+
			# MIDDLE OF LTR RETRO         |
			#-----------------------------+
			if ($feature_type) {
			    $feature = $feature_type;
			}
			else {
			    $feature = "LTR_retrotransposon";
			}

			print GFFOUT 
			    "$seqname\t".                # Seqname
			    "$source\t".                 # Source
			    "$feature\t".                # Feature type name
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

		    #-----------------------------+
		    # LTR FRAGMENT                |
		    #-----------------------------+
		    if ($feature_type) {
			$feature = $feature_type;
		    }
		    else {
			$feature = "LTR_retrotransposon";
		    }

		    print GFFOUT 
			"$seqname\t".                # Seqname
			"$source\t".                 # Source
			"$feature\t".                 # Feature type name
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

		    # Appending nltr_num to sol type will give a
		    # unique name for every occurrence in the gff file
		    my $nltr_name = "nonltr_".$nltr_type."_".$nltr_number;
		    # Print output to gff file

		    #-----------------------------+
		    # TRANSPOSABLE ELEMENT        |
		    #-----------------------------+
		    if ($feature_type) {
			$feature = $feature_type;
		    }
		    else {
			$feature = "transposable_element";
		    }

		    print GFFOUT 
			"$seqname\t".                # Seqname
			"$source\t".                 # Source
			"$feature\t".                 # Feature type name
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

cnv_tenest2gff.pl - Convert TENest output to GFF

=head1 VERSION

This documentation refers to program version Release 1.0

=head1 SYNOPSIS

=head2 Usage

    cnv_tenest2gff.pl -i infile.txt -o outfile.gff

=head2 Required Arguments

    -i,--infile   # Path to the TE Nest result to convert
                  # Expects input from standard input otherwise
    -o,--outfile  # Path to the gff format outfile
                  # Writes output to standard output otherwise

=head1 DESCRIPTION

Converts TE Nest output to gff format output. The gff file will label 
these features as 'exon' to make them compatible with the apollo genome
annotation curation program.

=head1 REQUIRED ARGUMENTS

The following arguments are the most useful in using the cnv_tenest2gff.pl
program. They are not actually required since the cnv_tenest2gff.p program
will attempt to use a default name or path that makes sense.

=over 2

=item -i,--infile

Path of the input file to convert. If an input file is not provided, the program
will expect input from STDIN.

=item -o,--outfile

Path of the gff formatted output file that . If an output path is not provided,
the program will write output to STDOUT.

=back

=head1 OPTIONS

=over 2

=item -s,--seqname

The sequence name to use in the GFF output file. Otherwise, this will
just use 'seq' as the sequence name.

=item --program

The program used to generate the annotation result. This data is written
to the second colum in the gff output file. Be default this value is
set to be "tenest".

=item -p, --param

The name of the paramter set used. This will be appened to the data in the
second column, and can be used to distinguish among parameter combinations
for multiple applications of TE Nest to the same sequence file.

=item --feature

The name to use for the feature type. This is the output value in the third
colum of the gff output file. Be default this value is set to 
"transposable_element".

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

The error message that can be generated will be listed here.

=head1 CONFIGURATION AND ENVIRONMENT

This program does not make use of a configuration file or variables set
in the user's environment.

=head1 DEPENDENCIES

=head2 Required Software

=over 2

=item * TE Nest

This program requires output from the TE Nest web server, or from a local
copy of the TE nest program. The TE Nest web query is available at:
http://www.plantgdb.org/prj/TE_nest/TE_nest.html. You can download the
TE nest program at: 
http://www.public.iastate.edu/~imagefpc/Subpages/software.html.

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

=item * Please report limitations to using this software

If you find a limitation that makes it difficult to use this program, please
file a bug report on the DAWG-PAWS
Sourceforge website: http://sourceforge.net/tracker/?group_id=204962

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

GNU LESSER GENERAL PUBLIC LICENSE

GNU General Public License, Version 3

L<http://www.gnu.org/licenses/gpl.html>

THIS SOFTWARE COMES AS IS, WITHOUT ANY EXPRESS OR IMPLIED
WARRANTY. USE AT YOUR OWN RISK.

=head1 AUTHOR

James C. Estill E<lt>JamesEstill at gmail.comE<gt>

=head1 HISTORY

STARTED: 08/27/2007

UPDATED: 03/24/2009

VERSION: Release 1.0

=cut

#-----------------------------------------------------------+
# HISTORY                                                   |
#-----------------------------------------------------------+
#
# 08/27/2007
# -Program started
#
# 01/21/2009
# -Moved POD documentation to the end of the file
# -Updating POD documentation
# -Adding new help subfunction to extract help messages
#  from POD documentation
# -Added option to specify a parameter set name, this will
#  be added to the second column of the gff output file and
#  can be used to distinguish among parameter sets.
# -Modified to accept input from STDIN when --infile not 
#  specified
# -Modiied to write otput to STDOUT when --outfile not
#  specified
#
# 03/30/2009
# -Added support for -s,--seqname
#
# 04/03/2009
# -Added support for --feature
# -Fixed name col for non_ltrs, these were previously mislabeled
#  as fragments of LTRs
#
