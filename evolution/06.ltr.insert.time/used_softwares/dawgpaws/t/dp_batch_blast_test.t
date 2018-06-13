#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_blast_test.pl - Test batch_blast.pl              |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/23/2009                                       |
# UPDATED: 04/23/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_trf.pl program.                           |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 7;

##-----------------------------+
## TESTING TRF BIN             |
##-----------------------------+
## The batch_blast program tries to use the env var TRF_BIN
## otherwise it will assume the binary is inclue in the user PATH
diag("Testing blast ENV options ...");
my $blast_bin;
if ($ENV{DP_BLAST_BIN}) {
    ok ( $blast_bin=$ENV{DP_BLAST_BIN},
	 "blastall binary defined in environment as: $blast_bin") 
    }
else {
    ok ( $blast_bin="blastall", 
	 "blastall binary not defined in environment\n".
	 "expecting blast in path as: $blast_bin");
}

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the blastall program ...");
# Basic command to see that blastall is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_blast_cmd = $blast_bin." 2>&1";

# Get the result
my $basic_blast_result = `$basic_blast_cmd`;

# Test for the ltr version string that comes as part of the help message
like ($basic_blast_result, qr/blastall(.*)arguments/, "blastall is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Test running the batch_trf.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

my $batch_blast_cmd = "batch_blast.pl -i data/fasta/".
    " -c data/config/batch_blast_tab.jcfg".
    " -o $out_dir".
    " -d data/blast/";

# This will exit zero when things work
ok ( system($batch_blast_cmd)==0 , "batch_blast.pl");

#-----------------------------+
# TEST BLAST m8 PROGRAM OUTPUT|
#-----------------------------+
diag ("Testing batch_blast m8 blast results ..");
my $exp_file = "data/exp/HEX2986I03_TREP9_nr.bln";
open (EXPECTED, "<".$exp_file);
my @blast_m8_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE m8 BLAST OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX2986I03/blast/HEX2986I03_TREP9_nr.bln";
ok ( (-e $obs_file) , "Blast m8 test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE m8 OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @blast_m8_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@blast_m8_obs, \@blast_m8_exp,
	    "blast m8 file contains correct data");

#-----------------------------+
# TEST BLAST GFF OUTPUT       |
#-----------------------------+
diag ("Testing batch_blast GFF results ..");
my $gff_exp_file = "data/exp/HEX2986I03_blast.gff";
open (EXPECTED, "<".$exp_file);
my @blast_gff_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE m8 BLAST OUTPUT FILE EXISTS
my $gff_obs_file = $out_dir."HEX2986I03/gff/HEX2986I03_blast.gff";
ok ( (-e $obs_file) , "Blast gff test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE m8 OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @blast_gff_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@blast_gff_obs, \@blast_gff_exp,
	    "blast gff file contains correct data");

exit;
