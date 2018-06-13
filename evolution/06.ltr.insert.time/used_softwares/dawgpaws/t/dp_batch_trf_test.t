#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_trf_test.pl - Test batch_trf.pl                  |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/18/2009                                       |
# UPDATED: 04/23/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_trf.pl program.                           |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 5;

#-----------------------------+
# TESTING TRF BIN             |
#-----------------------------+
# The batch_trf program tries to use the env var TRF_BIN
# otherwise it will assume the binary is inclue in the user PATH
diag("Testing trf ENV options ...");
my $trf_bin;
if ($ENV{TRF_BIN}) {
    ok ( $trf_bin=$ENV{TRF_BIN},
	 "trf binary defined in environment as: $trf_bin");
    }
else {
    ok ( $trf_bin="trf400.linux.exe", 
	 "trf binary not defined in environment\n".
	 "expecting trf in path as: $trf_bin");
}

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the trf program ...");
# Basic command to see that ltr_finder is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_trf_cmd = $trf_bin." 2>&1";

# Get the result
my $basic_trf_result = `$basic_trf_cmd`;

# Test for the ltr version string that comes as part of the help message
like ($basic_trf_result, qr/Tandem Repeats Finder/, "trf is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Test running the batch_trf.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

my $batch_trf_cmd = "batch_trf.pl".
    " -i data/fasta".
    " -o $out_dir";

# This will exit zero when things work
ok ( system($batch_trf_cmd)==0 , "batch_trf.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing batch_trf gff results ..");
my $exp_file = "data/exp/HEX3045G05_trf.gff";
open (EXPECTED, "<".$exp_file);
my @trf_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX3045G05/gff/HEX3045G05_trf.gff";
ok ( (-e $obs_file) , "trf GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");


# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @trf_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@trf_obs, \@trf_exp,
	    "trf GFF file contains correct data");

exit;
