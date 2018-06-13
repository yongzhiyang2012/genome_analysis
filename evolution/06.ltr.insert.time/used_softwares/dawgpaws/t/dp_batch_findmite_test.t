#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_findmite_test.pl - Test batch_findmite.pl        |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/17/2009                                       |
# UPDATED: 04/18/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_findmite.pl program.                      |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 5;

#-----------------------------+
# TESTING FINDMITE BIN        |
#-----------------------------+
# The batch_findmite program tries to use the env var FINDMITE_BIN
# otherwise it will assume the binary is inclue in the user PATH
diag("Testing findmite ENV options ...");
my $findmite_bin;
if ($ENV{FINDMITE_BIN}) {
    ok ( $findmite_bin=$ENV{FINDMITE_BIN},
	 "findmite binary defined in environment as: $findmite_bin") 
    }
else {
    ok ( $findmite_bin="FINDMITE1New.bin", 
	 "findmite binary not defined in environment\n".
	 "expecting findmite in path as: $findmite_bin");
}

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the findmite program ...");
# Basic command to see that ltr_finder is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_findmite_cmd = $findmite_bin." 2>&1";

# Get the result
my $basic_findmite_result = `$basic_findmite_cmd`;

# Test for the ltr version string that comes as part of the help message
like ($basic_findmite_result, qr/Usage/, "findmite is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Test running the batch_findmite.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

#batch_findmite.pl -i data/fasta -o /home/jestill/sandbox/ -c data/config/findmite_params.jcfg
my $batch_findmite_cmd = "batch_findmite.pl".
    " -i data/fasta".
    " -o $out_dir".
    " -c data/config/findmite_params.jcfg";

# This will exit zero when things work
ok ( system($batch_findmite_cmd)==0 , "batch_findmite.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing batch_findmite gff results ..");
my $exp_file = "data/exp/HEX3045G05_TA_12mite.gff";
open (EXPECTED, "<".$exp_file);
my @findmite_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX3045G05/gff/HEX3045G05_TA_12mite.gff";
ok ( (-e $obs_file) , "findmite GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @findmite_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@ltrseq_obs, \@ltrseq_exp,
	    "ltr_finder GFF file contains correct data");

exit;
