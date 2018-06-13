#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_repmask_test.pl - Test batch_trf.pl                  |
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

use Cwd;
use Test::More tests => 5;

#-----------------------------+
# TESTING TRF BIN             |
#-----------------------------+
# The batch_trf program tries to use the env var TRF_BIN
# otherwise it will assume the binary is inclue in the user PATH
diag("Testing RepeatMasker ENV options ...");
my $rm_bin;
if ($ENV{DP_RM_BIN}) {
    ok ( $rm_bin=$ENV{DP_RM_BIN},
	 "RepeatMasker binary defined in environment as: $rm_bin"); 
    }
else {
    ok ( $rm_bin="RepeatMasker", 
	 "RepeatMasker binary not defined in environment\n".
	 "expecting RepeatMasker in path as: $rm_bin");
}

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the RepeatMasker program ...");
# Basic command to see that RepeatMasker is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_rm_cmd = $rm_bin." 2>&1";

# Get the result
my $basic_rm_result = `$basic_rm_cmd`;

# Test for the ltr version string that comes as part of the help message
like ($basic_rm_result, qr/RepeatMasker version/, "RepeatMasker is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Test running the batch_repmask.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

# Get the working dir of the dp_batch_repmask_test.t program
# use this to set the dir for the repmask databases
my $cwd = getcwd($0);
my $db_dir = $cwd."/data/repmask/";
my $batch_rm_cmd = "batch_repmask.pl -i data/fasta".
    " -o $out_dir".
    " -c data/config/batch_mask.jcfg".
    " --rm-dir $db_dir".
    " --engine wublast";

# This will exit zero when things work
ok ( system($batch_rm_cmd)==0 , "batch_repmask.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing batch_repmask gff results ..");
my $exp_file = "data/exp/HEX3045G05_rm_TREP9.gff";
open (EXPECTED, "<".$exp_file);
my @rm_exp = <EXPECTED>;
close (EXPECTED);


# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX3045G05/gff/HEX3045G05_rm_TREP9.gff";
ok ( (-e $obs_file) , "RepeatMasker GFF test files appear to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @rm_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@rm_obs, \@rm_exp,
	    "RepeatMasker GFF file contains correct data");

exit;
