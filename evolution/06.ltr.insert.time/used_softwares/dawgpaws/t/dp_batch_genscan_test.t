#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_genscan_test.pl - Test batch_genscan.pl          |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/28/2009                                       |
# UPDATED: 04/28/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_genscan.pl program.                       |
#                                                           |
#-----------------------------------------------------------+

use Cwd;
use Test::More tests => 6;

#-----------------------------+
# TESTING TENEST BIN          |
#-----------------------------+
diag("Testing genscan ENV options ...");
my $genscan_bin;
if ($ENV{'DP_GENSCAN_BIN'}) {
    ok ( $genscan_bin=$ENV{'DP_GENSCAN_BIN'},
	 "genscan binary defined in environment as: $genscan_bin") 
    }
else {
    ok ( $genscan_bin="genscan", 
	 "genscan binary not defined in environment\n".
	 "expecting genscan in PATH as: $genscan_bin");
}

#-----------------------------+
# TESTING LIB PATH            |
#-----------------------------+
diag("Testing genscan ENV options ...");
my $lib_path;
if ($ENV{'DP_GENSCAN_LIB'}) {
    ok ( $lib_path=$ENV{'DP_GENSCAN_LIB'},
	 "genscan library defined in environment as: $lib_path") 
    }
else {
    ok ( $lib_path="Maize.smat", 
	 "genscan library not defined in environment\n".
	 "default library set to: $lib_path");
}

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the genscan program ...");
# Basic command to see that genscan is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_genscan_cmd = $genscan_bin." -h 2>&1";

# Get the result
my $basic_genscan_result = `$basic_genscan_cmd`;

# Test that TEnest.pl is working by finding the 'gimme some commands'
# response when running TEnest without command options
like ($basic_genscan_result, qr/usage: genscan/, "genscan is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+

diag ("Test running the batch_genscan.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

# Get the working dir of the dp_batch_repmask_test.t program
# use this to set the dir for the hmm databases
my $cwd = getcwd($0);

#batch_genscan.pl -i data/fasta -o /home/jestill/sandbox/
my $batch_genscan_cmd = "batch_genscan.pl -i data/fasta".
    " -o $out_dir";
#    " --lib-path Maize.smat";

# This will exit zero when things work
ok ( system($batch_genscan_cmd)==0 , "batch_genscan.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing batch_genscan gff results ..");
my $exp_file = "data/exp/HEX2903P03.genscan.gff";
open (EXPECTED, "<".$exp_file);
my @genscan_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX2903P03/gff/HEX2903P03.genscan.gff";
ok ( (-e $obs_file) , "batch_genscan GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @genscan_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@genscan_obs, \@genscan_exp,
	    "genscan GFF file contains correct data");

exit;
