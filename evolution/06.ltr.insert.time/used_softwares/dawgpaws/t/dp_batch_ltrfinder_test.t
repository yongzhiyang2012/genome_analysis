#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_ltrfinder_test.pl - Test batch_geneid.pl         |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/16/2009                                       |
# UPDATED: 04/16/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_ltrfinder.pl program.                     |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 5;

#-----------------------------+
# TESTING BINARY DEFINED      |
#-----------------------------+
diag("Testing that ltr_finder is available ...");
my $ltrfinder_bin;
ok ( $ltrfinder_bin = $ENV{LTR_FINDER},
     "ltr_finder location defined as $ltrfinder_bin") || 
    diag ("location of ltr_finder should be defined with FIND_LTR_PATH".
	 " or make available as ltr_finder in your PATH");

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
# Basic command to see that ltr_finder is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_ltrfinder_cmd = $ltrfinder_bin." -h 2>&1";
# Get the result
my $basic_ltrfinder_result = `$basic_ltrfinder_cmd`;

# Test for the ltr version string that comes as part of the help message
like ($basic_ltrfinder_result, qr/ltr_finder v/, "ltr_finder is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Test running the batch_ltrfinder.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

my $batch_ltrfinder_cmd = "batch_ltrfinder.pl -i data/fasta/ -o $out_dir".
    " -c data/config/batch_ltrfinder.jcfg --gff";

# This will exit zero when things work
ok ( system($batch_ltrfinder_cmd)==0 , "batch_ltrfinder.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing ltrfinder gff results ..");
my $exp_file = "data/exp/HEX3045G05_def_ltr_finder.gff";
open (EXPECTED, "<".$exp_file);
my @ltrfinder_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX3045G05/gff/HEX3045G05_def_ltr_finder.gff";
ok ( (-e $obs_file) , "ltr_finder GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @ltrfinder_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@ltrfinder_obs, \@ltrfinder_exp,
	    "ltr_finder GFF file contains correct data");

exit;
