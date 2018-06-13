#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_eugene_test.pl - Test batch_eugene.pl            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/16/2009                                       |
# UPDATED: 04/16/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_eugene.pl program.                        |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 4;

#-----------------------------+
# TESTING BINARY PATH         |
#-----------------------------+
my $eugene_bin = "eugene";
my $get_sites =  "egn_getsites4eugene.pl";
my $param_file = "hex_eugene.par";

# Basic command to see that eugene is installed
my $basic_eug_cmd = "eugene";
ok ( system($basic_eug_cmd), "Eugene is Installed" ) ||
    diag ("You must have eugene installed, and in your PATH.");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Testing running the batch_eugene.pl program, this is very slow ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";

my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

my $batch_eug_cmd = "batch_eugene.pl -i data/fasta/ -o $out_dir".
    " -p data/config/hex_eugene.par";

# This will exit zero when things work
ok ( system($batch_eug_cmd)==0 , "batch_eugene.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing batch_eugene results ..");
my $exp_file = "data/exp/HEX3045G05_eugene.gff";

open (EXPECTED, "<".$exp_file);
my @eugene_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX3045G05/gff/HEX3045G05_eugene.gff";
ok ( (-e $obs_file) , "Gaps GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE RESULTS MEET EXPECTATIONS
open (OBSERVED, "<".$obs_file);
my @eugene_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@eugene_obs, \@eugene_exp,
	    "Gaps GFF file contains correct data");

exit;
