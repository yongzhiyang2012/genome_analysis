#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_ltrseq_test.pl - Test batch_ltrseq.pl            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/17/2009                                       |
# UPDATED: 04/17/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_ltrfinder.pl program.                     |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 6;

#-----------------------------+
# TESTING CONFIG FILE DIR     |
#-----------------------------+
diag("Testing ltr_seq ENV options ...");
my $ltrseq_dir;
# The batch_ltrseq program tries to use the env var LTR_SEQ_DIR
# otherwise it will 
ok ( $ltrseq_dir=$ENV{LTR_SEQ_DIR} || $ENV{HOME},
     "ltr_seq config file directory set to: $ltrseq_dir") ||
    diag ("The directory holding the ltr_seq config files should be\n".
	  "defined with LTR_SEQ_DIR");

my $ltrseq_bin;
ok ( $ltrseq_bin = $ENV{LTR_SEQ_BIN} || "LTR_seq",
     "ltr_seq location defined as: $ltrseq_bin") || 
    diag ("location of LTR_seq program should be defined with FIND_LTR_BIN".
	 " or make available as LTR_seq in your PATH");

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the ltr_seq program ...");
# Basic command to see that ltr_finder is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_ltrseq_cmd = $ltrseq_bin." 2>&1";

# Get the result
my $basic_ltrseq_result = `$basic_ltrseq_cmd`;

# Test for the ltr version string that comes as part of the help message
like ($basic_ltrseq_result, qr/Usage : LTR_seq/, "ltr_seq is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Test running the batch_ltrseq.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

my $batch_ltrseq_cmd = "batch_ltrseq.pl -i data/fasta".
    " -o $out_dir".
    " --config-dir data/config".
    " -c data/config/batch_ltrseq.jcfg";

# This will exit zero when things work
ok ( system($batch_ltrseq_cmd)==0 , "batch_ltrseq.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing ltr_seq gff results ..");
my $exp_file = "data/exp/HEX3045G05_ltrseq_default.gff";
open (EXPECTED, "<".$exp_file);
my @ltrseq_exp = <EXPECTED>;
close (EXPECTED);


# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX3045G05/gff/HEX3045G05_ltrseq_default.gff";
ok ( (-e $obs_file) , "ltr_seq GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @ltrseq_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@ltrseq_obs, \@ltrseq_exp,
	    "ltr_finder GFF file contains correct data");

exit;
