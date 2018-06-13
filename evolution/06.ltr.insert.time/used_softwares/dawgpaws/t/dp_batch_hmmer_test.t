#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_hmmer_test.pl - Test batch_hmmer.pl              |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/24/2009                                       |
# UPDATED: 04/27/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_hmmer.pl program.                         |
#                                                           |
#-----------------------------------------------------------+

use Cwd;
use Test::More tests => 5;

#-----------------------------+
# TESTING HMMSERACH BIN       |
#-----------------------------+
diag("Testing hmmsearch ENV options ...");
my $hmm_bin;
if ($ENV{DP_HMMSEARCH_BIN}) {
    ok ( $hmm_bin=$ENV{DP_HMMSEARCH_BIN},
	 "hmmsearch binary defined in environment as: $hmm_bin") 
    }
else {
    ok ( $hmm_bin="hmmsearch", 
	 "hmmsearch binary not defined in environment\n".
	 "expecting hmmsearch in PATH as: $hmm_bin");
}

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the hmmsearch program ...");
# Basic command to see that RepeatMasker is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_hmm_cmd = $hmm_bin." -h 2>&1";

# Get the result
my $basic_hmm_result = `$basic_hmm_cmd`;

# Test for the ltr version string that comes as part of the help message
like ($basic_hmm_result, qr/HMMER/, "hmmsearch is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+


#batch_hmmer.pl -i fasta -o sandbox/ -c config/batch_hmmer.jcfg


diag ("Test running the batch_hmmer.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}


# Get the working dir of the dp_batch_repmask_test.t program
# use this to set the dir for the hmm databases
my $cwd = getcwd($0);
my $db_dir = $cwd."/data/";

my $batch_hmm_cmd = "batch_hmmer.pl -i data/fasta".
    " -o $out_dir".
    " --db-parent $db_dir".
    " -c data/config/batch_hmmer.jcfg";

# This will exit zero when things work
ok ( system($batch_hmm_cmd)==0 , "batch_hmmer.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing batch_hmmer gff results ..");
my $exp_file = "data/exp/HEX2903P03_hmm_dp_test.gff";
open (EXPECTED, "<".$exp_file);
my @hmmer_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX2903P03/gff/HEX2903P03_hmm_dp_test.gff";
ok ( (-e $obs_file) , "batch_hmmer GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");


# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @hmmer_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@hmmer_obs, \@hmmer_exp,
	    "RepeatMasker GFF file contains correct data");

exit;
