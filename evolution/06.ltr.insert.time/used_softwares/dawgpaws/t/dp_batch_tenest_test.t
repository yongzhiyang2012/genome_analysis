#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_hmmer_test.pl - Test batch_hmmer.pl              |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/28/2009                                       |
# UPDATED: 04/28/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_tenest.pl program.                        |
#                                                           |
#-----------------------------------------------------------+

use Cwd;
use Test::More tests => 5;

#-----------------------------+
# TESTING TENEST BIN          |
#-----------------------------+
diag("Testing TEnest ENV options ...");
my $tenest_bin;
if ($ENV{'TE_NEST_BIN'}) {
    ok ( $tenest_bin=$ENV{'TE_NEST_BIN'},
	 "TEnest binary defined in environment as: $tenest_bin") 
    }
else {
    ok ( $tenest_bin="TEnest.pl", 
	 "TEnest binary not defined in environment\n".
	 "expecting hmmsearch in PATH as: $tenest_bin");
}

#-----------------------------+
# TESTING TENEST DATA DIR     |
#-----------------------------+
#my $tenest_dir;
#if ($ENV{'TE_NEST_DIR'}) {
#    ok ( $tenest_dir=$ENV{'TE_NEST_DIR'},
#	 "TEnest database directory defined in environment as: $tenest_dir") 
#    }
#else {
#    ok ( $tenest_dir=$ENV{'HOME'}, 
#	 "TEnest binary not defined in environment\n".
#	 "expecting databases in home directory as: $tenest_dir");
#}

#-----------------------------+
# TESTING TENEST WUBLAST DIR  |
#-----------------------------+
#my $wublast_dir;
#if ($ENV{'DP_WUBLAST_DIR'}) {
#    ok ( $wublast_dir=$ENV{'TE_NEST_DIR'},
#	 "Wublast directory defined in environment as: $wublast_dir") 
#    }
#else {
#    ok ( $wublast_dir="/usr/local/genome/wu_blast/", 
#	 "Wubalst directory not defined in environment\n".
#	 "expecting wublast at: /usr/local/genome/wu_blast");
#}


#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
diag ("Test running the TEnest.pl program ...");
# Basic command to see that RepeatMasker is installed
# Have to do a trick here to capture STDERR by sending it to STDOUT
my $basic_tenest_cmd = $tenest_bin." -h 2>&1";

# Get the result
my $basic_tenest_result = `$basic_tenest_cmd`;

# Test that TEnest.pl is working by finding the 'gimme some commands'
# response when running TEnest without command options
like ($basic_tenest_result, qr/gimme some commands/, "TEnest.pl is working");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+

# batch_tenest.pl -i data/fasta -o /home/jestill/sandbox/ --org /home/jestill/svnloc/dawg-paws/trunk/t/data/WHEAT_TEST

diag ("Test running the batch_tenest.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

# Get the working dir of the dp_batch_repmask_test.t program
# use this to set the dir for the hmm databases
my $cwd = getcwd($0);
my $db_dir = $cwd."/data/WHEAT_TEST";

my $batch_tenest_cmd = "batch_tenest.pl -i data/fasta".
    " -o $out_dir".
    " --org ".$db_dir;

# This will exit zero when things work
ok ( system($batch_tenest_cmd)==0 , "batch_tenest.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing batch_tenest gff results ..");
my $exp_file = "data/exp/HEX2903P03.tenest.gff";
open (EXPECTED, "<".$exp_file);
my @tenest_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX2903P03/gff/HEX2903P03.tenest.gff";
ok ( (-e $obs_file) , "batch_tenest GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @tenest_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@tenest_obs, \@tenest_exp,
	    "TEnest GFF file contains correct data");

exit;
