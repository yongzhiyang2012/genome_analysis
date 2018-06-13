#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_findltr_test.pl - Test batch_geneid.pl           |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/16/2009                                       |
# UPDATED: 04/16/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_findltr.pl program.                       |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 5;

#-----------------------------+
# TESTING BINARY DEFINED      |
#-----------------------------+
diag("Testing that find_ltr.pl is available");
my $geneid_bin;
my $findltr_bin;
ok ( $findltr_bin = $ENV{FIND_LTR_PATH},
     "find_ltr defined") || 
    diag ("location of find_ltr.pl should be defined with FIND_LTR_PATH");


#-----------------------------+
# TESTING BINARY PATH EXISTS  |
#-----------------------------+
# Basic command to see that geneid is installed
ok ( (-e $findltr_bin), "find_ltr.pl is installed at expected location" ) ||
    diag ("find_ltr.pl does not exist at $findltr_bin");


#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Test running the batch_findltr.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

my $batch_findltr_cmd = "batch_findltr.pl -i data/fasta/ -o $out_dir".
    " -c data/config/batch_findltr.jcfg";

# This will exit zero when things work
ok ( system($batch_findltr_cmd)==0 , "batch_findltr.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing geneid gff results ..");
my $exp_file = "data/exp/HEX2903P03_findltr_Def.gff";
open (EXPECTED, "<".$exp_file);
my @findltr_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX2903P03/gff/HEX2903P03_findltr_Def.gff";
ok ( (-e $obs_file) , "find_ltr GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @findltr_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@findltr_obs, \@findltr_exp,
	    "find_ltr GFF file contains correct data");

exit;
