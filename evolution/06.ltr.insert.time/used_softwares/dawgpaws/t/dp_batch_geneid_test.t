#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_batch_geneid_test.pl - Test batch_geneid.pl            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/16/2009                                       |
# UPDATED: 04/16/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_geneid.pl program.                        |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 5;

#-----------------------------+
# TESTING BINARY DEFINED      |
#-----------------------------+
my $geneid_bin;
ok ( $geneid_bin = $ENV{GENEID_BIN},
     "geneid binary defined") || 
    diag ("geneid binary should be defined with GENEID_BIN\n".
	  "otherwise will assume geneid in PATH");

#-----------------------------+
# TESTING BINARY PATH WORKS   |
#-----------------------------+
print STDERR "Geneid at: $geneid_bin\n";

# Basic command to see that geneid is installed
my $basic_geneid_cmd = $geneid_bin." -h";
system ($basic_geneid_cmd);
ok ( system($basic_geneid_cmd)==0, "geneid is installed" ) ||
    diag ("You must have geneid installed, and in your PATH ".
	  "or the path must be defined with GENEID_BIN");

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Testing running the batch_geneid.pl program ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $can_clean = 0;  # The tmp dir did not already exist
unless (-e $out_dir) {
    mkdir $out_dir;
}

my $batch_geneid_cmd = "batch_geneid.pl -i data/fasta/ -o $out_dir".
    " -p data/config/geneid_wheat.param";

# This will exit zero when things work
ok ( system($batch_geneid_cmd)==0 , "batch_geneid.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing geneid gff results ..");
my $exp_file = "data/exp/HEX3045G05.geneid.gff";
open (EXPECTED, "<".$exp_file);
my @geneid_exp = <EXPECTED>;
close (EXPECTED);


# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX3045G05/gff/HEX3045G05.geneid.gff";
ok ( (-e $obs_file) , "Gaps GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

# TEST THAT THE OUTPUT MATCHES THE EXPECTED RESULT
open (OBSERVED, "<".$obs_file);
my @geneid_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@geneid_obs, \@geneid_exp,
	    "geneid GFF file contains correct data");

exit;
