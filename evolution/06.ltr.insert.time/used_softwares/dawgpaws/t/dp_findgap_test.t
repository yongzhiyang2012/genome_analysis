#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_findgap_test.pl - Test the findgap program.            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/16/2009                                       |
# UPDATED: 04/16/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the batch_findgaps.pl program.                      |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 3;

#-----------------------------+
# TEST RUNNING THE PROGRAM    |
#-----------------------------+
diag ("Testing running the batch_findgaps.pl program, this is slow ...");
my $out_dir = $ENV{HOME}."/dp_temp_test/";
my $findgaps_cmd = "batch_findgaps.pl -i data/fasta/ -o $out_dir";

# This will exit zero when things work
ok ( system($findgaps_cmd)==0 , "batch_findgaps.pl");

#-----------------------------+
# TEST PROGRAM OUTPUT         |
#-----------------------------+
diag ("Testing gap results ..");
my $exp_file = "data/exp/HEX2903P03_gaps.gff";

open (EXPECTED, "<".$exp_file);
my @gaps_exp = <EXPECTED>;
close (EXPECTED);

# TEST THE THE GFF OUTPUT FILE EXISTS
my $obs_file = $out_dir."HEX2903P03/gff/HEX2903P03_gaps.gff";
ok ( (-e $obs_file) , "Gaps GFF test files appears to exist") ||
    diag("I expected to see the file\n $obs_file");

open (OBSERVED, "<".$obs_file);
my @gaps_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@gaps_obs, \@gaps_exp,
	    "Gaps GFF file contains correct data");
