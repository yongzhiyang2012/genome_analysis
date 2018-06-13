#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_cnv_gff3_test.t - Test convert from game.xml to gff3   |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/16/2009                                       |
# UPDATED: 04/16/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test of the conversion from game.xml to gff3 format.     |
#                                                           |
#-----------------------------------------------------------+

use Test::More tests => 4;

diag ("Testing Apollo dependent game.xml converters ...");

#-----------------------------+
# TESTING APOLLO BIN DEFINED  |
#-----------------------------+
my $ap_bin;
ok ( $ap_bin = $ENV{DP_APOLLO_BIN} || "apollo", 
     "Apollo Binary Path Defined as $ap_bin" );

#-----------------------------+
# APOLLO TEST COMMAND         |
#-----------------------------+
my $ap_test_cmd = "$ap_bin --help";
ok ( my $apollo_result = `$ap_test_cmd`,
    "Apollo test command");

#-----------------------------+
# cnv_game2gff3.pl            |
#-----------------------------+
diag ("Testing apollo mediated conversion ..");
open (EXPECTED, '<data/exp/HEX3045G05_gff3_withseq.gff');
my @gff3_exp = <EXPECTED>;
close (EXPECTED);

my $tmp_out_dir = $ENV{HOME}."/";
my $tmp_out_file = $tmp_out_dir."hex_gff3test.gff";

my $num_lines = @game_xml_expected;

# OBSERVED
my $cnv_cmd = "cnv_game2gff3.pl".
    " -i data/HEX3045G05/HEX3045G05_ann_test.game.xml".
    " -o $tmp_out_file";
# Need to make the cnv_game2gff3.pl program return true
# When everything works, this will exit 0 .. no errors
ok ( system($cnv_cmd)==0 , "cnv_game2gff3.pl - run");
#system ($cnv_cmd);

my $obs_file = $tmp_out_file;
open (OBSERVED, '<'.$obs_file) ||
    die "Can not open the result at $obs_file\n";
my @gff3_obs = <OBSERVED>;
close (OBSERVED);

is_deeply ( \@gff3_obs, \@gff3_exp,
	    "cnv_game2gff3.pl parse - result expected");
