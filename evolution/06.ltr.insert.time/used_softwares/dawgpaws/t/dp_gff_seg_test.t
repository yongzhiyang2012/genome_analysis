#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_gff_seg_test.pl - Testing gff segmentation program     |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/15/2009                                       |
# UPDATED: 04/15/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test of gff segmentation program.                        |
#-----------------------------------------------------------+

use Test::More tests => 1;

#-----------------------------+
# PROGRAM CROSS TAB           |
#-----------------------------+
diag ("Testing vennseq program crosstab output ...");
open (EXPECTED, '<data/exp/gff_seg_50x_expected.txt');
my @gffseg_exp = <EXPECTED>;
close (EXPECTED);

my @gffseg_obs = `gff_seg.pl -i data/oligo_count/gff_seg_test_in.gff --thresh 50`;

is_deeply ( \@gffseg_obs, \@gffseg_exp, "gff_seg.pl");
