#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_cnv_test.pl - Test DAWGPAWS conversion programs.       |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/14/2009                                       |
# UPDATED: 04/14/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test of the DAWGPAWS conversion programs.                |
#-----------------------------------------------------------+

use Test::More tests => 10;

diag ("Testing DAWGPAWS converters ...");

#-----------------------------+
# FGENESH TXT                 |
#-----------------------------+
open (EXPECTED, '<data/exp/fgenesh_cnv_expected.txt');
my @fgenesh_exp = <EXPECTED>;
close (EXPECTED);

my @fgenesh_obs =  `cnv_fgenesh2gff.pl -i data/HEX3045G05/fgenesh/fgenesh.txt`;

# Test the observed fgenesh parser result against the expected result
is_deeply ( \@fgenesh_obs, \@fgenesh_exp, 
	    "cnv_fgenesh2gff.pl text conversion" );

#-----------------------------+
# FGENESH HTML                |
#-----------------------------+
open (EXPECTED, '<data/exp/fgenesh_cnv_expected.txt');
my @fgenesh_html_exp = <EXPECTED>;
close (EXPECTED);

my @fgenesh_html_obs =  `cnv_fgenesh2gff.pl -i data/HEX3045G05/fgenesh/fgenesh_html.txt --html`;

# Test the observed fgenesh parser result against the expected result
is_deeply ( \@fgenesh_html_obs, \@fgenesh_html_exp, 
	    "cnv_fgenesh2gff.pl html conversion" );

#-----------------------------+
# GENEMARK                    |
#-----------------------------+
open (EXPECTED, '<data/exp/genemark_cnv_expected.txt');
my @genemark_exp = <EXPECTED>;
close (EXPECTED);

my @genemark_obs = `cnv_genemark2gff.pl -i data/HEX3045G05/genemark/HEX3045G05_genemark_ta.out`;

is_deeply ( \@genemark_obs, \@genemark_exp, "cnv_genemark2gff.pl");

#-----------------------------+
# FIND LTR                    |
#-----------------------------+
open (EXPECTED, '<data/exp/findltr_cnv_expected.txt');
my @findltr_exp = <EXPECTED>;
close (EXPECTED);

my @findltr_obs = `cnv_findltr2gff.pl -i data/HEX3045G05/find_ltr/HEX3045G05_findltr_def.txt`;

is_deeply ( \@findltr_obs, \@findltr_exp, "cnv_findltr2gff.pl");

#-----------------------------+
# LTR FINDER                  |
#-----------------------------+
open (EXPECTED, '<data/exp/ltrfinder_cnv_expected.txt');
my @ltrfinder_exp = <EXPECTED>;
close (EXPECTED);

my @ltrfinder_obs = `cnv_ltrfinder2gff.pl -i data/HEX3045G05/ltr_finder/HEX3045G05_def_ltr_finder.txt`;

is_deeply ( \@ltrfinder_obs, \@ltrfinder_exp, "cnv_ltrfinder2gff.pl");

#-----------------------------+
# LTR_seq                     |
#-----------------------------+
open (EXPECTED, '<data/exp/ltrseq_cnv_expected.txt');
my @ltrseq_exp = <EXPECTED>;
close (EXPECTED);

my @ltrseq_obs = `cnv_ltrseq2gff.pl -i data/HEX3045G05/ltr_seq/HEX3045G05_ltrseq.txt -s HEX304G05`;

is_deeply ( \@ltrseq_obs, \@ltrseq_exp, "cnv_ltrseq2gff.pl" );

#-----------------------------+
# REPSEEK                     |
#-----------------------------+
open (EXPECTED, '<data/exp/repseek_cnv_expected.txt');
my @repseek_exp = <EXPECTED>;
close (EXPECTED);

my @repseek_obs = `cnv_repseek2gff.pl -i data/HEX3045G05/repseek/repseek_l20.txt`;

is_deeply ( \@repseek_obs, \@repseek_exp, "cnv_repseek2gff.pl");

#-----------------------------+
# BLAST                       |
#-----------------------------+
open (EXPECTED, '<data/exp/blast_cnv_expected.txt');
my @blast_exp = <EXPECTED>;
close (EXPECTED);

my @blast_obs = `cnv_blast2gff.pl -i data/HEX3045G05/blast/HEX3045G05_AcTA_1.bln`;

is_deeply ( \@blast_obs, \@blast_exp, "cnv_blast2gff.pl" );

#-----------------------------+
# REPEATMASKER                |
#-----------------------------+
open (EXPECTED, '<data/exp/repeatmasker_cnv_expected.txt');
my @repmask_exp = <EXPECTED>;
close (EXPECTED);

my @repmask_obs = `cnv_repmask2gff.pl -i data/HEX3045G05/rm/HEX3045G05_TREP9.rm.out`;

is_deeply ( \@repmask_obs, \@repmask_exp, "cnv_repmask2gff.pl" );

#-----------------------------+
# TE NEST                     |
#-----------------------------+
open (EXPECTED, '<data/exp/tenest_cnv_expected.txt');
my @tenest_exp = <EXPECTED>;
close (EXPECTED);

my @tenest_obs = `cnv_tenest2gff.pl -i data/HEX3045G05/tenest/HEX3045G05.LTR`;

is_deeply ( \@tenest_obs, \@tenest_exp, "cnv_tenest2gff.pl" );

