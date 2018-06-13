#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_test_all.pl - Run all DAWGPAWS tests.                  |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/15/2009                                       |
# UPDATED: 04/15/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Run all of the test scripts for the DAWGPAWS programs.   |
#-----------------------------------------------------------+

use TAP::Harness;

print STDERR "Running All DAWGPAWS tests\n";

my %args = (verbosity =>1, color=>1);

my $harness = TAP::Harness->new( \%args );

$harness->runtests( 
		    ['dp_module_test.t', 'TESTING DAWGPAWS REQUIRED MODULES'],
		    ['dp_env_test.t', 'TESTING DAWGPAWS ENVIRONMENT'],
		    ['dp_cnv_test.t', 'TESTING DAWGPAWS CONVERTER SCRIPTS'],
		    ['dp_cnv_game_test.t', 'TESTING GFF2GAME CONVERTER'],
		    ['dp_venn_test.t', 'TESTING VENN PROGRAM'],
		    ['dp_gff_seg_test.t', 'TESTING GFF SEGMENTATION'],
		    ['dp_batch_blast_test.t', 'TESTING BATCHBLAST'],
		    ['dp_batch_eugene_test.t', 'TESTING BATCH EUGENE'],
		    ['dp_batch_findltr_test.t', 'TESTING BATCH FINDLTR'],
		    ['dp_batch_findmite_test.t', 'TESTING BATCH FINDMITE'],
		    ['dp_batch_geneid_test.t', 'TESTING BATCH GENEID'],
		    ['dp_batch_hmmer_test.t', 'TESTING BATCH HMMER'],
		    ['dp_batch_ltrfinder_test.t', 'TESTING BATCH LTRFINDER'],
		    ['dp_batch_ltrseq_test.t', 'TESTING BATCH LTRSEQ'],
		    ['dp_batch_repmask_test.t', 'TESTING BATCH REPMASK'],
		    ['dp_batch_tenest_test.t', 'TESTING BATCH TENEST'],
		    ['dp_batch_trf_test.t', 'TESTING BATCH TRF'],
		    ['dp_cnv_gff3_test.t', 'TESTING GFF3 CONVERSION'],
		    ['dp_findgap_test.t', 'TESTING BATCH FINDGAPS'],
		    ['dp_batch_genscan_test.t', 'TESTING BATCH GENSCAN'],
    );
