#!/usr/bin/perl -w
#-----------------------------------------------------------+
#                                                           |
# dp_venn_test.pl - Test the vennseq.pl program.            |
#                                                           |
#-----------------------------------------------------------+
#  AUTHOR: James C. Estill                                  |
# CONTACT: jestill_at_sourceforge.net                       |
# STARTED: 04/14/2009                                       |
# UPDATED: 04/14/2009                                       |
#                                                           |
# SHORT DESCRIPTION:                                        |
#  Test the vennseq.pl program. This tests the basic        |
#  crosstab results, as well as checks that the vennmaster  |
#  program is installed and runs correctly.                 |
#-----------------------------------------------------------+

use Test::More tests => 7;

#-----------------------------+
# PROGRAM CROSS TAB           |
#-----------------------------+
diag ("Testing vennseq program crosstab output ...");
open (EXPECTED, '<data/exp/venn_test_result.txt');
my @vennseq_exp = <EXPECTED>;
close (EXPECTED);

my @vennseq_obs = `vennseq.pl -d data/venn_test/ -i data/venn_test/HEX3045G05_TREP9.masked.fasta -o ~/test_venn.txt --no-vm`;

is_deeply ( \@vennseq_obs, \@vennseq_exp, "vennseq.pl");

#-----------------------------+
# CHECKING ENV OPTIONS        |
# VARS WITH EVN OPTIONS       |
#-----------------------------+
diag ("Checking vennseq ENV options ...");

# VennMaster Directory Defined
ok ( my $vmaster_dir = $ENV{VMASTER_DIR} 
     || "/Applications/VennMaster-0.37.3/" , 
     "VennMaster Dir Defined"  );
ok ( -e $vmaster_dir, "VennMaster Dir Exits");

# VennMaster Java Binary Defined
my $java_bin;
ok ( $java_bin = $ENV{VMASTER_JAVA_BIN} || 'java', 
     "VennMaster Java Binary Defined as $java_bin" );

# Java binary valid
my $java_test_cmd = "$java_bin -showversion";
ok ( my $java_result = system ($java_test_cmd) , 
     "Java Binary Valid");

# Java memory okay
my $java_mem;
ok ( $java_mem = $ENV{VMASTER_JAVA_MEM} || '512', 
     "VennMaster Java Memory Allotment Defined as $java_mem" );

#-----------------------------+
# CHECK FOR VMASTER PROGRAM   |
#-----------------------------+
my $vm_loc = "$vmaster_dir"."venn.jar";
my $vm_result = `$java_bin -Xms256m -Xmx256m -jar $vm_loc --version`;
ok ($vm_result =~ m/VennMaster version(.*)/, 
    "VennMaster Working Version: $1");

