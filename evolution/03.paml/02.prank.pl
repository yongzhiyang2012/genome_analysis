use strict;
use warnings;

## created by Yongzhi Yang. 2017/3/20 ##

my $tool="/home/share/users/yangyongzhi2012/tools/prank/prank/bin/prank";

my @in=<align/*/cds>;
for my $in (@in){
    print "$tool -d=$in -o=$in -f=fasta -codon +F\n";
}
