use strict;
use warnings;

## created by Yongzhi Yang. 2017/3/20 ##

my $Gblocks="/home/share/users/yangyongzhi2012/tools/gblocks/Gblocks_0.91b/Gblocks";
my $tools="scripts/Gblocks2Paml.pl";
my @in=<align/*/cds.best.fas>;
for my $in (@in){
    my $out=$in;
    $out=~s/cds.best.fas$/cds.paml/;
    print "$Gblocks $in -t=c ; perl $tools $in-gb $out\n";
}

