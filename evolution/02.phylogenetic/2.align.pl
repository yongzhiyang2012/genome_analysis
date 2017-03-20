use strict;
use warnings;

## created by Yongzhi Yang. 2017/3/20 ##

my $muscle="/home/share/users/yangyongzhi2012/tools/Muscle/muscle";
my $pal2nal="/home/share/users/yangyongzhi2012/tools/pal2nal/pal2nal.v14/pal2nal.pl";

my @in=<align/*/pep>;
for my $in (@in){
    my $in2=$in;
    $in2=~s/pep/cds/;
    print "$muscle -in $in -out $in.best.fas ; $pal2nal $in.best.fas $in2 -output fasta > $in2.best.fas\n";
}
