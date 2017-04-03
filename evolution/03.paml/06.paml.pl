use strict;
use warnings;

## created by Yongzhi Yang. 2017/3/20 ##

my @in=<align/*/cds.paml>;
my @ctl=<ctl/*ctl>;
for my $ctl (@ctl){
    for my $in (@in){
        $in=~/^(\S+)\/cds.paml/;
        print "cd $1 ; codeml ../../$ctl ; cd ../../\n";
    }
}
