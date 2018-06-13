use strict;
use warnings;

my $RepeatModeler=shift or die "perl $0 <path-to-RepeatModeler>;\nthe path should be the dir contained RepeatClassifier\n";

my $in="genome.ltrfinder.seq.fa";
die "give the right inpu genome.ltrfinder.seq.fa\n" if ! -e "$in";
my $pwd=`pwd`;
chomp $pwd;
open (O,">$0.sh");

`mkdir $in.split` if (! -e "$in.split");
chdir("$in.split");
`$pwd/used_softwares/fasta-splitter.pl --n-parts 100 ../$in`;

while (<./*part-*>) {
    chomp;
    if ($_=~/^(\S+\.part-(\d+))/){
        `mkdir $2 ; mv $1* $2 `;
        print O "cd $in.split/$2 ; $RepeatModeler/RepeatClassifier -engine abblast -consensi $1.fa ; cd ../../\n";
    }
}

chdir("../");
close O;
