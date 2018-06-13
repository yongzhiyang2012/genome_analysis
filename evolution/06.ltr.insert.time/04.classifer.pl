use strict;
use warnings;

my $RepeatModeler=shift or die "perl $0 <path-to-RepeatModeler> <consensi-lib-file>;\nthe path should be the dir contained RepeatClassifier\n";
my $consensi=shift or die "perl $0 <path-to-RepeatModeler> <consensi-lib-file>;\nthe path should be the dir contained RepeatClassifier\n";

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
        print O "cd $in.split/$2 ; $RepeatModeler/RepeatClassifier -consensi $consensi -engine abblast $1 ; cd ../../\n";
    }
}

chdir("../");
close O;
