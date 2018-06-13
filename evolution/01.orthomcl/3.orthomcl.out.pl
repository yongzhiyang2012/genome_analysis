use strict;
use warnings;

my $in="groups.txt";
open (F,"$in")||die"$!";
open (O,">AllClusters.txt");
print O "#Fam_id\tGenes\tSpecies\tOrthologous genes\n";
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    my $Fam_id=shift(@a);
    $Fam_id=~s/\://;
    my %h;
    for my $id (@a){
        $id=~/^(\S+)\|\S+$/;
        $h{$1}++;
    }
    print O "$Fam_id\t",scalar(@a),"\t",scalar(keys %h),"\t",join(",",@a),"\n";
}
close F;
close O;
