use strict;
use warnings;

## created by Yongzhi Yang. 2017/3/20 ##


open (F,"groups.txt")||die"$!";
my %count;
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    $a[0]=~s/\://;
    for (my $i=1;$i<@a;$i++){
        $a[$i]=~/^(\w+)\|/;
        my $id=$1;
        next unless($id=~/Omu|Ore|Ppe|Jre|Fve/);
        $count{$id}{$a[0]}++;
    }
}
close F;
for my $k (sort keys %count){
    open (O,">$0.$k.list");
    for my $k2 (sort keys %{$count{$k}}){
        print O "$k2\n";
    }
    close O;
}
