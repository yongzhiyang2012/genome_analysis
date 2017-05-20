use strict;
use warnings;

## created by Yongzhi Yang. 2017/3/20 ##

my %group;
my %sp;
open (F,"groups.txt")||die"$!";
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    for (my $i=1;$i<@a;$i++){
        $a[$i]=~/^(\w+)\|/;
        my $id=$1;
        $sp{$id}++;
        $group{$a[0]}{$id}{$a[$i]}++;
    }
}
close F;
my $spnum=scalar(keys %sp);
my %genenum;
open (F,"compliantFasta/genenum.txt")||die"$!";
while (<F>) {
    chomp;
    /(\w+)\.fasta\:(\d+)/;
    $genenum{$1}=$2;
}
close F;
my %count;
for my $k1 (sort keys %group){
    my @k2=sort keys %{$group{$k1}};
    my $type1="paralogs";
    if (scalar(@k2) == $spnum){
        $type1="orthologs";
    }elsif(scalar(@k2) > 1){
        $type1="other";
    }elsif(scalar(@k2) == 1){
        $type1="paralogs";
    }else{
        die "wrong type\n";
    }
    for my $k2 (@k2){
        my @k3=sort keys %{$group{$k1}{$k2}};
        my $type2="single";
        $type2="multiple" if (scalar(@k3) > 1);
        $type2=$type1 if $type1 eq 'other';
        $type2=$type1 if $type1 eq 'paralogs';
        $count{$k2}{$type1}{$type2} += scalar(@k3);
        $count{$k2}{all} += scalar(@k3);
    }
}

open (O,">$0.plot");
print O "ID\ttype\tnum\n";
for my $k1 (sort keys %count){
    my $unique=$genenum{$k1} - $count{$k1}{all};
    #print O "$k1\tUnclustered_genes\t$unique\n";
    $count{$k1}{orthologs}{single}=0 if ! exists $count{$k1}{orthologs}{single};
    print O "$k1\tSingle_copy_orthologs\t$count{$k1}{orthologs}{single}\n";
    $count{$k1}{orthologs}{multiple}=0 if ! exists $count{$k1}{orthologs}{multiple};
    print O "$k1\tMultiple_copy_orthologs\t$count{$k1}{orthologs}{multiple}\n";
    $count{$k1}{paralogs}{paralogs}=0 if ! exists $count{$k1}{paralogs}{paralogs};
    print O "$k1\tUnique_paralogs\t$count{$k1}{paralogs}{paralogs}\n";
    $count{$k1}{other}{other}=0 if ! exists $count{$k1}{other}{other};
    print O "$k1\tOther_orthologs\t$count{$k1}{other}{other}\n";
    print O "$k1\tUnclustered_genes\t$unique\n";
}
close O;

    
