use strict;
use warnings;

my %h;
my %max;
open (F,"compliantFasta.genenum.txt");
while (<F>) {
    chomp;
    /^(\w+)\.fasta\:(\d+)$/ or die "$_\n";
    $h{$1}{gene}=$2;
}
close F;

open (F,"groups.txt");
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    my %flt;
    for (my $i=1;$i<@a;$i++){
        $a[$i]=~/^(\w+)\|(\S+)$/;
        $flt{$1}{$2}++;
    }
    my @k1=sort keys %flt;
    if (scalar(@k1) == 1){
        my @k2=keys %{$flt{$k1[0]}};
        my $num=scalar(@k2);
        $max{$k1[0]}{$num}++;
        $h{$k1[0]}{GF}{uniq}++;
        $h{$k1[0]}{GF}{GF}++;
        #print "\$h{$k1[0]}{GF}{gene} += ",scalar(@k2),"\n";
        $h{$k1[0]}{GF}{gene} += scalar(@k2);
    }else{
        for my $k1 (@k1){
            my @k2=keys %{$flt{$k1}};
            my $num=scalar(@k2);
            $max{$k1}{$num}++;
            $h{$k1}{GF}{GF}++;
            $h{$k1}{GF}{gene} += scalar(@k2);
        }
    }
}
close F;

print "Species\tTotal_genes\tGenes_in_families\tUnclustered_genes\tFamilies\tUnique_families\tGenes_per_family\tMaximum_gene_family_size\n";
for my $k (sort keys %h){
    my $allgenenum=$h{$k}{gene};
    my $allGF=$h{$k}{GF}{GF};
    my $allGFgene=$h{$k}{GF}{gene};
    my $uniqGF=$h{$k}{GF}{uniq};
    my $unclusterd=$allgenenum-$allGFgene;
    my $geneperGF=$allGFgene/$allGF;
    my @max=sort{$a<=>$b} keys %{$max{$k}};
    print "$k\t$allgenenum\t$allGFgene\t$unclusterd\t$allGF\t$uniqGF\t$geneperGF\t$max[-1]\n";
}
