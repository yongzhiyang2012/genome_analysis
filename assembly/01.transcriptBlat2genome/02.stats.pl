use strict;
use warnings;

my %blocklen;
my %count;
my %all;
my %select;
my $in=shift or die "perl $0 <blat out file> <trancripts info>\n";
my $inlen=shift or die "perl $0 <blat out file> <trancripts info>\n";

my $out="$in.stats";
$in=~/[\/]*([^\/]+)\.out$/ or die "$in\n";
my $type=$1;
open (F,"$inlen")||die"$!";
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    $select{$a[0]}++;
    if ($a[1]>200){
        $all{200}{num}++;
        $all{200}{len} += $a[1];
    }
    if ($a[1]>500){
        $all{500}{num}++;
        $all{500}{len} += $a[1];
    }
    if ($a[1]>1000){
        $all{1000}{num}++;
        $all{1000}{len} += $a[1];
    }
}
close F;

open (F,"$in")||die"$!";
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    next unless /^\d+\s+\d+/;
    next unless exists $select{$a[9]};
    my $qlen=$a[10];
    my @b=split(/,/,$a[18]);
    my @c=split(/,/,$a[19]);
    $blocklen{$a[9]}{qlen}=$a[10];
    for (my $i=0;$i<scalar(@b);$i++){
        for (my $j=$c[$i];$j<$c[$i]+$b[$i];$j++){
            $blocklen{$a[9]}{blocklen}{$a[13]}{$j}++;
        }
    }
}
close F;
for my $k1 (sort keys %blocklen){
    my $qlen=$blocklen{$k1}{qlen};
    my $blocklen=0;
    for my $k2 (sort keys %{$blocklen{$k1}{blocklen}}){
        my @k=keys %{$blocklen{$k1}{blocklen}{$k2}};
        $blocklen=scalar(@k) if scalar(@k)>$blocklen;
    }
    my $per=$blocklen/$qlen;
    for my $lencutoff (200,500,1000){
        if ($qlen>$lencutoff){
            $count{$lencutoff}{all}++;
            $count{$lencutoff}{50}++ if $per>0.5;
            $count{$lencutoff}{90}++ if $per>0.9;
        }
    }
}
open (O,">$out");
print O "Data_set\tLength_type\tNumber\tTotal_length\tCovered_by_assembly_(\%)\t\tWith_>_90\%\_sequence_in_one_scaffold\t\tWith_>_50\%\_sequence_in_one_scaffold\n";
print O "\t\t\t\t\tNumber\tPercentage(\%)\tNumber\tPercentage(\%)\n";
for my $k (sort{$a<=>$b} keys %count){
    my $allqnum=$all{$k}{num};
    my $allqlen=$all{$k}{len};
    my $targetqall=$count{$k}{all};
    my $targetq50=$count{$k}{50};
    my $targetq90=$count{$k}{90};
    print O "$type\t>$k","bp\t$allqnum\t$allqlen\t",$targetqall/$allqnum,"\t$targetq90\t",$targetq90/$allqnum,"\t$targetq50\t",$targetq50/$allqnum,"\n";
}
close O;
