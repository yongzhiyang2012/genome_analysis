use strict;
use warnings;
use Bio::SeqIO;

## created by Yongzhi Yang. 2017/3/20 ##

my %fourDegenerateSite=(
	        'TC'=>'Ser',
	        'CT'=>'Leu',
	        'CC'=>'Pho',
	        'CG'=>'Arg',
	        'AC'=>'Thr',
	        'GT'=>'Val',
	        'GC'=>'Phe',
	        'GG'=>'Gln',
	       );

my %seq;
my %count;
#my $in=shift or die "perl $0 \$inFasta output\n";
#my $out=shift or die "perl $0 \$inFasta output\n";
my $in="3.connect.pl.connect.cds.fa";
my $out="$0.connect4Dsites.fa";
my $fa=Bio::SeqIO->new(-file=>"$in",-format=>"fasta");
while (my $seq=$fa->next_seq) {
    my $id=$seq->id;
    my $seq=$seq->seq;
    my @seq=split(//,$seq);
    for (my $i=0;$i<@seq;$i=$i+3){
        my $key=$seq[$i].$seq[$i+1];
        my $value=$seq[$i+2];
        if (exists $fourDegenerateSite{$key}){
            $seq{$id}{$i}=$value;
            print "$id\t$i\t$value\t",$seq[$i-2],$seq[$i-1],$seq[$i],$seq[$i+1],$seq[$i+2],"\n" if ($value eq '-');
            $count{$i}++;
        }
    }
}


open (O,">$out");
open (O1,">$out.list");
my @k1=sort keys %seq;
for my $k1 (@k1){
    print O ">$k1\n";
    print O1 "$k1\n";
    for my $k2 (sort{$a<=>$b} keys %{$seq{$k1}}){
        my $count=$count{$k2};
        if ($count==scalar(@k1)){
            print O "$seq{$k1}{$k2}";
            print O1 $k2+2," ";
            #print "$k1\t$k2\t$seq{$k1}{$k2}\n" if $seq{$k1}{$k2} eq '-';
        }
    }
    print O "\n";
    print O1 "\n";
}
close O;
close O1;
