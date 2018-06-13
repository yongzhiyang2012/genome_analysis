use strict;
use warnings;
use Bio::SeqIO;

my $genome_fa=shift or die "perl $0 <genome_fa>\n";

my %gff;
my @gff=<LTR_finder.batch.result/*/gff/*_def_ltr_finder.gff>;
for my $gff (@gff){
    open (F,"$gff")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\t/,$_);
        if ($a[2] =~ /five_prime_LTR|three_prime_LTR|LTR_retrotransposon/){
            $a[8]=~/ltr_finder_(\d+)/;
            my $key="$a[0]"."_LTRFinder_def$1";
            $gff{$a[0]}{$key}{$a[2]}{start}=$a[3];
            $gff{$a[0]}{$key}{$a[2]}{end}=$a[4];
            $gff{$a[0]}{$key}{$a[2]}{strand}=$a[6];
        }
    }
    close F;
}

open (O,">genome.ltrfinder.seq.fa");
open (O1,">genome.ltrfinder.five_prime_LTR.fa");
open (O2,">genome.ltrfinder.three_prime_LTR.fa");

my $fa=Bio::SeqIO->new(-file=>"$genome_fa",-format=>"fasta");
while (my $seq=$fa->next_seq) {
    my $id=$seq->id;
    my $seq=$seq->seq;
    if (exists $gff{$id}){
        for my $k (sort keys %{$gff{$id}}){

            my ($start0,$end0,$strand0)=($gff{$id}{$k}{LTR_retrotransposon}{start},$gff{$id}{$k}{LTR_retrotransposon}{end},$gff{$id}{$k}{LTR_retrotransposon}{strand});
            die "$k\n" if $end0 < $start0;
            my $len0=$end0 - $start0 + 1;
            $start0=$start0 - 1;
            my $newseq0=substr($seq,$start0,$len0);
            if ($strand0 eq '-'){
	$newseq0=reverse $newseq0;
	$newseq0=~tr/[ATCG]/[TAGC]/;
            }
            print O ">$k\n$newseq0\n";

            my ($start1,$end1,$strand1)=($gff{$id}{$k}{five_prime_LTR}{start},$gff{$id}{$k}{five_prime_LTR}{end},$gff{$id}{$k}{five_prime_LTR}{strand});
            die "$k\n" if $end1 < $start1;
            my $len1=$end1 - $start1 + 1;
            $start1=$start1 - 1;
            my $newseq1=substr($seq,$start1,$len1);
            if ($strand1 eq '-'){
	$newseq1=reverse $newseq1;
	$newseq1=~tr/[ATCG]/[TAGC]/;
            }
            print O1 ">$k\n$newseq1\n";

            my ($start2,$end2,$strand2)=($gff{$id}{$k}{three_prime_LTR}{start},$gff{$id}{$k}{three_prime_LTR}{end},$gff{$id}{$k}{three_prime_LTR}{strand});
            die "$k\n" if $end2 < $start2;
            my $len2=$end2 - $start2 + 1;
            $start2=$start2 - 1;
            my $newseq2=substr($seq,$start2,$len2);
            if ($strand2 eq '-'){
	$newseq2=reverse $newseq2;
	$newseq2=~tr/[ATCG]/[TAGC]/;
            }
            print O2 ">$k\n$newseq2\n";
        }
    }
}
close O;
close O1;
close O2;

