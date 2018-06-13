use strict;
use warnings;

#my $depth_list="00.Ore.all.depth.txt";
my $depth_list="00.all.depth.txt";
my $outdir='psmc';
#my $ref="/home/pool/users/yangyongzhi2012/Ostrya_resequnce/00.ref/Orehderiana.genome.2K.fa";
my $ref="/home/pool/users/yangyongzhi2012/Ostrya_resequnce/00.ref/Ore.final.assembly.flt.2k.fa";
my $samtools="samtools";
my $bcftools="bcftools";
my $vcfutils_pl="vcfutils.pl";
my $bamdir="/home/pool/users/yangyongzhi2012/Ostrya_resequnce/02.map-LZU-withCP/3.realign";
`mkdir $outdir` if (! -e "$outdir");
###Read depth file
my %dp;
open (F,"$depth_list")||die"$!";
while (<F>) {
    chomp;
    next if /^(#|ID)/;
    my @a=split(/\s+/,$_);
    next if $a[5]<20;
    $dp{$a[0]}=$a[5];
}
close F;

my @file=<$bamdir/*realn.bam>;
open(OUT,"> $0.01.sh");
open(OUT2,"> $0.02.bootstrap.sh");
foreach my $f(@file){
    chomp $f;
    $f =~ /\/([^\/]+)\.realn\.bam/;
    my $id=$1;
    next if $id ne 'Omu01';
    next if ! exists $dp{$id};
    next if $id=~/Ore/;
    my $min=int($dp{$id}/3)-1;
    my $max=int($dp{$id}*2)+1;
    print OUT "$samtools mpileup -C50 -uf $ref $f | $bcftools view -c -  | $vcfutils_pl vcf2fq -d $min -D $max -Q 10 -l 5 | gzip > $outdir/$id.fa.gz; ";
    print OUT "fq2psmcfa -q20 $outdir/$id.fa.gz > $outdir/$id.psmcfa; ";
    print OUT "psmc -N25 -t15 -r5 -p \"4+25*2+4+6\" -o $outdir/$id.psmc $outdir/$id.psmcfa\n";
    print OUT2 "splitfa $outdir/$id.psmcfa > $outdir/$id.split.psmcfa; mkdir -p $outdir/$id.bt; ";
    print OUT2 "seq 100 | xargs -P 20 -i psmc -N25 -t15 -r5 -b -p \"4+25*2+4+6\" -o $outdir/$id.bt/$id.{}.psmc $outdir/$id.split.psmcfa | sh; ";
    print OUT2 "cat $outdir/$id.psmc $outdir/$id.bt/$id.*.psmc > $outdir/$id.bt/$id-bootstrap.psmc\n";

}
close OUT;
close OUT2;
