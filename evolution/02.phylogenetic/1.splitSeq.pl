use strict;
use warnings;
use Bio::SeqIO;

## created by Yongzhi Yang. 2017/3/20 ##

my ($dir,$singlecoypegenelist)=@ARGV;
die "Usage:\nperl $0 dir single_coype_gene_list
dir: containing cds and pep seq with the name of XXX.cds.clean.fa and XXX.pep.clean.fa. XXX represent the species ID.
single_coype_gene_list: clusterid\tGene_num\tSpecies_num\tOrthologous_genes_split_comma
" if (! $singlecoypegenelist);

my %list;
my %seq;
my @seq=<$dir/*clean.fa>;
for my $seqfile (@seq){
    $seqfile=~/\/(\w+)\.(\w+)\.clean.fa$/;
    my $type=$2;
    my $name=$1;
    my $fa=Bio::SeqIO->new(-file=>"$seqfile",-format=>'fasta');
    while (my $seq=$fa->next_seq) {
        my $id=$seq->id;
        my $seq=$seq->seq;
        my $key="$name|$id";
        $seq{$type}{$key}=$seq;
    }
}

open (F,"$singlecoypegenelist")||die"$!";
while (<F>) {
    chomp;
    next if /^#/;
    my @a=split(/\t/,$_);
    my @b=split(/,/,$a[3]);
    for (my $i=0;$i<@b;$i++){
        $list{$a[0]}{$b[$i]}++;
    }
}
close F;

`mkdir align` if (! -e "align");
for my $k (sort keys %list){
    my $dir="align/$k";
    `mkdir $dir` if (! -e "$dir");
    for my $type ("cds","pep"){
        open (O,">$dir/$type");
        for my $k2 (sort keys %{$list{$k}}){
            die "$k2\n" if ! exists $seq{$type}{$k2};
            print O ">$k2\n$seq{$type}{$k2}\n";
        }
        close O;
    }
}
