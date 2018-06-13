use strict;
use warnings;
use Bio::SeqIO;

my $mucle="/home/share/users/yangyongzhi2012//anaconda2/bin/muscle"; # path-to-mucle;
my $clustalw2="/home/share/software/ClustalW/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2"; # path-to-clustalw2;
my $dnadist="/home/share/software/phylip/phylip-3.696/exe/dnadist"; # path-to-dnadist;
die "please check the path of the software used in this script\n" if ((! -e "$mucle") || (! -e "$clustalw2") || (! -e "$dnadist"));

my %h;
my @inclassified=<genome.ltrfinder.seq.fa.split/*/*classified>;
for my $inclassified (@inclassified){
    my $fa=Bio::SeqIO->new(-file=>"$inclassified",-format=>"fasta");
    while (my $seq=$fa->next_seq) {
        my $id=$seq->id;
        next unless $id=~/\_LTRFinder\_/;
        $id=~/(\S+)\#(\S+)/;
        my ($id_name,$id_family)=($1,$2);
        $id_family="LTR_unknown" unless $id_family=~/^LTR/;
        $id_family=~s/\//_/g;
        $h{$id_name}=$id_family;
    }
}
my $dir="align_soloLTRs";
`mkdir $dir` if (! -e "$dir");

my %seq;
&readfasta("genome.ltrfinder.three_prime_LTR.fa");
&readfasta("genome.ltrfinder.five_prime_LTR.fa");

open (SH,">$0.sh");
for my $id_name (sort keys %h){
    my $dir2="$dir/$id_name";
    `mkdir $dir2` if (! -e "$dir2");
    my $id_family=$h{$id_name};
    open (O,">$dir2/$id_name.3.5.ltr.fa");
    print O ">$id_name#$id_family#five_prime_LTR\n$seq{$id_name}{five_prime_LTR}\n";
    print O ">$id_name#$id_family#three_prime_LTR\n$seq{$id_name}{three_prime_LTR}\n";
    close O;
    open (O,">$dir2/$id_name.3.5.ltr.fa.align.phylip");
    print O "$id_name.3.5.ltr.fa.align.phy\nD\nY\n";
    close O;
    print SH "cd $dir2; $mucle -in $id_name.3.5.ltr.fa -out $id_name.3.5.ltr.fa.align ; $clustalw2 -INFILE=$id_name.3.5.ltr.fa.align  -CONVERT -TYPE=DNA -OUTFILE=$id_name.3.5.ltr.fa.align.phy -OUTPUT=PHYLIP ; $dnadist < $id_name.3.5.ltr.fa.align.phylip ; mv outfile $id_name.3.5.ltr.fa.align.dnadist; cd ../../\n";
}
close SH;

sub readfasta{
    my ($infile_fa)=@_;
    die "no $infile_fa\n" if (! -e "$infile_fa");
    $infile_fa=~/(five_prime_LTR|three_prime_LTR)/;
    my $type=$1;
    my $infile_fa_read=Bio::SeqIO->new(-file=>"$infile_fa",-format=>"fasta");
    while (my $seq=$infile_fa_read->next_seq){
        my $id=$seq->id;
        my $seq=$seq->seq;
        $seq{$id}{$type}=$seq;
    }
}
