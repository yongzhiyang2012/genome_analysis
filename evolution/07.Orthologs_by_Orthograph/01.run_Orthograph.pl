use strict;
use warnings;

my %par=&read_par("00.data.list");
`mkdir needed_parameters` if (! -e "needed_parameters");
`mkdir input_data` if (! -e "input_data");
open (SH,">$0.sh");
print SH "orthograph-manager -c orthograph.conf --create < needed_parameters/00.create.par\n";
open (P,">needed_parameters/00.create.par")||die"$!";
print P "y\n";
close P;

my @sp=sort keys %{$par{pep}};
for my $sp (@sp){
    my $inputsppep=$par{pep}{$sp};
    open (P,">needed_parameters/01.load.pep.$sp.par");
    print P "$sp\nv1.0\n";
    close P;
    print SH "orthograph-manager -c orthograph.conf  --load-ogs-peptide $inputsppep < needed_parameters/01.load.pep.$sp.par\n";
}

open (F,$par{Orthomcl_output})||die"no such file: $par{Orthomcl_output}\n";
open (O,">input_data/Orthomcl_output.tables");
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    $a[0]=~s/\://;
    for (my $i=1;$i<@a;$i++){
        $a[$i]=~/^(\w+)\|/ or die "wrong Orthomcl_output format\n";
        print O "$a[0]\t$a[$i]\t$1\n";
    }
}
close F;
close O;
open (P,">needed_parameters/02.load.ortho.table.par")||die"$!";
print P "OGStable\nNA\n";
for (my $i=1;$i<=@sp;$i++){
    print P "$i\n";
}
close P;
print SH "orthograph-manager -c orthograph.conf input_data/Orthomcl_output.tables < needed_parameters/02.load.ortho.table.par\n";
print SH "orthograph-analyzer -c orthograph.conf --prepare\northograph-analyzer -c orthograph.conf\northograph-reporter -c orthograph.conf\n";
close SH;

my @k=sort keys %{$par{classify_protein}};
@k=("input_data/need_classify_protein.fa",@k);
&sub_write_protein(@k);

open (CONFIG,">orthograph.conf")||die"$!";
print CONFIG "alignment-program    = mafft-linsi
blast-program        = blastp
database-backend   = sqlite
exonerate-program    = exonerate
hmmbuild-program     = hmmbuild
hmmsearch-program    = hmmsearch
input-file         = input_data/need_classify_protein.fa
makeblastdb-program  = makeblastdb
num-threads                = 30
ortholog-set       = OGStable
output-directory   = Orthograph_output/output
species-name       = all_SP
sqlite-database    = Orthograph_output/orthodb.sqlite
sqlite-program     = /usr/bin/sqlite3
translate-program    = fastatranslate
";
close CONFIG;

sub sub_write_protein{
    my @inputkey=@_;
    my $output_file=shift @inputkey;
    open (O,">$output_file")||die"$!";
    for my $infa (@inputkey){
        use Bio::SeqIO;
        my $build_fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$infa");
        while (my $seq=$build_fa->next_seq) {
            my $id=$seq->id;
            my $seq=$seq->seq;
            print O ">$id\n$seq\n";
        }
    }
    close O;
}
sub read_par{
    my %r;
    my ($inpar)=@_;
    open (F,"$inpar")||die"$!";
    while (<F>) {
        chomp;
        my @a=split(/\s+=\s+/,$_);
        if ($a[0]=~/pep:(\S+)/){
            $r{pep}{$1}=$a[1];
        }elsif($a[0] eq 'classify_protein'){
            $r{$a[0]}{$a[1]}++;
        }else {
            $r{$a[0]}=$a[1];
        }
    }
    close F;
    return %r;
}
