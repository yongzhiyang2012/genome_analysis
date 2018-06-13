use strict;
use warnings;
use Bio::SeqIO;

my $genome=shift or die "perl $0 input_genome\n";
my $dir2="splitBYscaff";
`mkdir $dir2` if (! -e "$dir2");
my $fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$genome");
while (my $seq=$fa->next_seq){
    my $id=$seq->id;
    my $seq=$seq->seq;
    `mkdir $dir2/$id` if (! -e "$dir2/$id");
    open (O,">$dir2/$id/$id.fa");
    print O ">$id\n$seq\n";
    close O;
}
