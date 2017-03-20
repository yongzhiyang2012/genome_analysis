#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $in=shift or die "perl $0 input_fasta\n";

open (O,">$in.name.txt");

my $i='aaaa';
my $fa=Bio::SeqIO->new(-format=>"fasta",-file=>"$in");
while (my $seq=$fa->next_seq) {
    my $id=$seq->id;
    my $seq=$seq->seq;
    print ">N$i","M\n$seq\n";
    print O "N$i","M\t$id\n";
    $i++;
}
