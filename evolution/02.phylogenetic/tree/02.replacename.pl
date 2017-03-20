#!/usr/bin/perl
use strict;
use warnings;

my ($nameinfo,$tree)=@ARGV;
die "give the nameinfo and tree\n" if (! $nameinfo);
my %h;
open (F,"$nameinfo");
while (<F>) {
    chomp;
    /^(\S+)\s+(\S+)$/;
    $h{$1}=$2;
}
close F;

open (F,"$tree");
while (<F>) {
    chomp;
    for my $k (sort keys %h){
        #s/$k\:/$h{$k}\:/;
        s/$k/$h{$k}/;
    }
    print "$_\n";
}
