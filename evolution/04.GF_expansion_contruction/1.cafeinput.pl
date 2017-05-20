use strict;
use warnings;

### phase gene family size >= (species_num/2) && >=3 species in this family ###
my $in=shift or die "perl $0 orthomclOUTPUT\n";
my $out="cafeinput.tab";
my %h;
my %id;
my %flt;
open (F,"$in")||die"$!";
while (<F>) {
    chomp;
    my @a=split(/\s+/,$_);
    $a[0]=~s/\://;
    for (my $i=1;$i<@a;$i++){
        $a[$i]=~/^(\w\w\w)\|/ or die "wrong type: $a[0]\t$a[$i]\n";
        my $species=$1;
        $flt{$a[0]}{sp}{$species}++;
        $flt{$a[0]}{num}{$a[$i]}++;
        $h{$a[0]}{$species}++;
        $id{$species}++;
    }
}
close F;

open (O,">$out");
my @id=sort keys %id;
print O "FAMILYDESC\tFAMILY\t",join("\t",@id),"\n";
for my $k (sort keys %h){
    next if scalar(keys %{$flt{$k}{sp}}) < 3;
    next if scalar(keys %{$flt{$k}{num}}) < int(scalar(@id)/2);
    print O "UNKNOWN\t$k";
    for my $k2 (@id){
        $h{$k}{$k2}=0 if ! exists $h{$k}{$k2};
        print O "\t$h{$k}{$k2}";
    }
    print O "\n";
}
close O;
