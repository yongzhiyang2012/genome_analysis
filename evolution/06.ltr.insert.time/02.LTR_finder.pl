use strict;
use warnings;

my $ltrfinder_software=shift or die "perl $0 <ltr_finder_abs_path>\n";
$ltrfinder_software=~/^(\S+)\/ltr_finder$/ or die "give the abs ltr_finder path and end with ltr_finder\n";
my $trna_data="$1/tRNAdb/Athal-tRNAs.fa";

my $ltrfinder_pl="used_softwares/dawgpaws/scripts/batch_ltrfinder.pl";
my @in=`ls splitBYscaff`;
my $outdir="LTR_finder.batch.result";
`mkdir $outdir` if (! -e "$outdir");

open (O,">$0.sh");
for my $in (@in){
    chomp $in;
    my $indir="splitBYscaff/$in";
    print O "$ltrfinder_pl -i $indir -o $outdir -c used_softwares/dawgpaws/scripts/config/batch_ltrfinder.jcfg --ltr-finder $ltrfinder_software -s $trna_data -g -f\n";
}
close O;
