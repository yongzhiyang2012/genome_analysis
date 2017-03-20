use strict;
use warnings;

## created by Yongzhi Yang. 2017/3/20 ##

my $inputfile=shift or die "Please give the input information file.\nThe example file is inputseq.txt\n";

`mkdir compliantFasta` if (! -e "compliantFasta");

open (O,">$0.sh");
open (F,"$inputfile");
while (<F>) {
    chomp;
    /^(\S+)\s+(\S+)/;
    print O "cd compliantFasta ; orthomclAdjustFasta $1 $2 1 ; cd ..\n";
}
close F;
close O;
