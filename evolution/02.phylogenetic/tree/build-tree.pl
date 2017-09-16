#!/usr/bin/perl
use strict;
use warnings;
use Cwd 'abs_path';

## build tree command creating scripts ##
## if you sequence id length larger than 10, please use 01.fixname.pl to rename you seqid in the fastafile and use 02.replacename.pl to retrevie seqid in the treefile ##
## required software ##
my $clustalw="/home/share/software/ClustalW/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2";
my $raxml="/home/share/software/RAxML/standard-RAxML-master/raxmlHPC-PTHREADS";
my $phyml="/home/share/software/phyml/PhyML-3.1/PhyML-3.1_linux64";
my $mrbayes="/home/share/software/mrbayes/mrbayes-3.2.4/src/mb";
## done ##

my ($input_alignment,$software,$bootstrap)=@ARGV;
die "Usage: perl $0 input_alignment <mrbayes|phyml|raxml|all> bootstrap\nThis scripts will creat the running dir and the running command\nNotice: for mrbayes bootstrap is also necessary, but don't read\n" if (! $bootstrap);

$input_alignment=&check_input_alignment;
&check_software;

sub check_input_alignment{
    my $abs_input_alignment=abs_path($input_alignment);
    return ($abs_input_alignment);
}
sub check_software{
    if ($software =~ /mrbayes|phyml|raxml|all/){
        for my $sf ("mrbayes","phyml","raxml"){
            next if (($software ne $sf) && ($software ne 'all'));
            `mkdir $sf` if (! -e "$sf");
            chdir($sf);
            `ln -s $input_alignment .`;
            chdir("..");
            $input_alignment=~/\/([^\/]+)$/;
            my $base_alignment_name=$1;
            my ($print_convertType_alignment,$out_base_alignment_name)=&print_convertType_alignment($base_alignment_name,$sf);
            my $print_software_command=&print_software_command($out_base_alignment_name,$sf);
            print "cd $sf ; $print_convertType_alignment ; $print_software_command ; cd ..\n";
        }
    }else{
        die "Please give the right software\n";
    }
}
sub print_convertType_alignment{
    my %sftype;
    $sftype{mrbayes}="NEXUS";
    $sftype{phyml}="PHYLIP";
    $sftype{raxml}="PHYLIP";
    my ($lninput,$sf)=@_;
    return ("$clustalw -INFILE=$lninput -CONVERT -TYPE=DNA -OUTFILE=$lninput.$sftype{$sf}  -OUTPUT=$sftype{$sf}","$lninput.$sftype{$sf}");
}
  sub print_software_command{
      my ($in,$sf)=@_;
      if ($sf eq 'mrbayes'){
          return ("mpirun -np 24 $mrbayes ; execute $in ; lset nst=6 rates=invgamma ; mcmc nchains=24 ngen=1000000 samplefreq=500 printfreq=500 diagnfreq=500 ; no ; sump relburnin=yes ; sumt relburnin=yes ; exit");
      }elsif($sf eq 'phyml'){
          return ("$phyml -i $in -d nt -b $bootstrap -m GTR -f e -a e --no_memory_check");
      }elsif($sf eq 'raxml'){
          return ("$raxml -s $in -n $in -m GTRGAMMAI -f a -x 12345 -N 100 -p 12345 -T 30");
      }
  }
