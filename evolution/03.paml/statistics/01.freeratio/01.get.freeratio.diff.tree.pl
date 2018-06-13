use strict;
use warnings;

my $dir=shift or die "perl $0 dir\nthe dir should containing paml result with the following format: dir/clusterid/signal.mlc\n";

my %tiplist=&read_tiplist();

open (O,">All.freeratio.result.all_tip.out");
print O "cluster\tlength\tspeices\tN\tS\tdN\tdS\tNdN\tSdS\tw\tlnl\n";

my @mlc=<$dir/*/branch.freeratio.mlc>; 
for my $mlc (@mlc){
    $mlc=~/\/([^\/]+)\/branch.freeratio.mlc$/;
    my ($cluster)=($1);
    my %tipfeature=&read_tree_feature("$dir/$cluster");
    my ($lnl,$N,$S,$len,$outsp,$dstree,$dntree,$wtree);

    open (F,"$mlc");
    my $line=0;
    while (<F>){
        chomp;
        $line++;
        if ($line == 1){
            /^\s+\d+\s+(\d+)$/;
            $len = $1;
        }
        if (/^lnL\(ntime:\s*\d+\s+np:\s*\d+\):\s*(\S+)/){
	    $lnl=$1;
        }
        if (/^\s+branch\s+t\s+N\s+S\s+/){
            while (<F>){
	chomp;
	next if /^$/;
	last if /^tree\s+/;
	s/^\s+//;
	my @a=split(/\s+/,$_);
	$N=$a[2];
	$S=$a[3];
            }
        }
        if (/^dS\s+tree:/){
            $dstree=<F>;
            chomp $dstree;
        }
        if (/^dN\s+tree:/){
            $dntree=<F>;
            chomp $dntree;
        }
        if (/^w/){
            $wtree=<F>;
            chomp $wtree;
        }
    }
    close F;

    my %out;
    for my $key (sort keys %tipfeature){
        my $v=$tipfeature{$key};
        my $time = $v =~ tr/\)/\)/;
        $v=~s/\)+//;
        next if $dstree !~ /\Q$v\E\W+/;
        if ($time == 0){
            $dstree=~/$v\:\s*([0-9\.]+)/ or die "ds0:\t$dstree\n$key\n$v\n";
            $out{$key}{ds}=$1;
            $dntree=~/$v\:\s*([0-9\.]+)/ or die "dn0:\t$dntree\n";
            $out{$key}{dn}=$1;
            $wtree=~/$v\s*\#([0-9\.]+)/ or die "w0:\t$wtree\n";
            $out{$key}{w}=$1;
        }else{
            $dstree=~/$v(\:\s*[0-9\.]+\s*\)){$time}\:\s*([0-9\.]+)/ or die "ds1:\t$dstree\n$key\n$v\n";
            $out{$key}{ds}=$2;
            $dntree=~/$v(\:\s+[0-9\.]+\s*\)){$time}\:\s*([0-9\.]+)/ or die "dn1:\t$dntree\n";
            $out{$key}{dn}=$2;
            $wtree=~/$v(\s*\#[0-9\.]+\s*\)){$time}\s*#([0-9\.]+)/ or die "w1:\t$wtree\n";
            $out{$key}{w}=$2;
        }
    }
    for my $k1 (sort keys %out){
	print  O "$cluster\t$len\t$k1\t$N\t$S\t$out{$k1}{dn}\t$out{$k1}{ds}\t",$N*$out{$k1}{dn},"\t",$S*$out{$k1}{ds},"\t$out{$k1}{w}\t$lnl\n";
    }
}
close O;

sub read_tree_feature{
    my ($in)=@_;
    my %r;
    for my $k (sort keys %tiplist){
        my $tree="$in/tree/tree.$k";
        if (-e "$tree"){
            open (F,"$tree");
            my $treeline=<F>;
            chomp $treeline;
            close F;
            if ($treeline=~/\#1/){
		$treeline=~/(\w+\)*)\s+\#1/; #  "$in\n$treeline\n$k\n";
		$r{$k}=$1;
	    }else{
		print "$in\t$k\n";
		next;
	    }
        }else{
            $r{$k}=$k;
        }
    }
    return %r;
}
sub read_tiplist{
    my %r;
    open (F,"need.tip.list")||die"no file need.tip.list\n";
    while (<F>) {
        chomp;
        $_=~/^(\S+)\s+/;
        $r{$1}++;
    }
    close F;
    return %r;
}
