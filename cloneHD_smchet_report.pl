#!/usr/bin/perl -w
use strict;
use Getopt::Std;

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

$,= " ";

my %OPTS;
getopts('i:o:',\%OPTS);

my $resFile = $OPTS{i};

my $out = $OPTS{o};
my $l=0;

my %frac=();
my %mutC=();
my @mut=();
my @probV=();

open my $IN_FILE, '<', $resFile or die "Cannot open input file";
foreach my $line (<$IN_FILE>) {
	chomp($line);
	
	$l++;
	
	my @fields=split(/\s+/,$line);
	
	if($l==1){
		my $c=1;
		for(my $i=2; $i<scalar(@fields); $i+=2){
			$frac{$c}=$fields[$i];
	    $mutC{$c}=0;
	    $c++;
		}
	}
	
	if(!($line =~ /\#/)){ 
		my $chr=$fields[0];
		my $coord=$fields[1];
		
		my @prob=();
		for(my $i=2; $i<scalar(@fields); $i++){
			push(@prob,$fields[$i]);
		}
		push(@probV,$line);
		
		my $clusterID=randomDraw(\@prob);
		$mutC{$clusterID}++;
		
		push(@mut,$clusterID);
	};
};
close $IN_FILE;

my $nClusters=scalar(keys(%frac));
my $fpCluster=0;

for(my $i=1; $i<=$nClusters; $i++){
	if($frac{$i}<0.00001){
		$fpCluster++;
	}
}

print $nClusters-$fpCluster,"\n";

my $purity=$frac{1};

$,= "\t";

my $fileOutPurity=$out.".1A.txt";
my $fileOutNc=$out.".1B.txt";
my $fileOutClusters=$out.".1C.txt";
my $fileOutAssignment=$out.".2A.txt";

open(my $fh, '>', $fileOutPurity) or die "Could not open file '$fileOutPurity' $!";
print $fh $purity,"\n";
close($fh);      
 
open($fh, '>', $fileOutNc) or die "Could not open file '$fileOutNc' $!";
print $fh $nClusters-$fpCluster,"\n";
close($fh);     

open($fh, '>', $fileOutClusters) or die "Could not open file '$fileOutClusters' $!";
for(my $i=1; $i<=$nClusters; $i++){
	print $fh $i,$mutC{$i},$frac{$i},"\n";
}
close($fh);

open($fh, '>', $fileOutAssignment) or die "Could not open file '$fileOutAssignment' $!";
for(my $i=0; $i<scalar(@mut); $i++){
	print $fh $mut[$i],"\n";
}
close($fh);

sub randomDraw{
	
	my ($refProb) = @_;
	
	my @prob=@$refProb;
	my @cprob=();
	
	push(@cprob,$prob[0]);
	for(my $i=1; $i<scalar(@prob); $i++){
		push(@cprob,$prob[$i]+$cprob[$i-1]);
	}
	
	for(my $i=0; $i<scalar(@cprob); $i++){
		my $rn=rand(1.0);
		if($rn<$cprob[$i]){
			return ($i+1);
		}
	}
	return 1;
}

sub CoClusterP{
	my ($line1, $line2) = @_;
	
	my @prob1=split(/\s+/,$line1);
	my @prob2=split(/\s+/,$line2);
	my $sum=0.0;
	for(my $i=2; $i<scalar(@prob1); $i++){
		$sum+=$prob1[$i]*$prob2[$i];
	}
	return($sum);
}