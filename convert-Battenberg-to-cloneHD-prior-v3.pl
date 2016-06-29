#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use POSIX;
#use List::Util qw(first max maxstr min minstr reduce shuffle sum);

$,= " ";

my %OPTS;
getopts('c:g:m:o:',\%OPTS);

my $copynumberFile = $OPTS{c};
#my $purity=1.0;
my $gender=$OPTS{g};
my $mode=$OPTS{m};
my $fileOut=$OPTS{o};

my $flag=0;
my @states;
my $nMax=8;

#my %cnSegments();
open(my $fh, '>', $fileOut) or die "Could not open file '$fileOut' $!";
if($mode eq "mean-tcn"){
                        print $fh "#chr first nLoci last meanCN\n";
}
else{
	print $fh "#chr first nLoci last  genotype-availability via baf\n";
}



my @totalCNV;
my @logRV;
my @infoV;


open my $IN_FILE, '<', $copynumberFile or die "Cannot open  input file";
foreach my $line (<$IN_FILE>) {
    chomp($line);
   
   my @fields=split(/\t+/,$line);
	
	if($fields[0] ne  "chr"){
		my $chr=$fields[0];
		if($chr eq "X"){
			$chr=23;
		}
		my $first=$fields[1];
		my $nLoci=1;
		my $last=$fields[2];
		
		my $major1=$fields[7];
		my $minor1=$fields[8];
		my $frac1=$fields[9];
		my $total1=$major1+$minor1;
		my $logR=$fields[5];
		my $major2=$fields[10];
		my $minor2=$fields[11];
		my $frac2=$fields[12];
		my $total2="NA";
		
		if($major2 ne "NA"){
			$total2=$major2+$minor2;
		}
	
		my $normal=2;
		if($gender eq "male" && $chr == 23){
			$normal=1;
		}

		my $total=$frac1*$total1;	
		my $major=$major1;
		my $totalR=$normal*2.0**$logR;
	
		if($major2 ne "NA"){
			$total=$frac1*$total1+$frac2*$total2;
			$major=ceil($frac1*$major1+$frac2*$major2);
		}
		
		
		#my $meanCN=$purity*$total+(1-$purity)*$normal;
		#print $major1, $minor1, $total1,$major2, $minor2, $total2,"\n";
		my $info=$chr." ".$first." ".$nLoci." ".$last;
		push(@totalCNV,$total);
		push(@logRV,$logR);
		push(@infoV,$info);
		
		#if($mode eq "mean-tcn"){
			#print $chr, $first, $nLoci, $last,$meanCN,"\n"
		#	print $chr, $first, $nLoci, $last,$totalR,"\n";
		#	}
		#}
		if($mode eq "avail-cn"){
		   printf $fh "%d %d %d %d ",$chr, $first, $nLoci, $last;
		   for(my $i=0; $i<=$nMax; $i++){
			if($i<= $major){
				print $fh "1.0 ";
			}
			else{
				print $fh "0.0 ";
			}
		   }			
			print $fh "\n";
		}

	}
 }
close $IN_FILE;

#close($fh);  
#open(my $fh, '>', $fileOut) or die "Could not open file '$fileOut' $!";

if($mode eq "mean-tcn"){

# now select the state that has cancer only copynumber closest to normal.
my $gauge=0.0;
my $min=1000.0;

for(my $i=0; $i<scalar(@infoV); $i++){
		
		if(abs($totalCNV[$i]-2.0)<$min){
			$gauge=$logRV[$i];
			$min=abs($totalCNV[$i]-2.0);
		} 
}

# Use that gauge to fix the d.o.f for the copynumber data 
for(my $i=0; $i<scalar(@infoV); $i++){
	
	my $logR=$logRV[$i];
	
	my @fields=split(/\s+/,$infoV[$i]);
	
	my $normal=2;
		if($gender eq "male" && $fields[0] == 23){
			$normal=1;
	}
	my $totalR=$normal*2.0**($logR-$gauge);
	print $fh $fields[0],$fields[1],$fields[2],$fields[3],$totalR,"\n";
}

}

close($fh); 

