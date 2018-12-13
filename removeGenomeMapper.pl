#!/usr/bin/env perl

use strict;
use PerlIO::gzip;

my $file1 = shift; #pre-tRNAs
my $file2 = shift; #sam
my $out   = shift; #filtered sam
my $flank = shift || 50;

open FILE1, "<:gzip(autopop)", "$file1" or die "can t open $file1\n";
open FILE2, "<:gzip(autopop)", "$file2" or die "can t open $file2\n";
open OUT,   "> $out" or die "can t open $out\n";

my %ID   = ();
my %SEQ  = ();
my @SAM  = ();


while(<FILE1>){
        chomp $_;

	if($_=~ m/>/){
		my @id = split/>/,$_;
		$ID{$id[1]} = 0;
	}
}

while(<FILE2>){
	chomp $_;
	
	my @LINE = split/\t/,$_;

	if($LINE[0] =~ m/^@/){
		print OUT $_."\n";
		next;
	}
	
	if(length $LINE[9] >= $flank || length $LINE[9] >= 30){ # Filtering for read length longer than flank or 30
		
		if(!exists $SEQ{$LINE[0]}){
			$SEQ{$LINE[0]} = 0;
		}
		
		if(!exists $ID{$LINE[2]}){
			$SEQ{$LINE[0]} = 1;
		}
	}
	else{
		if(!exists $SEQ{$LINE[0]}){
			$SEQ{$LINE[0]} = 0;
		}
	}
	push @SAM,$_;
}


foreach (@SAM){
	my @LINE = split/\t/,$_;
	if($SEQ{$LINE[0]} eq 0){
		print OUT $_."\n";
	}
}
