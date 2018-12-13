#!/usr/bin/perl
use strict;
use PerlIO::gzip;

my $file1 = shift; #pre-tRNAs bed12
my $file2 = shift; #sam
my $file3 = shift; #fq
my $flank = shift || 50; #length of added flanking region

open FILE1, "<:gzip(autopop)", "$file1" or die "can t open $file1\n";
open FILE2, "<:gzip(autopop)", "$file2" or die "can t open $file2\n";
open(FILE3, "gunzip -c $file3|") or die "Cannot open $file3\n";

my %tRNA  = ();
my %read  = ();
my @entry = ();

my $flankEnd = $flank - 6; # consider CCACCA tail 

while(<FILE1>){
	chomp $_;	
	my @id = split/\t/,$_;
	$tRNA{$id[3]} = $id[10];  #length
}

while(<FILE2>){
	chomp $_;

	my @id = split/\t/,$_;
	my $start  = $id[3];
	my $len    = length($id[9]); 
	my @cigar = split(/(I|M|X|D)/, $id[5]);
	my $I = 0;
	my $D = 0;
	for(my $i=0; $i < scalar(@cigar); $i++){
		if($cigar[$i] eq "I"){ $I += $cigar[$i-1];}
		if($cigar[$i] eq "D"){ $D += $cigar[$i-1];}
	}
	my $end    = $start + $len -1 + $D - $I;

	if($start > $flank && $end <= $tRNA{$id[2]}-$flankEnd){
		$read{$id[0]} = 0;
	}	
}

while(<FILE3>){
	chomp;
	push @entry, $_;

	if (scalar(@entry) == 4) {
		my ($id, $seq, $plusLine, $qual) = @entry;
		@entry = ();
		my @ID = split/\@/,$id;
		if(exists $read{$ID[1]}){
			print "$id\n$seq\n$plusLine\n$qual\n";
		}
	}
}
