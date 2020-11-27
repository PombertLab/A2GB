#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'splitMakerGFF.pl';
my $version = '0.1'; 

use strict; use warnings;

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
SYNOPSIS	Splits the Augustus, GeneMark and Repeatmasker entries into distinct ones; useful for loading them as separate tracks in Apollo.
USAGE		splitMakerGFF.pl *.gff
OPTIONS
die "$usage\n" unless @ARGV;

while(my $file = shift@ARGV){
	open IN, "<$file";
	$file =~ s/.gff$//;
	open OUT1, ">$file.augustus.gff";
	open OUT2, ">$file.genemark.gff";
	open OUT3, ">$file.repeats.gff";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^\S+\taugustus_masked/){
			print OUT1 "$line\n";
		}
		elsif ($line =~ /^\S+\tgenemark/){
			print OUT2 "$line\n";
		}
		elsif ($line =~ /^\S+\trepeatmasker/){
			print OUT3 "$line\n";
		}
	}
}
