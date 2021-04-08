#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'get_queries.pl';
my $version = '0.1';
my $updated = '2021-04-08';

use strict; use warnings;

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Generates simple lists of sequences present in FASTA files

USAGE		get_queries.pl *.fasta
OPTIONS
die "\n$usage\n" unless @ARGV;

while (my $fasta = shift@ARGV){

	open FASTA, "<", "$fasta" or die "Can't open $fasta: $!\n";
	$fasta =~ s/\.\w+$//;
	open OUT, ">", "$fasta.queries" or die "Can't create $fasta.queries: $!\n";

	while (my $line = <FASTA>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			my $query = $1;
			print OUT "$query\n";
		}
	}
} 