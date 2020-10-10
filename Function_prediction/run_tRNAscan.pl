#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'run_tRNAscan.pl';
my $version = '0.1';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
SYNOPSIS	Runs tRNAscan-SE on FASTA file(s)
USAGE		$name -f *.fasta
OPTIONS:
-f (--fasta)	## FASTA file(s)
OPTIONS
die "$usage\n" unless @ARGV;

my @fasta;
GetOptions('f|fasta=s@{1,}' => \@fasta);

while (my $fasta = shift@fasta){
	system "tRNAscan-SE $fasta > $fasta.tRNAs";
}
