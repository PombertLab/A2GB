#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'annot_multifasta.pl';
my $version = '0.1';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Adds predicted functions FASTA headers
USAGE		annot_multifasta.pl -a BEOM2.annotations.curated -f BEOM2.proteins.fasta -s 'Hamiltosporidium magnivora' -i BE-OM-2

-a	Annotations list (tab-delimited)
-f	Protein/mRNA multifasta file
-s	Species
-i	Isolate

OPTIONS
die "$options" unless @ARGV;

my $annots; my $fasta;
my $species; my $isolate;
GetOptions(
	'a=s' => \$annots,
	'f=s' => \$fasta,
	's=s' => \$species,
	'i=s' => \$isolate
);

## Creating DB of products
my %prod = ();
open PROD, "<$annots";
while (my $line = <PROD>){
	chomp $line;
	my @cols = split("\t", $line);
	$prod{$cols[0]} = $cols[1];
}

## Working on FASTA
open FAS, "<$fasta";
$fasta =~ s/.fasta//;
open TMP, ">$fasta.annots.fasta";
while (my $line = <FAS>){
	chomp $line;
	if ($line =~ /^>(\S+)\s+(.*)$/){
		print TMP ">$1 \[organism=$species\]\[isolate=$isolate\]\[product=$prod{$1}\]$2\n";}
	else {print TMP "$line\n";}
}