#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'add_metadata_to_fasta.pl';
my $version = '0.1';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	This script adds metadata to fasta headers. This metadata is required for submission to NCBI GenBank.
		
USAGE		$NAME -o 'Chloropicon primus RCC138' -s RCC138 -g 1 -f *.fasta

OPTIONS:
-o (--organism)	Full organism name; e.g. 'Chloropicon primus RCC138'
-s (--strain)	Strain definition; e.g. RCC138
-l (--lineage)	NCBI taxonomic lineage; e.g. 'cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;'
-g (--gcode)	NCBI genetic code [Default: 1]
-m (--moltype)	NCBI moltype descriptor [Default: genomic]
-c (--chromosome)	Annotate contigs as chromosomes

OPTIONS
die "$usage\n" unless @ARGV;

## NCBI Fasta headers 
my @fasta;
my $organism;
my $strain;
my $lineage;
my $gcode = 1;
my $moltype = 'genomic';
my $chromosome = 0;

GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|organism=s' => \$organism,
	's|strain=s' => \$strain,
	'l|lineage=s' => \$lineage,
	'g|gcode=i' => \$gcode,
	'm|moltype=s' => \$moltype,
	'c|chromosome=i' => \$chromosome
);

while (my $file = shift@fasta){
	open IN, "<$file";
	open OUT, ">$file.headers";
	$chromosome++;
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			print OUT ">$1 [organism=$organism][strain=$strain][lineage=$lineage][gcode=$gcode][moltype=$moltype][chromosome=$chromosome]\n";
		}
		else {print OUT "$line\n";}
	}
	close OUT;
	system "mv $file.headers $file"; ## Overwrites original file 
}
