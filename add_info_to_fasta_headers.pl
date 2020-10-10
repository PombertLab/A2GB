#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'add_info_to_fasta_headers.pl';
my $version = '0.1';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	XXX
		
USAGE		$NAME *.fasta

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
my $organism;
my $strain;
my $lineage;
my $gcode = 1;
my $moltype = 'genomic';
my $chromosome = 0;

GetOptions(
	'o|organism=s' => \$organism,
	's|strain=s' => \$strain,
	'l|lineage=s' => \$lineage,
	'g|gcode=i' => \$gcode,
	'm|moltype=s' => \$moltype,
	'c|chromosome=i' => \$chromosome
);

while (my $file = shift@ARGV){
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
