#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'add_metadata_to_fasta.pl';
my $version = '0.3';
my $updated = '02/24/21';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	This script adds metadata to fasta headers. This metadata is required for submission to NCBI GenBank.
		
USAGE		${name} -o 'Chloropicon primus RCC138' -s RCC138 -g 1 -f *.fasta

OPTIONS:
-f (--fasta)		Specifies which FASTA files to add metadata to
-o (--organism)		Full organism name; e.g. 'Chloropicon primus RCC138'
-s (--strain)		Strain definition; e.g. RCC138
-i (--isolate)		Isolate name; e.g. 'Pacific Isolate'
-l (--lineage)		NCBI taxonomic lineage; e.g. 'cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;'
-g (--gcode)		NCBI genetic code [Default: 1]
-m (--moltype)		NCBI moltype descriptor [Default: genomic]

OPTIONS
die "$usage\n" unless @ARGV;

## NCBI Fasta headers 
my @fasta;
my %meta = (
	"organism" => undef,
	"strain" => undef,
	"isolate" => undef,
	"lineage" => undef,
	"gcode" => 1,
	"moltype" => 'genomic',
);

GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|organism=s' => \$meta{"organism"},
	's|strain=s' => \$meta{"strain"},
	'i|isolate=s' => \$meta{"isolate"},
	'l|lineage=s' => \$meta{"lineage"},
	'g|gcode=i' => \$meta{"gcode"},
	'm|moltype=s' => \$meta{"moltype"},
);
die("[E] Organism required.\n") unless(exists(${meta{"organism"}}));

while (my $file = shift@fasta){
	open IN, "<$file";
	open OUT, ">$file.headers";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			print OUT ">$1";
			for my $key (keys %meta){
				if($meta{$key}){
					print OUT " [$key=$meta{$key}]";
				}
			}
			print OUT "\n";
		}
		else {print OUT "$line\n";}
	}
	close OUT;
	system "mv $file.headers $file"; ## Overwrites original file 
}
