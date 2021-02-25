#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'add_metadata_to_fasta.pl';
my $version = '0.3a';
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
-c (--chromosome)	Tab-delimited contig/chromosome assignment file

OPTIONS
die "$usage\n" unless @ARGV;

my @fasta; my $chromosomes;

## NCBI Fasta headers 
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
	'c|chromosome=s' => \$chromosomes
);
die("[E] Organism required.\n") unless(exists(${meta{"organism"}}));
die("[E] Fasta files required.\n") unless(@fasta);

## Populating database of contigs and their assigned chromosomes, if desired
my %chromo;
if ($chromosomes){
	open CHR, "<", "$chromosomes" or die "Can't open chromosome file $chromosomes\n";
	while (my $line = <CHR>){
		chomp $line;
		if ($line =~ /^(\S+)\s+(.*)$/){
			my $contig = $1;
			my $chromo_assig = $2;
			$chromo_assig =~ s/\s+$//; ## Removing trailing spaces, if any
			$chromo{$contig} =  $chromo_assig;
		}
	} 
}

## Working on FASTA files
while (my $file = shift@fasta){
	open IN, "<$file";
	open OUT, ">$file.headers";
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			my $contig = $1;
			print OUT ">$contig";
			for my $key (keys %meta){
				if($meta{$key}){
					print OUT " [$key=$meta{$key}]";
				}
			}
			if ($chromosomes){
				if (exists $chromo{$contig}){
					print OUT " [chromosome=$chromo{$contig}]";
				}
			}
			print OUT "\n";
		}
		else {print OUT "$line\n";}
	}
	close OUT;
	system "mv $file.headers $file"; ## Overwrites original file 
}
