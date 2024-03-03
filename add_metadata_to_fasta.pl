#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'add_metadata_to_fasta.pl';
my $version = '0.4a';
my $updated = '2021-04-09';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	This script adds metadata to fasta headers. This metadata is required for submission to NCBI GenBank.
		
USAGE 1		${name} -f *.fasta -o 'Chloropicon primus RCC138' -s RCC138 -g 1	## Using CMD line switches
USAGE 2		${name} -f *.fasta -k metakeys_NCBI.tsv -c chromosomes.tsv		## Using metadata files

OPTIONS:
-f (--fasta)		Specifies which FASTA files to add metadata to

## Single metadata keys
-o (--organism)		Full organism name; e.g. 'Chloropicon primus RCC138'
-s (--strain)		Strain definition; e.g. RCC138
-i (--isolate)		Isolate name; e.g. 'Pacific Isolate'
-g (--gcode)		NCBI genetic code ## https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
-m (--moltype)		NCBI moltype descriptor (e.g. genomic)

## Metadata files
-k (--keys)		Tab-delimited NCBI metadata key -> value file
-c (--chromosome)	Tab-delimited contig -> chromosome assignment file
OPTIONS
die "$usage\n" unless @ARGV;

my @fasta;
my $chromosomes;
my $metakeys;

## NCBI Fasta headers 
my %meta = (
	"organism" => undef,
	"strain" => undef,
	"isolate" => undef,
	"gcode" => undef,
	"moltype" => undef,
);

GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'k|keys=s' => \$metakeys,
	'c|chromosome=s' => \$chromosomes,
	'o|organism=s' => \$meta{"organism"},
	's|strain=s' => \$meta{"strain"},
	'i|isolate=s' => \$meta{"isolate"},
	'g|gcode=i' => \$meta{"gcode"},
	'm|moltype=s' => \$meta{"moltype"},
);
die "[E] Fasta files required.\n" unless @fasta;

## Populating database of metadata keys and their values, if desired
if ($metakeys){
	open META, '<', $metakeys or die "Can't open metadata file $metakeys\n";
	while (my $line = <META>){
		chomp $line;
		if ($line =~ /^#/){ next; } ## Ignoring comments
		elsif ($line =~ /^(\S+)\s+(.*)$/){
			my $metakey = $1;
			my $metavalue = $2;
			$metavalue =~ s/\s+$//; ## Removing trailing spaces, if any
			$meta{$metakey} =  $metavalue;
		}
	} 
}

## Populating database of contigs and their assigned chromosomes, if desired
my %chromo;
if ($chromosomes){
	open CHR, '<', $chromosomes or die "Can't open chromosome file $chromosomes\n";
	while (my $line = <CHR>){
		chomp $line;
		if ($line =~ /^#/){ next; } ## Ignoring comments
		elsif ($line =~ /^(\S+)\s+(.*)$/){
			my $contig = $1;
			my $chromo_assig = $2;
			$chromo_assig =~ s/\s+$//; ## Removing trailing spaces, if any
			$chromo{$contig} =  $chromo_assig;
		}
	} 
}

## Working on FASTA files
while (my $file = shift@fasta){

	open IN, '<', $file or die "Can't open $file: $!\n";
	open OUT, '>', "$file.headers" or die "Can't create $file.headers: $!\n";

	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^>(\S+)/){
			my $contig = $1;
			print OUT ">$contig";
			for my $key (keys %meta){
				if ($meta{$key}){
					print OUT " [$key=$meta{$key}]";
				}
			}
			if ($chromosomes){
				if (exists $chromo{$contig}){
					print OUT " [location=chromosome] [chromosome=$chromo{$contig}]";
				}
			}
			print OUT "\n";
		}
		else { print OUT "$line\n"; }
	}

	close OUT;

	system "mv $file.headers $file"; ## Overwrites original file 

}
