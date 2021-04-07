#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'annot_multifasta.pl';
my $version = '0.2';
my $updated = '2021-04-07';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $options = <<"OPTIONS";

NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Adds metadata to FASTA headers from a tab-delimited list. Useful to  
		add annotation metadata to FASTA files of genes and/or proteins.

USAGE		${name} \\
		  -f proteins.fasta \\
		  -o proteins.annotated.fasta \\
		  -a annotations.tsv \\
		  -k 'species=Encephalitozoon hellem' 'strain=Swiss'

-f (--fasta)	Protein/mRNA multifasta file
-o (--output)	Output file name
-a (--annots)	Annotations list (tab-delimited)
-k (--keys)	Metadata key(s) to add

OPTIONS
die "$options" unless @ARGV;

my $fasta;
my $output;
my $annots;
my @keys;
GetOptions(
	'f|fasta=s' => \$fasta,
	'o|output=s' => \$output,
	'a|annots=s' => \$annots,
	'k|keys=s@{1,}' => \@keys
);

## Checking for options
unless ($fasta){
	die "\nERROR: Please provide an input FASTA file with -f\n\n";
}
unless ($output){
	die "\nERROR: Please provide an output name with -o\n\n";
}

## Creating DB of products, if provided
my %prod;
if ($annots){
	open PRODUCTS, "<", "$annots" or die "Can't open $annots: $!\n";
	while (my $line = <PRODUCTS>){
		chomp $line;
		my @cols = split("\t", $line);
		my $locus = $cols[0];
		my $description = $cols[1];
		$prod{$locus} = $description;
	}
}

## Working on metadata keys, if provided
my $metadata_header;
if (@keys){
	while (my $key = shift@keys){
		my ($metakey, $metadata) = split ("=", $key);
		my $entry = '['."$key".']';
		$metadata_header .= $entry;
	}
}
else{
	$metadata_header = '';
}

## Working on FASTA file
open FASTA, "<", "$fasta" or die "Can't open $fasta: $!\n";
$fasta =~ s/.\w+//;
open ANNOT, ">", "$output" or die "Can't write to $output: $!\n";

while (my $line = <FASTA>){
	chomp $line;

	## Working on header
	my $product;
	if ($line =~ /^>(\S+)\s+(.*)$/){
		my $locus = $1;
		my $remains = $2;

		if ($annots){
			if (exists $prod{$locus}){
				$product = '[product='."$prod{$locus}".']';
			}
			else{
				$product = '';
				print "product for $locus is missing!...\n";
			}
		}
		else {
			$product = '';
		}

		print ANNOT ">$locus ${product}${metadata_header}\n";}
	
	## Print sequence
	else {
		print ANNOT "$line\n";
	}
}
