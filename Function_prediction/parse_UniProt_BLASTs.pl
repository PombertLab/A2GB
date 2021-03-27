#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'parse_UniProt_BLASTs.pl';
my $version = '0.3a';
my $updated = '27/03/2021';

use strict; use warnings; use Getopt::Long qw(GetOptions);

## Defining options
my $options = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Generates a tab-delimited list of products found/not found with BLAST/DIAMOND searches
REQUIREMENTS	BLAST/DIAMOND outfmt 6 format
		Tab-separated accession number/product list ## Can be created with get_uniprot_products.pl

USAGE		${name} \\
		  -b blast_output \\
		  -e 1e-10 \\
		  -q query.list \\
		  -u uniprot_list \\
		  -o parsed.tsv

OPTIONS:
-b (--blast)	BLAST/DIAMOND tabular output (outfmt 6)
-e (--evalue)	E-value cutoff [Default: 1e-10]
-q (--query)	List of proteins queried against UniProt
-u (--uniprot)	Tab-delimited list of UniProt accesssion numbers/products 
-o (--output)	Desired output name
OPTIONS
die "\n$options\n" unless @ARGV;

my $blast;
my $eval = '1e-10';
my $query;
my $uniprot;
my $output;
GetOptions(
	'b|blast=s' => \$blast,
	'e|evalue=s' => \$eval,
	'q|query=s' => \$query,
	'u|uniprot=s' => \$uniprot,
	'o|output=s' => \$output
);

open QUERIES, "<$query";
open PRODUCTS, "<$uniprot";
open BLAST, "<$blast";
open OUT, ">$output";

my %products = (); ## Empty product hash; WARNING will take a decent chunk of RAM for 
my %hits = (); ## Keeping track of BLAST hits
my %evalues = (); 

## Filling the product hash
print "\nLoading the products from $uniprot. I/O limited step; this might take a while...\n";
my $stime = time;
while (my $line = <PRODUCTS>){
	chomp $line;
	if ($line =~ /^(\S+)\t(.*)$/){$products{$1}=$2;}
}
my $runtime = time - $stime;
print "Time to load $uniprot: $runtime seconds\n";
close PRODUCTS;

## Parsing BLAST
print "\nParsing BLAST/DIAMOND outfmt 6 file...\n";
my $ptime = time;
while (my $line = <BLAST>){
	chomp $line;
	my @columns = split("\t", $line);
	my $query = $columns[0];
	my $hit = $columns[1];
	my $evalue = $columns[10];
	if (exists $hits{$query}){next;}
	elsif ($products{$hit} =~ /uncharacterized/i){next;} ## Discarding uninformative BLAST hists
	elsif ($products{$hit} =~ /hypothetical/i){next;} ## Discarding uninformative BLAST hists
	elsif ($products{$hit} =~ /predicted protein/i){next;} ## Discarding uninformative BLAST hists
	else{
		if ($evalue <= $eval){ ## Checking against desired E-value cutoff
			$hits{$query}=$products{$hit};
			$evalues{$query}=$evalue;
		}
	}
}
close BLAST;

## Writing parsed list
print OUT "\###Protein\tE-value\tProduct\n";
while (my $line = <QUERIES>){
	chomp $line;
	if ($line =~ /^(\S+)/){
		my $query = $1;
		if (exists $hits{$query}){print OUT "$query\t$evalues{$query}\t$hits{$query}\n";}
		else{print OUT "$query\tN\/A\thypothetical protein\n";}
	}
}
my $pruntime = time - $stime;
print "Time to parse $blast: $pruntime seconds\n\n";
close OUT;
close QUERIES;