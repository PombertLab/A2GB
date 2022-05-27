#!/usr/bin/perl
my $name = 'BLAST_to_GFF3.pl';
my $version = '0.3';
my $updated = '2021-04-09';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

## Defining options
my $usage = << "OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Converts the output of blastn/tblastn searches (outfmt 6) to GFF format for loading into Apollo

USAGE		${name} \\
		  -b *.tblastn \\
		  -p product_list.tsv \\
		  -t tblastn

-b (--blast)	BLAST output file(s) in outfmt 6 format
-p (--products)	Tab-delimited list of queries and their products (can be generated with getProducts.pl)
-t (--type)	BLAST type: blastn or tblastn [Default: blastn]
-o (--outdir)	Output directory [Default: ./]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @blast;
my $products;
my $type = 'blastn';
my $outdir = './';
GetOptions(
	'b|blast=s@{1,}' => \@blast,
	'p|products=s' => \$products,
	't|type=s' => \$type,
	'o|outdir=s' => \$outdir
);

## Creating output directory
unless (-d $outdir){
	mkdir ($outdir,0755) or die "Can't create $outdir: $!\n";
}

## Setting up variables
my %products;
my $query;
my $target;
my $identity;
my $len;
my $mis;
my $gap;
my $qstart;
my $qend;
my $tstart;
my $tend;
my $evalue;
my $bit;
my $hit;

## Workign on BLAST output files
while (my $blast = shift@blast){

	open BLAST, "<", "$blast" or die "Can't open $blast: $!\n";
	open PROD, "<", "$products" or die "Can't open $products: $!\n";

	my ($basename) = fileparse($blast);
	$basename =~ s/.tblastn.6$//;
	$basename =~ s/.blastn.6$//;
	open GFF, ">", "$outdir/$basename.gff" or die "Can't create $outdir/$basename.gff: $!\n";
	
	$hit = 0;
	
	## Filling the products database
	%products = ();
	while (my $line = <PROD>){
		chomp $line;
		if ($line =~ /^(\S+)\t(.*)$/){
			my $locus = $1;
			my $product = $2;
			$products{$locus} = $product;
		}
	}
	
	## Converting BLAST to GFF3
	while (my $line = <BLAST>){
		chomp $line;

		## Discarding comments, if any
		if ($line =~ /^\#/){ next; }

		## Printing matches
		else {
			my @columns = split ("\t", $line); 
			$query = $columns[0];    ## protein accession number
			$target = $columns[1];   ## location of the hit (contig or chromosome)
			$identity = $columns[2]; ## identity %
			$len = $columns[3];      ## alignment length
			$mis = $columns[4];      ## mismatches
			$gap = $columns[5];      ## gaps
			$qstart = $columns[6];   ## query start
			$qend = $columns[7];     ## query end
			$tstart = $columns[8];   ## target start
			$tend = $columns[9];     ## target end
			$evalue = $columns[10];  ## evalue
			$bit = $columns[11];     ## bit score
			$hit++;

			data_print();
		}
	}
}
## Subroutines
sub data_print {
	my $start;
	my $end;
	my $strand;

	## Forward strand
	if ($tstart < $tend){
		$start = $tstart;
		$end = $tend;
		$strand = '+';
	}

	## Reverse strand
	else {
		$start = $tend;
		$end = $tstart;
		$strand = '-';
	}

	print GFF "$target"."\t"."$type"."\t"."match"."\t"."$start"."\t"."$end"."\t"."$evalue"."\t"."$strand";
	print GFF "\t".'.'."\t"."ID=hit$hit".';'."Name=hit$hit".';'."Note=$query".':'."$products{$query}"."\n";

	print GFF "$target"."\t"."$type"."\t"."match_part"."\t"."$start"."\t"."$end"."\t"."$evalue"."\t"."$strand";
	print GFF "\t"."0"."\t"."gene_id=hit$hit".';'."Parent=hit$hit".';'."transcript_id=hit$hit.t1".';';
	print GFF "Note=$query".':'."$products{$query}"."\n";
}