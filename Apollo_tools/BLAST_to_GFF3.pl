#!/usr/bin/perl
my $name = 'BLAST_to_GFF3.pl';
my $version = '0.2a';
my $updated = '2021-04-02';

use strict; use warnings; use Getopt::Long qw(GetOptions);

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
-t (--type)		BLAST type: blastn or tblastn [Default: blastn]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @blast;
my $products;
my $type = 'blastn';
GetOptions(
	'b|blast=s@{1,}' => \@blast,
	'p|products=s' => \$products,
	't|type=s' => \$type
);

while (my $blast = shift@blast){
	open BLAST, "<", "$blast" or die "Can't open $blast: $!\n";
	$blast =~ s/.tblastn.6$//;
	$blast =~ s/.blastn.6$//;
	open PROD, "<", "$products" or die "Can't open $products: $!\n";
	open GFF, ">", "$blast.gff" or die "Can't create $blast.gff: $!\n";
	
	my $hit = 0;
	
	## Filling the products database
	my %products = ();
	while (my $line = <PROD>){
		chomp $line;
		if ($line =~ /^(\S+)\t(.*)$/){
			my $locus = $1;
			my $product = $2;
			$products{$locus}=$product;
		}
	}
	
	## Converting BLAST to GFF3
	while (my $line = <BLAST>){
		chomp $line;
		if ($line =~ /^\#/){ ## Discarding comments, if any
			next;
		}

		else { ## Printing matches
			my @columns = split ("\t", $line); 
			my $query = $columns[0]; 	## protein accession number
			my $target = $columns[1];	## location of the hit (contig or chromosome)
			my $identity = $columns[2]; ## identity %
			my $len = $columns[3];		## alignment length
			my $mis = $columns[4];		## mismatches
			my $gap = $columns[5];		## gaps
			my $qstart = $columns[6];	## query start
			my $qend = $columns[7];		## query end
			my $tstart = $columns[8];	## target start
			my $tend = $columns[9];		## target end
			my $evalue = $columns[10];	## evalue
			my $bit = $columns[11];		## bit score
			$hit++;

			if ($tstart < $tend){ ## Forward strand
				print GFF "$target"."\t"."$type"."\t"."match"."\t"."$tstart"."\t"."$tend"."\t"."$evalue"."\t".'+';
				print GFF "\t".'.'."\t"."ID=hit$hit".';'."Name=hit$hit".';'."Note=$query".':'."$products{$query}"."\n";
				print GFF "$target"."\t"."$type"."\t"."match_part"."\t"."$tstart"."\t"."$tend"."\t"."$evalue"."\t".'+';
				print GFF "\t"."0"."\t"."gene_id=hit$hit".';'."Parent=hit$hit".';'."transcript_id=hit$hit.t1".';';
				print GFF "Note=$query".':'."$products{$query}"."\n";
			}

			elsif ($tstart > $tend){ ## Reverse strand
				print GFF "$target"."\t"."$type"."\t"."match"."\t"."$tend"."\t"."$tstart"."\t"."$evalue"."\t".'-';
				print GFF "\t".'.'."\t"."ID=hit$hit".';'."Name=hit$hit".';'."Note=$query".':'."$products{$query}"."\n";
				print GFF "$target"."\t"."$type"."\t"."match_part"."\t"."$tend"."\t"."$tstart"."\t"."$evalue"."\t".'-';
				print GFF "\t"."0"."\t"."gene_id=hit$hit".';'."Parent=hit$hit".';'."transcript_id=hit$hit.t1".';';
				print GFF "Note=$query".':'."$products{$query}"."\n";
			}
		}
	}
}
