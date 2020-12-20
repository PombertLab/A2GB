#!/usr/bin/perl
## Converts the output of a tblastn search (outfmt 6) to the proper GFF3 format for loading into webApollo
## Requires the product list from the NCBI faa file (obtained with getProducts.pl)

use strict;
use warnings;

my $usage = 'USAGE = perl TBLASTN_to_GFF3.pl *.tblastn';
die "$usage\n" unless @ARGV;

while (my $file = shift@ARGV){
	open BLAST, "<$file";
	$file =~ s/.tblastn.6$//;
	$file =~ s/.blastn.6$//;
	open PROD, "<$file.products"; ## product list obtained with getProducts.pl
	open OUT, ">$file.gff";
	my $hit = 0;
	my %products = ();
	## Filling the products database
	while (my $line = <PROD>){
		chomp $line;
		if ($line =~ /^(\S+)\t(.*)$/){
			$products{$1}=$2;
		}
	}
	
	## Converting BLAST to GFF3
	while (my $line = <BLAST>){
		chomp $line;
		## Discarding comments, if any
		if ($line =~ /^\#/){next;}
		## Printing genes
		elsif ($line =~ /^(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)\t(\S+)/){
			my $query = $1; 	## protein accession number
			my $target = $2;	## location of the hit (contig or chromosome)
			my $identity = $3;  ## identity %
			my $len = $4;		## alignment length
			my $mis = $5;		## mismatches
			my $gap = $6;		## gaps
			my $qstart = $7;	## query start
			my $qend = $8;		## query end
			my $tstart = $9;	## target start
			my $tend = $10;		## target end
			my $evalue = $11;	## evalue
			my $bit = $12;		## bit score
			$hit++;
			## Forward strand
			if ($tstart < $tend){
				print OUT "$target"."\t"."TBLASTN"."\t"."match"."\t"."$tstart"."\t"."$tend"."\t"."$evalue"."\t".'+'."\t".'.'."\t"."ID=hit$hit".';'."Name=hit$hit".';'."Note=$query".':'."$products{$query}"."\n";
				print OUT "$target"."\t"."TBLASTN"."\t"."match_part"."\t"."$tstart"."\t"."$tend"."\t"."$evalue"."\t".'+'."\t"."0"."\t"."gene_id=hit$hit".';'."Parent=hit$hit".';'."transcript_id=hit$hit.t1".';'."Note=$query".':'."$products{$query}"."\n";
			}
			elsif ($tstart > $tend){
				print OUT "$target"."\t"."TBLASTN"."\t"."match"."\t"."$tend"."\t"."$tstart"."\t"."$evalue"."\t".'-'."\t".'.'."\t"."ID=hit$hit".';'."Name=hit$hit".';'."Note=$query".':'."$products{$query}"."\n";
				print OUT "$target"."\t"."TBLASTN"."\t"."match_part"."\t"."$tend"."\t"."$tstart"."\t"."$evalue"."\t".'-'."\t"."0"."\t"."gene_id=hit$hit".';'."Parent=hit$hit".';'."transcript_id=hit$hit.t1".';'."Note=$query".':'."$products{$query}"."\n";
			}
		}
	}
}
