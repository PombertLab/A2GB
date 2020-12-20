#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'get_uniprot_products.pl';
my $version = '0.2';

use strict; use warnings; use PerlIO::gzip;

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Creates tab-delimited lists of products from the trEMBL and SwissProt FASTA files
USAGE		get_uniprot_products.pl uniprot_sprot.fasta.gz uniprot_trembl.fasta.gz
OPTIONS
die "$usage\n" unless @ARGV;

while (my $file = shift@ARGV){
	my $fh; my $format;
	my $stime  = time;
	if ($file =~ /.gz$/){ ## Autodecting if file is gzipped from the file extension
		open $fh, "<:gzip", "$file" or die "Could not open $file for reading: $!\n";
		$format = 'gzip';
		$file =~ s/.fasta.gz$//;
	}
	else{
		open $fh, "<$file" or die "Could not open $file for reading: $!\n"; 
		$format = 'fasta';
		$file =~ s/.\w+$//;
	}
	
	open OUT, ">$file.list";
	print "Extracting information from $file. This might take a while...\n";
	while (my $line = <$fh>){
		chomp $line;
		if ($line =~ /^>(\S+)\s(.*)\sOS/){print OUT "$1\t$2\n";}
	}
	if ($format eq 'gzip'){binmode $fh, ":gzip(none)";}
	my $runtime = time - $stime;
	print "Time to extract products from file $file: $runtime seconds.\n";
}
