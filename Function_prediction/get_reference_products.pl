#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'get_reference_products.pl';
my $version = '0.1';
my $updated = '2021-04-08';

use strict; use warnings; use PerlIO::gzip; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Creates a single tab-delimited list of products NCBI protein (.faa) FASTA files

USAGE		${name} \\
		  -f *.gz \\
		  -l reference.list

OPTIONS:
-f (--fasta)	NCBI protein FASTA files 
-l (--list)	Desired output list name
OPTIONS
die "\n$usage\n" unless @ARGV;

my $list;
my @fasta;
GetOptions(
	'l|list=s' => \$list,
	'f|fasta=s@{1,}' => \@fasta
);

open OUT, ">", "$list" or die "Can't create $list: $!\n"; 
while (my $file = shift@fasta){

	my $fh;
	my $format;
	my $stime = time;

	if ($file =~ /.gz$/){ ## Autodecting if file is gzipped from the file extension
		open $fh, "<:gzip", "$file" or die "Can't open $file: $!\n";
		$format = 'gzip';
		$file =~ s/.faa.gz$//;
		$file =~ s/.fasta.gz$//;
		$file =~ s/.gz$//;
	}
	else {
		open $fh, "<", "$file" or die "Can't open $file: $!\n"; 
		$format = 'fasta';
		$file =~ s/.\w+$//;
	}

	print "Extracting information from $file. This might take a while...\n";
	while (my $line = <$fh>){
		chomp $line;
		if ($line =~ /^>(\S+)\s(.*)\s\[/){
			my $locus = $1;
			my $product = $2;
			print OUT "$locus\t$product\n";
		}
	}
	
	if ($format eq 'gzip'){ binmode $fh, ":gzip(none)"; }
	my $runtime = time - $stime;
	print "Time to extract products from file $file: $runtime seconds.\n";
}
