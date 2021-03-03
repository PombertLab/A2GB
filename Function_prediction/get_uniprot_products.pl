#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'get_uniprot_products.pl';
my $version = '0.2a';
my $updated = '03/03/2021';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename; use PerlIO::gzip;

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
SYNOPSIS	Creates tab-delimited lists of products from UniProt trEMBL/SwissProt FASTA files
UPDATED		${updated}

USAGE		${name} -f uniprot_*.fasta.gz -o ./

OPTIONS:
-f (--fasta)	FASTA files from Uniprot (.fasta or .fasta.gz)
-o (--output)	Output directory [Default: ./]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $odir = './';
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	'o|output=s' => \$odir
);

## Creating output folder, if needed 
unless (-e $odir) { mkdir $odir or die "\nError: Can't create folder $odir: $!\n\n"; }

## Workign on Fasta files
while (my $file = shift@fasta){
	my $fh; my $format;
	my $stime  = time;
	my $basename;
	if ($file =~ /.gz$/){ ## Autodetecting if file is gzipped from the file extension
		open $fh, "<:gzip", "$file" or die "Could not open $file for reading: $!\n";
		$format = 'gzip';
		$file =~ s/.fasta.gz$//;
		$basename = fileparse($file);
	}
	else{
		open $fh, "<", "$file" or die "Could not open $file for reading: $!\n"; 
		$format = 'fasta';
		$file =~ s/.\w+$//;
		$basename = fileparse($file);
	}
	
	open OUT, ">", "$odir/$basename.list";
	print "\nOutput directory = $odir\n";
	print "Extracting information from $file. This might take a while...\n\n";
	while (my $line = <$fh>){
		chomp $line;
		if ($line =~ /^>(\S+)\s(.*)\sOS/){
			my $locus = $1; my $product = $2;
			print OUT "$locus\t$product\n";
		}
	}
	if ($format eq 'gzip'){binmode $fh, ":gzip(none)";}
	my $runtime = time - $stime;
	print "Time to extract products from file $file: $runtime seconds.\n";
}
