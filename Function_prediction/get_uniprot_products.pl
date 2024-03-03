#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'get_uniprot_products.pl';
my $version = '0.2c';
my $updated = '2024-03-03';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use File::Basename;
use PerlIO::gzip;

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Creates tab-delimited lists of products from UniProt trEMBL/SwissProt FASTA files

USAGE		${name} \\
		  -f uniprot_*.fasta.gz \\
		  -o ./

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
unless (-d $odir) { 
	mkdir ($odir,0755) or die "\nError: Can't create $odir: $!\n";
}

## Workign on Fasta files
while (my $file = shift@fasta){

	my $fh;
	my $format;
	my $stime = time;
	my $basename;

	if ($file =~ /.gz$/){ ## Autodetecting if file is gzipped from the file extension
		open $fh, "<:gzip", $file or die "Can't open $file: $!\n";
		$format = 'gzip';
		$file =~ s/.fasta.gz$//;
		$basename = fileparse($file);
	}
	else{
		open $fh, "<", $file or die "Can't open $file: $!\n"; 
		$format = 'fasta';
		$file =~ s/.\w+$//;
		$basename = fileparse($file);
	}

	my $outfile = "$odir/$basename.list";
	open OUT, ">", $outfile or die "Can't create $outfile: $!\n";

	print "\nOutput directory = $odir\n";
	print "Extracting information from $file. This might take a while...\n\n";

	while (my $line = <$fh>){
		chomp $line;
		if ($line =~ /^>(\S+)\s(.*)\sOS/){
			my $locus = $1;
			my $product = $2;
			print OUT "$locus\t$product\n";
		}
	}

	if ($format eq 'gzip'){ binmode $fh, ":gzip(none)"; }
	my $runtime = time - $stime;
	print "Time to extract products from file $file: $runtime seconds.\n";

}
