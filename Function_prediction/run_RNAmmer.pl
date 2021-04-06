#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'run_RNAmmer.pl';
my $version = '0.5a';
my $updated = '2021-03-27'; 

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Finds ribosomal RNAs using RNAmmer 1.2+

USAGE		${name} \\
		  -k kingdom \\
		  -f *.fasta \\
		  -d rRNA/

OPTIONS:
-k (--kingdom)	Kingdom; arc, bac or euk [Default: euk]
-f (--fasta)	FASTA file(s) to annotate
-d (--dir) 	Output directory [Default: ./]
OPTIONS
die "\n$usage\n" unless @ARGV; 

my $kingdom = 'euk';
my @fasta;
my $odir = './';
GetOptions(
	'k|kingdom=s' => \$kingdom,
	'f|fasta=s@{1,}' => \@fasta,
	'd|dir=s' => \$odir
);

## Creating output directory
unless (-d $odir){
	mkdir ($odir,0755) or die "Can't create folder $odir: $!\n";
}
print "\nOutput files will be located in directory $odir\n";

## Iterating through files
while (my $file = shift@fasta){
	my ($fasta, $dir) = fileparse($file);
	print "Working on file $fasta located in $dir\n";
	system "rnammer \\
	   -S $kingdom \\
	   -m tsu,ssu,lsu \\
	   -gff ${odir}/$fasta.gff2 \\
	   -h ${odir}/$fasta.hmm \\
	   -f ${odir}/$fasta.rRNAs \\
	   < $file";
}
