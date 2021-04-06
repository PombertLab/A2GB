#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'run_InterProScan.pl';
my $version = '0.1';
my $updated = '2021-03-27';

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs InterProScan 5 on FASTA file(s)

USAGE		${name} \\
		  -c 10 \\
		  -f *.fasta \\
		  -ip \\
		  -go \\
		  -pa \\
		  -d Interproscan/ \\
		  -l interproscan.log

OPTIONS:
-c (--cpu)		Number of CPU cores to use [Default: 10]
-f (--fasta)		FASTA file(s) to query
-ip (--iprlookup)	Use InterPro's pre-calculated match lookup service
-go (--goterms)		Gene ontology search (requires --iprlookup)
-pa (--pathways)	KEGG pathways (requires --iprlookup)
-d (--dir)		Output directory [Default: ./]
-l (--log)		Log name [Default: interproscan.log]
OPTIONS
die "\n$usage\n" unless @ARGV;

my $cpu = 10;
my @fasta;
my $odir = './';
my $ipr;
my $go;
my $pa;
my $log = 'interproscan.log';
GetOptions(
	'c|cpu=i' => \$cpu,
	'f|fasta=s@{1,}' => \@fasta,
	'd|dir=s' => \$odir,
	'ip|iprlookup' => \$ipr,
	'go|goterms' => \$go,
	'pa|pathways' => \$pa,
	'l|log=s' => \$log
);

my $iprlookup; my $goterms; my $pathways;
if (defined $ipr){$iprlookup = '-iprlookup';}
if ((defined $go) and (defined $ipr)){$goterms = '-goterms';}
if ((defined $pa) and (defined $ipr)){$pathways = '-pa';}

## Creating output directory
unless (-d $odir){
	mkdir ($odir,0755) or die "Can't create folder $odir: $!\n";
}
print "\nOutput files will be located in directory $odir\n";
open LOG, ">", "${odir}/$log" or die "Can't create file ${odir}/$log: $!\n";

## Working on fasta files
while (my $file = shift@fasta){
	
	my ($fasta, $dir) = fileparse($file);
	print "Working on file $fasta located in $dir\n";
	
	my $sdate = `date`;
	print LOG "InterProScan search on $fasta started on: $sdate\n";
	
	system "interproscan.sh \\
	   -cpu $cpu \\
	   -i $file \\
	   $iprlookup \\
	   $goterms \\
	   $pathways \\
	   -b ${odir}/$fasta.interpro";
	
	my $edate = `date`;
	print LOG "InterProScan search on $fasta completed on: $edate\n";
}
