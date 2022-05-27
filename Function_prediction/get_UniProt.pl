#!/usr/bin/perl
## Pombert Lab, 2020
my $name = 'get_UniProt.pl';
my $version = '0.2e';
my $updated = '2022-05-27';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Downloads the SwissProt and/or trEMBL databases from UniProt

EXAMPLE		${name} \\
		  -s \\
		  -t \\
		  -f ./ \\
		  -n 20 \\
		  -l download.log 

OPTIONS:
-s (--swiss)		Download Swiss-Prot
-t (--trembl)		Download trEMBL
-f (--folder)		Download folder [Default: ./]
-n (--nice)		Linux Process Priority [Default: 20] ## Runs downloads in the background
-l (--log)		Print download information to log file
-d (--decompress)	Decompresss downloaded files with gunzip ## trEMBL files will be huge, off by default
OPTIONS
die "\n$usage\n" unless @ARGV;

my $nice = 20;
my $swiss;
my $trembl;
my $folder = './';
my $log;
my $decomp;
GetOptions(
	'n|nice=i' => \$nice,
	's|swiss' => \$swiss,
	't|trembl' => \$trembl,
	'f|folder=s' => \$folder,
	'l|log=s' => \$log,
	'd|decompress' => \$decomp
);

## Checking output directory + creating download log
unless (-d $folder){ 
	mkdir ($folder,0755) or die "Can't create $folder: $!\n";
}
print "\nOutput files will be located in directory $folder\n";
if ($log){ open LOG, ">", "${folder}/${log}"; }

my $url = "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete";

## Downloading SwissProt
if ($swiss){
	print "\nDownloading SwissProt...\n\n";
	if ($log){
		my $date = `date`;
		print LOG "Downloading SwissProt on $date";
	}

	system ("nice \\
	  -n $nice \\
	  wget \\
	  -c \\
	  -P $folder \\
	  $url/uniprot_sprot.fasta.gz") == 0 or checksig();

	if ($log){
		my $size = `du -h $folder/uniprot_sprot.fasta.gz`;
		print LOG "$size";
	}
}

## Downloading trEMBL
if ($trembl){
	print "\nDownloading trEMBL database. This will take a while...\n\n";
	if ($log){
		my $date = `date`;
		print LOG "Downloading trEMBL on $date";
	}

	system ("nice \\
	  -n $nice \\
	  wget \\
	  -c \\
	  -P $folder \\
	  $url/uniprot_trembl.fasta.gz") == 0 or checksig();
	
	if ($log){
		my $size = `du -h $folder/uniprot_trembl.fasta.gz`;
		print LOG "$size";
	}
}

## Decompressing the downloaded databases
if ($decomp){
	if ($swiss){
		print "\nDecompressing downloaded SwissProt database with gunzip...\n\n";
		system "gunzip $folder/uniprot_sprot.fasta.gz";
	}
	if ($trembl){
		print "\nDecompressing downloaded trEMBL database with gunzip. Will take a while...\n\n";
		system "gunzip $folder/trembl.fasta.gz";
	}
}
## Updating logs
if ($log){
	my $end = `date`;
	print LOG "Finished downloading database(s) on $end";
}

### Subroutine(s)
sub checksig {

	my $exit_code = $?;
	my $modulo = $exit_code % 255;

	print "\nExit code = $exit_code; modulo = $modulo \n";

	if ($modulo == 2) {
		print "\nSIGINT detected: Ctrl+C => exiting...\n";
		exit(2);
	}
	elsif ($modulo == 131) {
		print "\nSIGTERM detected: Ctrl+\\ => exiting...\n";
		exit(131);
	}

}
