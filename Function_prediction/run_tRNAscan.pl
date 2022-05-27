#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'run_tRNAscan.pl';
my $version = '0.2b';
my $updated = '2022-05-27';

use strict;
use warnings;
use File::Basename;
use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Runs tRNAscan-SE on FASTA file(s) using Infernal

USAGE		${name} \\
		  -f *.fasta \\
		  -t E \\
		  -d tRNA/

OPTIONS:
-f (--fasta)	FASTA file(s)
-t (--type)	tRNA type: eukaryotic (E); bacterial (B); archaeal (A) [Default: E]
-d (--dir) 	Output directory [Default: ./]
OPTIONS
die "\n$usage\n" unless @ARGV;

my @fasta;
my $type = 'E';
my $odir = './';
GetOptions(
	'f|fasta=s@{1,}' => \@fasta,
	't|type=s' => \$type,
	'd|dir=s' => \$odir
);

$type = uc($type); ## Making sure that the letters are uppercase
unless ($type =~ /[ABE]/){die "Unrecognized tRNA type: $type. Please use A, B or E...\n";}

## Creating output directory
unless (-d $odir){
	mkdir ($odir,0755) or die "Can't create folder $odir: $!\n";
}
print "\nOutput files will be located in directory $odir\n";

while (my $file = shift@fasta){
	my ($fasta, $dir) = fileparse($file);
	print "Working on file $fasta located in $dir\n";
	system ("tRNAscan-SE \\
	   -$type $file \\
	   -o ${odir}/$fasta.tRNAs \\
	   -f ${odir}/$fasta.struct \\
	   -m ${odir}/$fasta.stats") == 0 or checksig();
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