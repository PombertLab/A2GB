#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'check_problems.pl';
my $version = '0.2';

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Checks multifasta protein files for problems to help detect annotation errors:
		- Proteins interrupted by stop codons
		- Proteins that do not start with methionines
		
USAGE		$name -s -m -f *.prot

OPTIONS:
-s (--stop)	Checks for stop codons
-m (--meth)	Checks for missing start methionine
-f (--fasta)	FASTA amino acids input file(s)
OPTIONS
die "$usage\n" unless @ARGV;

my $stop; my $meth; my @fasta;
GetOptions(
	's|stop' => \$stop,
	'm|meth' => \$meth,
	'f|fasta=s@{1,}' => \@fasta
);
my $count = 0;
if ($stop){ ## Stop codons
	foreach my $file (@fasta){
		$count = 0;
		my ($fasta, $dir) = fileparse($file);
		print "\nChecking for internal stop codons in $fasta located in $dir\n";
		my $locus; my %prot;
		open IN, "<$file";
		while (my $line = <IN>){
			chomp $line;
			if ($line =~ /^>(\S+)/){$locus = $1;}
			else {$prot{$locus} .= $line;}
		}
		for (keys %prot){
			my $key = $_;
			my $seq = $prot{$_};
			if ($seq =~ /\*/){
				print "Protein $key contains one or more stop codon(s) in $file\n";
				$count++;
			}
		}
	if ($count == 0){print "No internal stop codon found\n"}
	}
	print "\n";
}
if ($meth){ ## Start methionines
	print "\nChecking for missing start methionines:";
	foreach my $file (@fasta){
		$count = 0;
		my ($fasta, $dir) = fileparse($file);
		print "\nChecking for missing start methionines in $fasta located in $dir\n";
		my $locus; my %prot;
		open IN, "<$file";
		while (my $line = <IN>){
			chomp $line;
			if ($line =~ /^>(\S+)/){$locus = $1;}
			else {$prot{$locus} .= $line;}
		}
		for (keys %prot){
			my $key = $_;
			my $seq = $prot{$_};
			if ($seq =~ /^([^M])/){
				print "Protein $key starts with $1 in $file\n";
				$count++;
			}
		}
	if ($count == 0){print "All proteins start with methionines...\n"}
	}
	print "\n";
}
