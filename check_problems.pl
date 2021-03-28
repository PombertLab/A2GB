#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'check_problems.pl';
my $version = '0.3a';
my $updated = '28/03/2021';

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Checks multifasta protein files for problems to help detect annotation errors:
		- Proteins interrupted by stop codons
		- Proteins that do not start with methionines
		
USAGE		${name} -s -m -f *.prot

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
		print STDOUT "\nChecking for internal stop codons in $fasta located in $dir\n";

		my $locus;
		my %prot;

		open IN, "<", "$file" or die "Can't open $file: $!\n";
		while (my $line = <IN>){
			chomp $line;
			if ($line =~ /^>(\S+)/){$locus = $1;}
			else {$prot{$locus} .= $line;}
		}

		for (keys %prot){
			my $key = $_;
			my $seq = $prot{$_};
			if ($seq =~ /\*/){
				print STDERR "ERROR: Protein $key contains one or more stop codon(s) in $file\n";
				$count++;
			}
		}

		if ($count == 0){ print STDOUT "OK: No internal stop codon found\n"; }
	}
	print STDOUT "\n";
}
if ($meth){ ## Start methionines
	foreach my $file (@fasta){
		$count = 0;
		my ($fasta, $dir) = fileparse($file);
		print STDOUT "\nChecking for missing start methionines in $fasta located in $dir\n";

		my $locus;
		my %prot;

		open IN, "<", "$file" or die "Can't open $file: $!\n";
		while (my $line = <IN>){
			chomp $line;
			if ($line =~ /^>(\S+)/){$locus = $1;}
			else {$prot{$locus} .= $line;}
		}

		for (keys %prot){
			my $key = $_;
			my $seq = $prot{$_};
			if ($seq =~ /^([^M])/){
				print STDERR "ERROR: Protein $key starts with $1 in $file\n";
				$count++;
			}
		}

		if ($count == 0){ print STDOUT "OK: All proteins start with methionines...\n"; }
	}
	print STDOUT "\n";
}
