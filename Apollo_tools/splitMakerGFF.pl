#!/usr/bin/perl
## Pombert Lab, IIT, 2017
my $name = 'splitMakerGFF.pl';
my $version = '0.1a';
my $updated = '27/03/2021'; 

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Splits the Maker entries into distinct ones by source;
		Useful for loading them as separate tracks in Apollo.

USAGE		${name} \\
		  -g *.gff \\
		  -o SPLIT_MAKER \\
OPTIONS
die "\n$usage\n" unless @ARGV;

my @GFF;
my $odir;
GetOptions(
	'g|gff=s@{1,}' => \@GFF,
	'o|outdir=s' => \$odir,
);

## Creating output directory
unless (-d $odir){
	mkdir ($odir,0755) or die "Can't create folder $odir: $!\n";
}

while (my $gff = shift@GFF){
	open GFF, "<", "$gff" or die "Can't read $gff: $!\n";
	$gff =~ s/.gff$//;
	my ($basename) = fileparse($gff);

	my %database;
	while (my $line = <GFF>){
		chomp $line;
		if ($line =~ /^#/){ next; }
		elsif ($line =~ /^\S+\t\./) { next;}
		elsif ($line =~ /^\S+\t(\w+)/){
			my $source = $1;
			push (@{$database{$source}}, $line);
		}
	}

	foreach my $key (keys %database){
		my $filename = "$odir/$basename.$key.gff";
		open FH, ">", "$filename" or die "Can't create $filename: $!\n";
		while (my $data = shift @{$database{$key}}){
			print FH "$data\n";
		}
		close FH;
	}

}
