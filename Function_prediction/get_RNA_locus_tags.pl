#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'get_RNA_locus_tags.pl';
my $version = '0.1'; 

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
SYNOPSIS	Creates tab-delimited lists of locus_tags and corresponding RNAs

USAGE		$name -f features.list -ti *.tRNAs -ri *.gff3 -to tRNA.annotations -ro rRNA.annotations

OPTIONS:
-f	Feature list generated by ApolloGFF3toEMBL.pl [Default: features.list]
-ti	tRNAscan input files (.tRNAs) produced by run_tRNAscan.pl
-ri	RNAmmer input files (.gff3) produced by run_RNAmmer.pl + RNAmmer_to_GFF3.pl
-to 	tRNA output file [Default: tRNA.annotations]
-ro 	rRNA output file [Default: rRNA.annotations]
OPTIONS
die "\n$usage\n" unless @ARGV;

my $feat = 'features.list';
my @ti;
my @ri;
my $to = 'tRNA.annotations';
my $ro = 'rRNA.annotations';
GetOptions(
	'f=s' => \$feat,
	'ti=s@{1,}' => \@ti,
	'ri=s@{1,}' => \@ri,
	'to=s' => \$to,
	'ro=s' => \$ro
);

## Creating database of features
open FEAT, "<", "$feat" or die "Can't open features list: $feat\n";
my %features;
while (my $line = <FEAT>){
	chomp $line;
	unless ($line =~ /^#/){
		my @cols = split("\t", $line);
		if (($cols[2] eq "rRNA") or ($cols[2] eq "tRNA")){
			#print "$line\n";
			for ($cols[4]..$cols[5]){
				$features{$cols[2]}{$cols[1]}{$_} = $cols[0];
			}
		}
	}
}

## Working on rRNA files
if (@ri){
	open ROUT, ">", "$ro" or die "Can't write to rRNA annotations to file: $ro\n";
	while (my $ri = shift @ri){
		open RRNA, "<", "$ri" or die "Can't open .rRNAs file: $ri\n";
		while (my $line = <RRNA>){
			chomp $line;
			my @cols = split("\t", $line);
			my ($rrna) = $cols[8] =~ /Note=(.*)$/;
			$rrna =~ s/_/ /g; $rrna =~ s/rRNA/ribosomal RNA/g;
			if (exists $features{$cols[2]}{$cols[0]}{$cols[3]}){
				print ROUT "$features{$cols[2]}{$cols[0]}{$cols[3]}"."\t"."$rrna"."\n";
			}
		}
	}
}

## Working on tRNA files
if (@ti){
	open TOUT, ">", "$to" or die "Can't write to tRNA annotations to file: $to\n";
	while (my $ti = shift @ti){
		open TRNA, "<", "$ti" or die "Can't open .tRNAs file: $ti\n";
		while (my $line = <TRNA>){
			chomp $line;
			unless ($line =~ /^(Sequence|Name|--------)/){
				my @cols = split("\t", $line);
				for (0..$#cols){ $cols[$_] =~ s/\s+$//;} ## Getting rid of padding spaces
				$cols[5] =~ tr/Tt/Uu/;
				my $tRNA = 'tRNA-'."$cols[4]"."\($cols[5]\)";
				if (exists $features{'tRNA'}{$cols[0]}{$cols[2]}){
					print TOUT "$features{'tRNA'}{$cols[0]}{$cols[2]}"."\t"."$tRNA"."\n";
				}
			}
		}
	}
}
