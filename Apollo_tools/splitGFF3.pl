#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'splitGFF3.pl';
my $version = '0.3a';
my $updated = '27/03/2021';

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

## Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Splits an Apollo-like GFF3 file into distinct GFF3 (.gff3) and FASTA (.fsa) files, one per contig/chromosome

USAGE		${name} \\
		  -g file.gff3 \\
		  -d splitGFF3/

OPTIONS:
-g (--gff3)	GFF3 file generated by Apollo
-d (--dir)	Output directory (Optional)
OPTIONS
die "\n$usage\n" unless @ARGV;

my $gff3;
my $odir;
GetOptions(
	'g|gff3=s' => \$gff3,
	'd|dir=s' => \$odir
);

## Checking output directory
unless (defined $odir){ $odir = './'; }
unless (-d $odir){
	mkdir ($odir,0755) or die "Can't create folder $odir: $!\n";
}
print "\nOutput files will be located in directory $odir\n";

## Parsing GFF3 files
my %contigs;
my $fasta;
my ($gff, $dir) = fileparse($gff3);

print "Working on file $gff located in $dir\n\n";
open IN, "<", "$gff3" or die "Can't read $gff3: $!\n";

while (my $line = <IN>){
	chomp $line;
	if($line =~ /^##FASTA/){ ## Annotations are listed before FASTA sequences in Apollo GFF3 files
		$fasta = 1;
	}
	elsif ($line =~ /^#/){ ## Skipping comments
		next;
	}
	elsif ($line =~ /^>(\S+)/){ ## Checking for presence of FASTA sequences at the end of Apollo GFF3 files
		$fasta = $1;
	}
	elsif (($line =~ /^(\S+)/) && (!defined $fasta)){
		my $contig = $1;
		push (@{$contigs{$contig}[0]}, $line);
	}
	else {
		$contigs{$fasta}[1] .= $line;
	}
}
for my $key (keys %contigs){

	open GFF3, ">", "$odir/$key.gff3" or die "Can't create file $odir/$key.gff3: $!\n";

	if (exists $contigs{$key}[1]){ ## Prints only if FASTA is present in the GFF3 file
		open FASTA, ">", "$odir/$key.fsa" or die "Can't create file $odir/$key.fsa: $!\n";
		print FASTA ">$key\n";
		my @seq = unpack ("(A60)*", $contigs{$key}[1]);
		while (my $tmp = shift@seq){
			print FASTA "$tmp\n";
		}
		close FASTA;
	}

	if (exists $contigs{$key}[0]){
		while (my $feature = shift@{$contigs{$key}[0]}){print GFF3 "$feature\n";}
		close GFF3;
	}
}
exit;
