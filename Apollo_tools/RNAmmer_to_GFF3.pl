#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'RNAmmer_to_GFF3.pl';
my $version = '0.7';
my $updated = '2021-04-09';

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Converts the output of RNAmmer GFF2 output to the proper GFF3 format for loading into Apollo

USAGE		${name} \\
		  -g *.gff2 \\
		  -d rRNA/

OPTIONS:
-g (--gff2)	RNAmmer gff2 input file
-d (--dir) 	Output directory (Optional)
OPTIONS
die "\n$usage\n" unless @ARGV;

my @gff2;
my $odir = './';
GetOptions(
	'g|gff2=s@{1,}' => \@gff2,
	'd|dir=s' => \$odir
);

## Checking output directory
unless (-d $odir){
	mkdir ($odir,0755) or die "Can't create folder $odir: $!\n";
}
print "Output files will be located in directory $odir\n";

## Setting up variables
my $location;
my $source;
my $feature;
my $start;
my $end;
my $score;
my $strand;
my $gene;
my $rRNA = 0;

## Workign on gff file
while (my $file = shift@gff2){

	open RRNA, "<", "$file" or die "Can't read $file: $!\n";
	my ($gff2, $dir) = fileparse($file);

	print "Working on file $gff2 located in $dir\n";
	$gff2 =~ s/.gff2$//;
	open GFF, ">", "$odir/$gff2.gff" or die "Can't create $odir/$gff2.gff: $!\n";
	open GFF3, ">", "$odir/$gff2.gff3" or die "Can't create $odir/$gff2.gff3: $!\n";
		
	## Converting rRNAs to GFF3
	while (my $line = <RRNA>){
		chomp $line;
		if ($line =~ /^#/){ next; } ## disregard comments
		else{
			my @cols = split("\t", $line);
			$location = $cols[0];
			$source = $cols[1];
			$feature = $cols[2];
			$start = $cols[3];
			$end = $cols[4];
			$score = $cols[5];
			$strand = $cols[6];
			$gene = $cols[8];
			$rRNA++;

			data_print('rRNA', $strand, \*GFF);
			data_print('gene', $strand, \*GFF3);
			data_print('rRNA', $strand, \*GFF3);
		}
	}
}

## Subroutines
sub data_print {
	my ($type, $strand, $fh) = @_;
	print $fh "$location"."\t"."RNAMMER"."\t"."$type"."\t"."$start"."\t"."$end"."\t"."$score"."\t";
	print $fh "$strand"."\t".'.'."\t"."ID=rRNA$rRNA".';'."Name=rRNA$rRNA".';'."Note=$gene"."\n";
}