#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'tRNAscan_to_GFF3.pl';
my $version = '0.8c';
my $updated = '02/17/2021';

use strict; use warnings; use File::Basename; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Converts the output of tRNAscan-SE (in tabular format) to the proper GFF3 format for loading into Apollo.
USAGE		tRNAscan_to_GFF3.pl -t *.tRNAs -d tRNA/

OPTIONS:
-t (--tRNA)	tRNA file(s) to convert to GFF3
-d (--dir) 	Output directory (Optional)
OPTIONS
die "$usage\n" unless @ARGV;

my @tRNA; my $odir;
GetOptions(
	't|tRNA=s@{1,}' => \@tRNA,
	'd|dir=s' => \$odir
);

## Checking output directory
unless (defined $odir){$odir = './';}
unless (-d $odir){system "mkdir $odir";}
print "\nOutput files will be located in directory $odir\n";


my $tRNA = 0;
while (my $file = shift@tRNA){
	open TRNA, "<", "$file";
	my ($RNA, $dir) = fileparse($file);
	print "Working on file $RNA located in $dir\n\n";
	open OUT, ">", "$odir/$RNA.gff";
	open GFF3, ">", "$odir/$RNA.gff3";
	## Converting tRNAs to GFF3
	while (my $line = <TRNA>){
		chomp $line;
		if ($line =~ /^(Sequence|Name|--------)/){next;}
		else{
			my @cols = split("\t", $line);
			my $location = $cols[0]; $location =~ s/\s+$//; ## Removing padding spaces, if any
			my $tnum = $cols[1]; $tnum =~ s/\s+$//;
			my $start = $cols[2]; $start =~ s/\s+$//;
			my $end = $cols[3]; $end =~ s/\s+$//;
			my $aa = $cols[4]; $aa =~ s/\s+$//;
			my $ac = $cols[5]; $ac =~ s/\s+$//;
			my $int1 = $cols[6]; $int1 =~ s/\s+$//; ## intron boundaries, if any
			my $int2 = $cols[7]; $int2 =~ s/\s+$//;## intron boundaries, if any
			my $cove = $cols[8]; $cove =~ s/\s+$//;
			my $codon = lc($ac);
			$codon =~ tr/Tt/Uu/;
			$tRNA++;
			## Forward strand
			if ($start < $end){
				if (($int1 == 0) && ($int2 == 0)){
					print OUT "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."gene"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
				}
				else{
					my $boundary1 = $int1 - 1;
					my $boundary2 = $int2 + 1;
					print OUT "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$boundary1"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print OUT "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."gene"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$boundary1"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";

				}
			}
			elsif ($start > $end){
				if (($int1 == 0) && ($int2 == 0)){
					print OUT "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."gene"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
				}
				else{
					my $boundary1 = $int2 - 1;
					my $boundary2 = $int1 + 1;
					print OUT "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print OUT "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$boundary1"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."gene"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$location"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$boundary1"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
				}
			}
		}
	}
}
