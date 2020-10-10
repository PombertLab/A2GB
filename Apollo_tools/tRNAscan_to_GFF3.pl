#!/usr/bin/perl
## Pombert Lab, IIT, 2019
my $name = 'tRNAscan_to_GFF3.pl';
my $version = '0.7';

use strict; use warnings;

my $usage = <<"OPTIONS";
NAME		$name
VERSION		$version
SYNOPSIS	Converts the output of tRNAscan-SE (in tabular format) to the proper GFF3 format for loading into Apollo.
USAGE		tRNAscan_to_GFF3.pl *.tRNAs
OPTIONS
die "$usage\n" unless @ARGV;

my %anticodon = (
	'CGA' => 'cga','CGG' => 'cgg','CGT' => 'cgu','CGG' => 'cgg',
	'GCA' => 'gca','GCG' => 'gcg','GCT' => 'gcu','GCC' => 'gcc',
	'TCT' => 'ucu','TCC' => 'ucc','TTA' => 'uua','TTG' => 'uug',
	'CTA' => 'cua','CTG' => 'cug','ACA' => 'aca','ACG' => 'acg',
	'CTT' => 'cuu','CTC' => 'cuc','GTT' => 'guu','GTC' => 'guc',
	'CCA' => 'cca','CCG' => 'ccg','CCC' => 'ccc','CCT' => 'ccu',
	'GTA' => 'gua','GTG' => 'gug','TAA' => 'uaa','TAG' => 'uag',
	'TAT' => 'uau','AAT' => 'aau','AAC' => 'aac','GAA' => 'gaa',
	'GAG' => 'gag','GAT' => 'gau','GAC' => 'gac','TTT' => 'uuu',
	'TTC' => 'uuc','TAC' => 'uac','AAA' => 'aaa','AAG' => 'aag',
	'GGA' => 'gga','GGG' => 'ggg','GGT' => 'ggu','GGC' => 'ggc',
	'AGA' => 'aga','AGG' => 'agg','AGT' => 'agu','AGC' => 'agc',
	'TCA' => 'uca','TCG' => 'ucg','ATG' => 'end','ATT' => 'end',
	'ACT' => 'end','TGA' => 'uga','TGG' => 'ugg','TGT' => 'ugu',
	'TGC' => 'ugc','ACC' => 'acc','ATA' => 'aua','ATG' => 'aug',
	'CAA' => 'caa','CAG' => 'cag','CAT' => 'cau','CAC' => 'cac'
);

my $tRNA = 0;

while (my $file = shift@ARGV){
	open TRNA, "<$file";
	open OUT, ">$file.gff";
	open GFF3, ">$file.gff3";
	## Converting tRNAs to GFF3
	while (my $line = <TRNA>){
		chomp $line;
		if ($line =~ /^(\w+)\s+\d+\s+(\d+)\s+(\d+)\s+(\w{3,6})\t(\w{3})\t(\d+)\t(\d+)\t(\d+.\d+)/){
			my $target = $1;
			my $start = $2;
			my $end = $3;
			my $aa = $4;
			my $ac = $5;
			my $int1 = $6; ## intron boundaries, if any
			my $int2 = $7; ## intron boundaries, if any
			my $cove = $8;
			my $codon = $ac;
			$codon =~ tr/Tt/Uu/;
			$tRNA++;
			## Forward strand
			if ($start < $end){
				if (($int1 == 0) && ($int2 == 0)){
					print OUT "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."gene"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
				}
				else{
					my $boundary1 = $int1 - 1;
					my $boundary2 = $int2 + 1;
					print OUT "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$boundary1"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print OUT "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."gene"."\t"."$start"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$start"."\t"."$boundary1"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$end"."\t"."$cove"."\t".'+'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";

				}
			}
			elsif ($start > $end){
				if (($int1 == 0) && ($int2 == 0)){
					print OUT "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."gene"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
				}
				else{
					my $boundary1 = $int2 - 1;
					my $boundary2 = $int1 + 1;
					print OUT "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print OUT "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$boundary1"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."gene"."\t"."$end"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$boundary2"."\t"."$start"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
					print GFF3 "$target"."\t"."TRNASCAN"."\t"."tRNA"."\t"."$end"."\t"."$boundary1"."\t"."$cove"."\t".'-'."\t".'.'."\t"."ID=tRNA$tRNA".';'."Name=tRNA$tRNA".';'."Note=tRNA-$aa($codon)"."\n";
				}
			}
		}
	}
}
