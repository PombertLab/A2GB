#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'EMBLtoTBL.pl';
my $version = '1.5b';
my $updated = '02/24/21';

use strict; use warnings; use Bio::SeqIO; use File::Basename; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Converts EMBL files to NCBI TBL format for TBL2ASN
REQUIREMENTS	BioPerl's Bio::SeqIO module
NOTE		The EMBL (*.embl) and FASTA (*.fsa) files must be in the same folder.
		Requires locus_tags to be defined in the EMBL files.
		
USAGE		${name} -id IITBIO -p product_list.txt -embl *.embl
OPTIONS:
-id		Desired institute ID [default: IITBIO]
-p		Tab-delimited list of locus_tags and their products
-embl		EMBL files to convert
OPTIONS
die "$usage\n" unless @ARGV;

my $instID = 'IITBIO'; ## 
my $products; ## protein_list.txt
my @embl;
GetOptions(
	'id=s' => \$instID,
	'p=s' => \$products,
	'embl=s@{1,}' => \@embl,
);


sub numSort {if ($a < $b){return -1;} elsif ($a == $b){return 0;} elsif ($a > $b) {return 1;}}

### Filling the products database
my %hash = ();
open HASH, "<$products";
while(my $dbkey = <HASH>){chomp $dbkey;if($dbkey =~ /^(\S+)\t(.*)$/){my $prot = $1;my $prod = $2;$hash{$prot}=$prod;}}

### Working on EMBL files
my $locus_tag;
while(my $file = shift@embl){
	open IN, "<", "$file" or die "Can't open EMBL file file: $file\n";
	$file =~ s/.embl$//;
	my ($head, $dir) = fileparse($file);
	open DNA, "<", "$file.fsa" or die "Can't open FASTA file file: $file.fsa\n";
	open TBL, ">", "$file.tbl" or die "Can't create TBL output file: $file.tbl\n";
	print TBL ">Feature $head\n"; ## Generate TBL header
	
	### Creating a single DNA string for codon verification
	my $DNAseq = undef;
	while (my $dna = <DNA>){chomp $dna;if ($dna =~ /^>/){next;}else{$DNAseq.= $dna;}}
	my $DNAsequence = lc($DNAseq); ## Changing to lower case to fit with the codon check
	my $contig_length = length($DNAsequence); ## Calculating the contig size
	$locus_tag = undef;

	while(my $line = <IN>){
		chomp $line;
		my @start = ();
		my @stop = ();
		my $asize = undef;
		my $num = undef;
		my $dum = undef;
		
		### Defining the locus tags
		if ($line =~ /FT\s+\/locus_tag="(\S+)"/){$locus_tag = $1;}
		
		### Working on tRNAs/rRNAs
		elsif ($line =~ /FT\s+(tRNA|rRNA)\s+(\d+)..(\d+)/){ ## Forward, single exon
			my $type = $1;
			my $start = $2;
			my $stop = $3;
			print TBL "$start\t$stop\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			print TBL "$start\t$stop\t$type\n";
			RNA();
		}
		elsif ($line =~ /FT\s+(tRNA|rRNA)\s+join\((.*)\)/){ ## Forward, multiple exons
			my $feat = $1;
			my @array = split(',',$2);
			my $mRNA = undef;
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $1;
					my $stp = $2;
					push (@start, $strt);
					push (@stop, $stp);
				}
			}
			$asize = scalar(@start);
			$num = ($asize -1);
			$dum = ($asize -2);

			### Printing gene info
			print TBL "$start[0]\t$stop[$num]\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			
			### Printing tRNA/rRNA info
			print TBL "$start[0]\t$stop[0]\t$feat\n"; ## printing the 1st exon
			if ($asize == 2){print TBL "$start[$num]\t$stop[$num]\n";} ## printing the last exon
			elsif ($asize >= 3){
				foreach my $subs (1..$dum){print TBL "$start[$subs]\t$stop[$subs]\n";}
				print TBL "$start[$num]\t$stop[$num]\n"; ## printing the last exon	
			}
			RNA();
		}
		elsif ($line =~ /FT\s+(tRNA|rRNA)\s+complement\((\d+)..(\d+)\)/){ ## Reverse, single exon
			my $type = $1;
			my $start = $3;
			my $stop = $2;
			print TBL "$start\t$stop\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			print TBL "$start\t$stop\t$type\n";
			RNA();
		}
		elsif ($line =~ /FT\s+(tRNA|rRNA)\s+complement\(join\((.*)/){ ## Reverse, mutiple exons
			my $feat = $1;
			my @array = split(',',$2);
			my $mRNA = undef;
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $2;
					my $stp = $1;
					unshift (@start, $strt);
					unshift (@stop, $stp);
				}
			}
			$asize = scalar(@start);
			$num = ($asize -1);
			$dum = ($asize -2);
			my @sstart = sort numSort @start;
			my @sstop = sort numSort @stop;
			
			### Printing gene info
			print TBL "$start[0]\t$stop[$num]\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			
			### Printing tRNA/rRNA info
			print TBL "$start[0]\t$stop[0]\t$feat\n"; ## printing the 1st exon
			if ($asize == 2){print TBL "$start[$num]\t$stop[$num]\n";}
			elsif ($asize >= 3){
				foreach my $subs (1..$dum){print TBL "$start[$subs]\t$stop[$subs]\n";}
				print TBL "$start[$num]\t$stop[$num]\n"; ## printing the last exon	
			}
			RNA();
		}
		
		### Working on CDS
		elsif ($line =~ /FT\s+\/codon_start=(\d)/){print TBL "\t\t\tcodon_start\t$1\n";} ## Looking for phased codons
		elsif ($line =~ /FT\s+CDS\s+(\d+)..(\d+)/){ ## Forward, single exon
			my $start = $1;
			my $stop = $2;
			my $startcodon = substr($DNAsequence, $start-1, 3);
			my $stopcodon = substr($DNAsequence, $stop-3, 3);
			print TBL "<$start\t>$stop\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			## Defining mRNA
			print TBL "<$start\t>$stop\tmRNA\n";
			mRNA();
			if (($startcodon eq 'atg') && (($stopcodon eq 'taa') || ($stopcodon eq'tag') || ($stopcodon eq'tga'))){ ## =, =
				# print "$file\t$start\t$stop\t$contig_length\n"; Debugging line
				if (($start == 1) && ($stop == $contig_length)){print TBL "<$start\t>$stop\tCDS\n";}
				elsif (($start == 1) && ($stop < $contig_length)){print TBL "<$start\t$stop\tCDS\n";}
				elsif (($start > 1) && ($stop == $contig_length)){print TBL "$start\t>$stop\tCDS\n";}
				else{print TBL "$start\t$stop\tCDS\n";}
			}
			elsif (($startcodon ne 'atg') && (($stopcodon eq 'taa') || ($stopcodon eq'tag') || ($stopcodon eq'tga'))){## !=, =
				if ($stop == $contig_length){print TBL "<$start\t>$stop\tCDS\n";}
				else{print TBL "<$start\t$stop\tCDS\n";} 
			}
			elsif ($startcodon eq 'atg'){ ## =, !=
				if ($start == 1){print TBL "<$start\t>$stop\tCDS\n";}
				else{print TBL "$start\t>$stop\tCDS\n";}
			}
			else{print TBL "<$start\t>$stop\tCDS\n";}## !=, !=
			CDS();
		}
		elsif ($line =~ /FT \s+CDS\s+complement\((\d+)..(\d+)\)/){ ## Reverse, single exon
			my $start = $2;
			my $stop = $1;
			my $startcodon = substr($DNAsequence, $start-3, 3);
			my $stopcodon = substr($DNAsequence, $stop-1, 3);
			print TBL "<$start\t>$stop\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			## defining mRNA
			print TBL "<$start\t>$stop\tmRNA\n";
			mRNA();
			if (($startcodon eq 'cat') && (($stopcodon eq 'tta') || ($stopcodon eq'cta') || ($stopcodon eq'tca'))){ ## =, =
				if (($start == $contig_length) && ($stop == 1)){print TBL "<$start\t>$stop\tCDS\n";}
				elsif (($start == $contig_length) && ($stop > 1)){print TBL "<$start\t$stop\tCDS\n";}
				elsif (($start < $contig_length) && ($stop == 1)){print TBL "$start\t>$stop\tCDS\n";}
				else{print TBL "$start\t$stop\tCDS\n";}
			}
			elsif (($startcodon ne 'cat') && (($stopcodon eq 'tta') || ($stopcodon eq'cta') || ($stopcodon eq'tca'))){## !=, =
				if ($stop == 1){print TBL "<$start\t>$stop\tCDS\n";}
				else{print TBL "<$start\t$stop\tCDS\n";}
			}
			elsif ($startcodon eq 'cat'){ ## =, !=
				if ($start == $contig_length){print TBL "<$start\t>$stop\tCDS\n";}
				else{print TBL "$start\t>$stop\tCDS\n";}
			}
			else{print TBL "<$start\t>$stop\tCDS\n";} ## !=, !=
			CDS();
		}	
		elsif ($line =~ /FT\s+CDS\s+join\((.*)\)/){ ## Forward, multiple exons
			my @array = split(',',$1);
			my $mRNA = undef;
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $1;
					my $stp = $2;
					push (@start, $strt);
					push (@stop, $stp);
				}
			}
			$asize = scalar(@start);
			$num = ($asize -1);
			$dum = ($asize -2);

			### Printing gene info
			print TBL "<$start[0]\t>$stop[$num]\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			
			### Printing mRNA info
			if ($asize == 1){print TBL "<$start[0]\t>$stop[0]\tmRNA\n";}
			elsif ($asize == 2){print TBL "<$start[0]\t$stop[0]\tmRNA\n"; print TBL "$start[1]\t>$stop[1]\n";}
			elsif ($asize >= 3){
				print TBL "<$start[0]\t$stop[0]\tmRNA\n";
				foreach my $subs (1..$dum){print TBL "$start[$subs]\t$stop[$subs]\n";}
				print TBL "$start[$num]\t>$stop[$num]\n";	
			}
			mRNA();
			
			### Printing CDS info
			my $startcodon = substr($DNAsequence, $start[0]-1, 3);
			my $stopcodon = substr($DNAsequence, $stop[$num]-3, 3);
			if ($startcodon eq 'atg'){
				if ($start[0] == 1){print TBL "<$start[0]\t$stop[0]\tCDS\n";}
				else{print TBL "$start[0]\t$stop[0]\tCDS\n";}
			} ## printing the 1st exon
			else{print TBL "<$start[0]\t$stop[0]\tCDS\n";} ## printing the 1st exon
			if ($asize == 2){ ## Do we have only two exons?
				if (($stopcodon eq 'taa') || ($stopcodon eq'tag') || ($stopcodon eq'tga')){
					if ($stop[$num] == $contig_length){print TBL "$start[$num]\t>$stop[$num]\n";}
					else{print TBL "$start[$num]\t$stop[$num]\n";} ## printing the last exon
				}
				else{print TBL "$start[$num]\t>$stop[$num]\n";} ## printing the last exon	
			}
			elsif ($asize >= 3){
				foreach my $subs (1..$dum){print TBL "$start[$subs]\t$stop[$subs]\n";}
				if (($stopcodon eq 'taa') || ($stopcodon eq'tag') || ($stopcodon eq'tga')){
					if ($stop[$num] == $contig_length){print TBL "$start[$num]\t>$stop[$num]\n";}
					else{print TBL "$start[$num]\t$stop[$num]\n";}
				} ## printing the last exon
				else{print TBL "$start[$num]\t>$stop[$num]\n";} ## printing the last exon	
			}
			CDS();
		}
		elsif ($line =~ /FT\s+CDS\s+complement\(join\((.*)/){ ## Reverse, mutiple exons
			my @array = split(',',$1);
			my $mRNA = undef;
			while (my $segment = shift@array){
				chomp $segment;
				if ($segment =~ /(\d+)..(\d+)/){
					my $strt = $2;
					my $stp = $1;
					unshift (@start, $strt);
					unshift (@stop, $stp);
				}
			}
			$asize = scalar(@start);
			$num = ($asize -1);
			$dum = ($asize -2);
			my @sstart = sort numSort @start;
			my @sstop = sort numSort @stop;
			
			### Printing gene info
			print TBL "<$start[0]\t>$stop[$num]\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			
			### Printing mRNA info
			if ($asize == 1){print TBL "<$start[0]\t>$stop[0]\tmRNA\n";}
			elsif ($asize == 2){print TBL "<$start[0]\t$stop[0]\tmRNA\n"; print TBL "$start[1]\t>$stop[1]\n";}
			elsif ($asize >= 3){
				print TBL "<$start[0]\t$stop[0]\tmRNA\n";
				foreach my $subs (1..$dum){print TBL "$start[$subs]\t$stop[$subs]\n";}
				print TBL "$start[$num]\t>$stop[$num]\n";	
			}
			mRNA();

			### Printing CDS info
			my $startcodon = substr($DNAsequence, $start[0]-3, 3);
			my $stopcodon = substr($DNAsequence, $stop[$num]-1, 3);
			if ($startcodon eq 'cat'){
				if ($start[0] == $contig_length){print TBL "<$start[0]\t$stop[0]\tCDS\n";}
				else{print TBL "$start[0]\t$stop[0]\tCDS\n";}
			} ## printing the 1st exon
			else{print TBL "<$start[0]\t$stop[0]\tCDS\n";} ## printing the 1st exon
			if ($asize == 2){ ## Do we have only two exons?
				if (($stopcodon eq 'tta') || ($stopcodon eq'cta') || ($stopcodon eq'tca')){
					if ($stop[$num] == 1){print TBL "$start[$num]\t>$stop[$num]\n";}
					else{print TBL "$start[$num]\t$stop[$num]\n";}
				} ## printing the last exon
				else{print TBL "$start[$num]\t>$stop[$num]\n";} ## printing the last exon	
			}
			elsif ($asize >= 3){
				foreach my $subs (1..$dum){print TBL "$start[$subs]\t$stop[$subs]\n";}
				if (($stopcodon eq 'tta') || ($stopcodon eq'cta') || ($stopcodon eq'tca')){
					if ($stop[$num] == 1){print TBL "$start[$num]\t>$stop[$num]\n";}
					else{print TBL "$start[$num]\t$stop[$num]\n";}
				} ## printing the last exon
				else{print TBL "$start[$num]\t>$stop[$num]\n";} ## printing the last exon	
			}
			CDS();
		}
	}
}
close IN;
close TBL;

### Subroutines
sub mRNA{
	print TBL "\t\t\tlocus_tag\t$locus_tag\n";
	if (exists $hash{$locus_tag}){ print TBL "\t\t\tproduct\t$hash{$locus_tag}\n"; }
	else{ print TBL "\t\t\tproduct\thypothetical protein\n"; }
	print TBL "\t\t\tprotein_id\tgnl|$instID|$locus_tag\n";
	print TBL "\t\t\ttranscript_id\tgnl|$instID|$locus_tag"."_mRNA\n";
}
sub CDS{
	print TBL "\t\t\tlocus_tag\t$locus_tag\n";
	if (exists $hash{$locus_tag}){ print TBL "\t\t\tproduct\t$hash{$locus_tag}\n"; }
	else{
		print STDERR "Cannot find database entry for locus_tag: $locus_tag\n";
		print TBL "\t\t\tproduct\thypothetical protein\n";
	}
	print TBL "\t\t\tprotein_id\tgnl|$instID|$locus_tag\n";
	print TBL "\t\t\ttranscript_id\tgnl|$instID|$locus_tag"."_mRNA\n";
}
sub RNA{
	print TBL "\t\t\tlocus_tag\t$locus_tag\n";
	if (exists $hash{$locus_tag}){ print TBL "\t\t\tproduct\t$hash{$locus_tag}\n"; }
	else{
		print STDERR "Cannot find database entry for locus_tag: $locus_tag\n";
		print TBL "\t\t\tproduct\thypothetical RNA\n";
	}
}
