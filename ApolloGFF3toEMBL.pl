#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'ApolloGFF3toEMBL.pl';
my $version = '3.9';

use strict; use warnings; use File::Basename; use Bio::SeqIO; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";

NAME		$name
VERSION		$version
SYNOPSIS	Converts Apollo GFF3 files to EMBL files and writes the proteins and RNAs to separate FASTA files with the .prot and .RNA extensions.
		Generates locus tags automatically based on the provided prefix.
REQUIREMENTS	BioPerl's Bio::SeqIO module
NOTE		The GFF3 (*.gff3) and corresponding FASTA (*.fsa) files must be in the same folder

USAGE		$name -p LOCUS_TAG_PREFIX -g *.gff3 -c 1

OPTIONS:
-p (--prefix)	## locus_tag prefix
-g (--gff)	## GFF3 files generated by Apollo
-c (--gcode)	## NCBI genetic code [Default: 1]
		1  - The Standard Code
		2  - The Vertebrate Mitochondrial Code
		3  - The Yeast Mitochondrial Code
		4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		11 - The Bacterial, Archaeal and Plant Plastid Code
		NOTE - For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
OPTIONS
die "$usage\n" unless @ARGV;

my $locus_tag_prefix;
my @gff3;
my $gc = 1;
GetOptions(
	'p|prefix=s' => \$locus_tag_prefix,
	'g|gff=s@{1,}' => \@gff3,
	'c|gcode=i' => \$gc
);

my $locus_id = 10; my $contig_number = 0;
my $protein = undef; my $mRNA =undef;
my %gcodes; gcodes();

my $gff; my $dir;
while (my $file = shift@gff3){
	($gff, $dir) = fileparse($file);
	print "Working on file $gff located in $dir\n";
	open IN, "<", "$file";
	$file =~ s/.gff3$//;
	open DNA, "<", "$file.fsa";
	open OUT, ">", "$file.embl";
	open PROT, ">", "$file.prot";
	open MRNA, ">", "$file.RNA";
	if ($file =~ /(\d+)$/){$contig_number = $1;}
	
	### Creating a single DNA string for protein translation
	my $DNAseq= undef;
	while (my $dna = <DNA>){
		chomp $dna;
		if ($dna =~ /^>/){next;}
		else{$DNAseq.= $dna;}
	}
	my $DNAsequence = uc($DNAseq); 			## Changing to upper case to fit with the translation hash
	my $contig_length = length($DNAsequence);	## Calculating the contig size
	
	### Init hashes, arrays and values
	my $geneid = 0;		## Initialize gene id
	my $CDS_counter = 0;	## Initialize CDS counter
	my %gene = ();		## Initialize gene hash
	my %exon = ();		## Initialize exon hash
	my %strands = ();	## Initialize strandedness hash
	my @todo = ();		## Initialize array of things to do
	
	my %feature;
	### Filling the hashes and arrays ###
	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^\S+\t\S+\tgene\t(\d+)\t(\d+)\t\S+\t([-+])/){
			my $start = $1;
			my $end = $2;
			my $strandedness = $3;
			$geneid = $start;
			push (@{$gene{$geneid}}, $start, $end);
			push (@todo, $geneid);
			$strands{$geneid} = $strandedness;
			$CDS_counter = 0;
		}
		elsif ($line =~ /^\S+\t\S+\t(CDS|rRNA|tRNA|tmRNA|ncRNA|snRNA|snoRNA)\t(\d+)\t(\d+)\t\S+\t([-+])\t(\d|\.)/){
			my $type = $1;
			my $start = $2;
			my $end = $3;
			my $strandedness = $4;
			my $phase = $5; ## 0,1,2 (CDS) or . for RNAs
			$feature{$geneid} = $type;
			if ($strandedness eq '+'){
				push (@{$exon{$geneid}}, $start, $end);
				$CDS_counter++;
			}
			elsif ($strandedness eq '-'){
				push (@{$exon{$geneid}}, $end, $start);
				$CDS_counter++;
			}
		}
	}
	
	### Working on the gene list. Using the @todo array so that we don't have to sort the hash ouput ###
	my @sorted_todo = sort {$a <=> $b} @todo;
	while (my $list = shift@sorted_todo){
		if ($strands{$list} eq '+'){ ## Looking for positive strandedness
			
			### Working on gene features
			my $locus_number = sprintf("%05d", $locus_id);
			my $contig = sprintf("%02d", $contig_number);
			print OUT 'FT   gene             '."$gene{$list}[0]".'..'."$gene{$list}[1]"."\n";
			print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'g'."$locus_number".'"'."\n";
			$locus_id += 10;
				
			### Working on CDS features
			my $cds_count = scalar(@{$exon{$list}});
			my $end = ($cds_count - 1);
			my $num = ($cds_count - 2);
			my $stopcodon = $exon{$list}[$end];
			
			if (scalar(@{$exon{$list}}) == 2){ ## Verifying if we have a single exon
				print OUT 'FT   '."$feature{$list}".'             '."$exon{$list}[0]..$stopcodon\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'g'."$locus_number".'"'."\n";
				if ($feature{$list} eq 'CDS'){print PROT ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";}
				print MRNA ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";
				my $start = (($exon{$list}[0]));
				my $stop = (($exon{$list}[1]));
				$protein = undef;
				$mRNA =  substr($DNAsequence, $start-1, (($stop-$start)+1));
				rna();
				if ($feature{$list} eq 'CDS'){
					for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
						my $codon = substr($mRNA, $i, 3);
						$protein .=  $gcodes{$gc}{$codon};
					}
					prot();
				}
			}
			elsif (scalar(@{$exon{$list}}) > 2){ ## Verifying if we have more than one exon
				print OUT 'FT   '."$feature{$list}".'             '."join($exon{$list}[0]..";
				if ($feature{$list} eq 'CDS'){print PROT ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";}
				print MRNA ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";
				foreach my $count (1..$num){
					if ($count % 2){print OUT "$exon{$list}[$count],";} # Working on odd numbers
					else{print OUT "$exon{$list}[$count]..";} # Working on even numbers
				}
				print OUT "$stopcodon)\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'g'."$locus_number".'"'."\n";
				
				$mRNA = undef;
				$protein = undef;
				my $tmp1 = undef;
				my $tmp2 = undef;
				foreach my $subs (0..$end){
					if ($subs % 2){# Working on odd numbers
						$tmp2 = ($exon{$list}[$subs]-1);
						$mRNA .= substr($DNAsequence, $tmp1, (($tmp2-$tmp1)+1));
					}
					else{	# Working on even numbers
						$tmp1 = ($exon{$list}[$subs]-1);
					}
				}
				rna();
				if ($feature{$list} eq 'CDS'){
					for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
						my $codon = substr($mRNA, $i, 3);
						$protein .= $gcodes{$gc}{$codon};
					}
					prot();
				}
			}
		}
		if ($strands{$list} eq '-'){ ## Looking for negative strandedness
			
			### Working on gene features
			my $locus_number = sprintf("%05d", $locus_id);
			my $contig = sprintf("%02d", $contig_number);
			print OUT 'FT   gene             complement('."$gene{$list}[0]".'..'."$gene{$list}[1]".')'."\n";
			print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'g'."$locus_number".'"'."\n";
			$locus_id += 10;
				
			### Working on CDS features
			my $cds_count = scalar(@{$exon{$list}});
			my $end = ($cds_count - 1);
			my $num = ($cds_count - 2);
			my @reversed = reverse(@{$exon{$list}}); ## Reversing the list of exons
			my $stopcodon = $reversed[$end];
			
			if (scalar(@reversed) == 2){ ## Verifying if we have a single exon
				print OUT 'FT   '."$feature{$list}".'             complement('."$reversed[0]..$reversed[1]".')'."\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'g'."$locus_number".'"'."\n";
				if ($feature{$list} eq 'CDS'){print PROT ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";}
				print MRNA ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";
				my $start = (($reversed[0]));
				my $stop = (($reversed[1]));
				$protein = undef;
				my $revmRNA =  substr($DNAsequence, $start-1, (($stop-$start)+1));
				$mRNA = reverse($revmRNA);
				$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
				rna();
				if ($feature{$list} eq 'CDS'){
					for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
						my $codon = substr($mRNA, $i, 3);
						$protein .= $gcodes{$gc}{$codon};
					}
					prot();
				}
			}
			elsif (scalar(@reversed) > 2){ ## Verifying if we have more than one exon
				print OUT 'FT   '."$feature{$list}".'             complement('."join($reversed[0]..";
				if ($feature{$list} eq 'CDS'){print PROT ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";}
				print MRNA ">$locus_tag_prefix".'_'."$contig".'g'."$locus_number\n";
				foreach my $count (1..$num){
					if ($count % 2){print OUT "$reversed[$count],";} # Working on odd numbers
					else{print OUT "$reversed[$count]..";} # Working on even numbers
				}
				print OUT "$stopcodon)\n";
				print OUT 'FT                   /locus_tag="'."$locus_tag_prefix".'_'."$contig".'g'."$locus_number".'"'."\n";
	
				my $revmRNA = undef;
				$protein = undef;
				my $tmp1 = undef;
				my $tmp2 = undef;
				foreach my $subs (0..$end){
					if ($subs % 2){# Working on odd numbers
						$tmp2 = ($reversed[$subs]-1);
						$revmRNA .= substr($DNAsequence, $tmp1, (($tmp2-$tmp1)+1));
					}
					else{	# Working on even numbers
						$tmp1 = ($reversed[$subs]-1);
					}
				}
				my $mRNA = reverse($revmRNA);
				$mRNA =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
				rna();
				if ($feature{$list} eq 'CDS'){
					for(my $i = 0; $i < (length($mRNA) - 5); $i += 3){ ## -2 if stop codon (*) desired, -5 if not
						my $codon = substr($mRNA, $i, 3);
						$protein .= $gcodes{$gc}{$codon};
					}
					prot();
				}
			}
		}
	}
	### Converting Fasta input to EMBL with BioPerl SeqIO
	my $in = Bio::SeqIO->new(-file => "$file.fsa", -format => "fasta");
	my $out = Bio::SeqIO->new(-file => ">$file.tmp1", -format => "embl");
	while (my $seq = $in->next_seq()){
		$out->write_seq($seq);
	}
	### Adding formatted sequence to the EMBL file
	open SEQ, "<$file.tmp1";
		while (my $read = <SEQ>){
			if ($read =~ /^AC|ID|XX|DE|FH/){next;}
			else{print OUT "$read"};
		}
	close SEQ;
}
close IN; close OUT;
system "rm ${dir}/*.tmp1";

## Subroutines
sub prot{my @PROTEIN = unpack ("(A60)*", $protein); while (my $tmp = shift@PROTEIN){print PROT "$tmp\n";}}
sub rna{my @RNA = unpack ("(A60)*", $mRNA);	while (my $tmp = shift@RNA){print MRNA "$tmp\n";}}
sub gcodes{ ## NCBI Genetic codes
	%gcodes = (
		1 => { ## The Standard Code (transl_table=1)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',   
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		2 => { ## The Vertebrate Mitochondrial Code (transl_table=2)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => '*',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => '*',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		3 => { ## The Yeast Mitochondrial Code (transl_table=3)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'T', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'T', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'T', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'T', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		4 => { ## The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code (transl_table=4)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G' 
		},
		5 => { ## The Invertebrate Mitochondrial Code (transl_table=5)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'S',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		6 => { ## The Ciliate, Dasycladacean and Hexamita Nuclear Code (transl_table=6)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C', 
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C', 
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => '*', 
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W', 
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R', 
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R', 
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R', 
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R', 
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S', 
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S', 
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R', 
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R', 
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G', 
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G', 
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G', 
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G' 
		},
		9 => { ## The Echinoderm and Flatworm Mitochondrial Code (transl_table=9)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		10 => { ## The Euplotid Nuclear Code (transl_table=10)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'C',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		11 => { ## The Bacterial, Archaeal and Plant Plastid Code (transl_table=11)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		12 => { ## The Alternative Yeast Nuclear Code (transl_table=12)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'S', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		13 => { ## The Ascidian Mitochondrial Code (transl_table=13)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'G',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'G',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
		},
		14 => { ## The Alternative Flatworm Mitochondrial Code (transl_table=14)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Y', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G',
		},
		16 => { ## Chlorophycean Mitochondrial Code (transl_table=16)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'L', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		21 => { ## Trematode Mitochondrial Code (transl_table=21)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'M', 'ACA' => 'T', 'AAA' => 'N', 'AGA' => 'S',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'S',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		22 => { ## Scenedesmus obliquus Mitochondrial Code (transl_table=22)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',  
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',  
			'TTA' => 'L', 'TCA' => '*', 'TAA' => '*', 'TGA' => '*',  
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'L', 'TGG' => 'W',  
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',  
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',  
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',  
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',  
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',  
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',  
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',  
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',  
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',  
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',  
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',  
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'  
		},
		23 => { ## Thraustochytrium Mitochondrial Code (transl_table=23)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => '*', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		24 => { ## Pterobranchia Mitochondrial Code (transl_table=24)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'S',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'K',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		25 => { ## Candidate Division SR1 and Gracilibacteria Code (transl_table=25)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => 'G',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		26 => { ## Pachysolen tannophilus Nuclear Code (transl_table=26)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => '*', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => '*', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'A', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		27 => { ## Karyorelict Nuclear (transl_table=27)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		28 => { ## Condylostoma Nuclear (transl_table=28)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Q', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Q', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		29 => { ## Mesodinium Nuclear (transl_table=29)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'Y', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'Y', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		30 => { ## Peritrich Nuclear (transl_table=30)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'E', 'TGA' => '*',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'E', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
		31 => { ## Blastocrithidia Nuclear (transl_table=31)
			'TTT' => 'F', 'TCT' => 'S', 'TAT' => 'Y', 'TGT' => 'C',
			'TTC' => 'F', 'TCC' => 'S', 'TAC' => 'Y', 'TGC' => 'C',
			'TTA' => 'L', 'TCA' => 'S', 'TAA' => 'E', 'TGA' => 'W',
			'TTG' => 'L', 'TCG' => 'S', 'TAG' => 'E', 'TGG' => 'W',
			'CTT' => 'L', 'CCT' => 'P', 'CAT' => 'H', 'CGT' => 'R',
			'CTC' => 'L', 'CCC' => 'P', 'CAC' => 'H', 'CGC' => 'R',
			'CTA' => 'L', 'CCA' => 'P', 'CAA' => 'Q', 'CGA' => 'R',
			'CTG' => 'L', 'CCG' => 'P', 'CAG' => 'Q', 'CGG' => 'R',
			'ATT' => 'I', 'ACT' => 'T', 'AAT' => 'N', 'AGT' => 'S',
			'ATC' => 'I', 'ACC' => 'T', 'AAC' => 'N', 'AGC' => 'S',
			'ATA' => 'I', 'ACA' => 'T', 'AAA' => 'K', 'AGA' => 'R',
			'ATG' => 'M', 'ACG' => 'T', 'AAG' => 'K', 'AGG' => 'R',
			'GTT' => 'V', 'GCT' => 'A', 'GAT' => 'D', 'GGT' => 'G',
			'GTC' => 'V', 'GCC' => 'A', 'GAC' => 'D', 'GGC' => 'G',
			'GTA' => 'V', 'GCA' => 'A', 'GAA' => 'E', 'GGA' => 'G',
			'GTG' => 'V', 'GCG' => 'A', 'GAG' => 'E', 'GGG' => 'G'
		},
	);
}
