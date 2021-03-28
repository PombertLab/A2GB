#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'EMBLtoFeatures.pl';
my $version = '2.0';
my $updated = '03/28/21';

use strict; use warnings; use File::Basename; use Bio::SeqIO; use Getopt::Long qw(GetOptions);

### Defining options
my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Writes features to separate FASTA files with the extensions .prot, .RNA, exons, etc.

REQUIREMENTS	BioPerl's Bio::SeqIO module
NOTE		The EMBL (*.embl) and corresponding FASTA (*.fsa) files must be in the same folder

USAGE		${name} \\
		  -e *.embl \\
		  -o FILES \\
		  -x \\
		  -i \\
		  -v \\
		  -c 1

OPTIONS:
-e (--embl)	EMBL files
-o (--outdir)	Output directory [Default: ./]
-x (--exon)	Create fasta files of exons (.exons) fo genes with introns
-i (--intron)	Create fasta files of exons (.introns) of genes with introns
-v (--verbose)	Add Verbosity
-c (--gcode)	NCBI genetic code [Default: 1]
		1  - The Standard Code
		2  - The Vertebrate Mitochondrial Code
		3  - The Yeast Mitochondrial Code
		4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		11 - The Bacterial, Archaeal and Plant Plastid Code
		NOTE - For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
OPTIONS
die "\n$usage\n" unless @ARGV;

my @embl;
my $outdir = './';
my $exon_feature;
my $intron_feature;
my $gc = 1;
my $verb;
GetOptions(
	'e|embl=s@{1,}' => \@embl,
	'o|outdir=s' => \$outdir,
	'x|exon' => \$exon_feature,
	'i|intron' => \$intron_feature,
	'c|gcode=i' => \$gc,
	'v|verbose' => \$verb
);

### Creating output directory
unless (-d $outdir){
	mkdir ($outdir,0755) or die "Can't create folder $outdir: $!\n";
}
if ($verb){ print "\nSetting output directory to: $outdir\n\n";}

### Loading NCBI genetic codes
my %gcodes; gcodes();

### Parsing EMBL files
sub numSort {
		if ($a < $b) { return -1; }
		elsif ($a == $b) { return 0;}
		elsif ($a > $b) { return 1; }
}

my $locus_tag;
my $locus;
my $protein;
my $mRNA;
my %codon_start;
my $exon_counter; ## To keep track of exon numbers
my $intron_counter; ## To keep track of intron numbers

while (my $file = shift@embl){
	
	open IN, "<", "$file" or die "Can't open $file: $!\n";
	my ($embl, $dir) = fileparse($file);
	if($verb) { print "Working on file $embl located in $dir\n"; }
	
	$file =~ s/.embl$//;

	open DNA, "<", "$file.fsa" or die "Can't open $file.fsa: $!\n";

	## Creating output files
	my $prot_file = "${outdir}/$file.prot";
	my $RNA_file = "${outdir}/$file.RNA";
	my $exon_file = "${outdir}/$file.exons";
	my $intron_file = "${outdir}/$file.introns";

	open PROT, ">", "$prot_file" or die "Can't create $prot_file: $!\n";
	open RNA, ">", "$RNA_file" or die "Can't create $RNA_file: $!\n";
	if ($exon_feature){
		open EXON, ">", "$exon_file" or die "Can't create $exon_file: $!\n";
	}
	if ($intron_feature){
		open INTRON, ">", "$intron_file" or die "Can't create $intron_file: $!\n";
	}

	### Creating a single DNA string for protein translation
	my $DNAseq;
	while (my $dna = <DNA>){
		chomp $dna;
		if ($dna =~ /^>/){ next; }
		else{$DNAseq.= $dna;}
	}
	my $DNAsequence = uc($DNAseq); ## Changing to upper case to fit with the translation hash
	my $contig_length = length($DNAsequence); ## Calculating the contig size
	
	### Create hashes to adjust for codon_start features 
	my %coord = ();
	my @feat = ();
	my $feature;
	%codon_start = ();

	my $regex = 'CDS|rRNA|tRNA|tmRNA|ncRNA|snRNA|snoRNA';

	while (my $line = <IN>){
		chomp $line;
		if ($line =~ /^FT   gene/){ $feature = 'gene'; }
		elsif (($line =~ /FT                   \/locus_tag="(.*)\"/) && ($feature eq 'gene')){
			$locus_tag = $1;
			push(@feat, $locus_tag);
		}
		elsif ($line =~ /FT   ($regex)/){
			$feature = '$1';
			$coord{$locus_tag} = $line;
		}
		elsif ($line =~ /FT\s+\/codon_start=(\d+)/){
			$codon_start{$locus_tag} = $1;
		}
	}
	
	while ($locus = shift@feat){

		chomp $locus;

		if (exists $coord{$locus}){

			my $line = $coord{$locus}; ## Fetch corresponding coordinates line
			$exon_counter = 0;
			$intron_counter = 0;

			### Forward, single exon
			if ($line =~ /FT   ($regex)\s+(\d+)..(\d+)/){

				my $type = $1;
				my $start = $2; $start--;
				my $stop = $3; $stop--;
				my $length = ($stop - $start + 1);

				### Creating mRNA sequence
				$mRNA = substr($DNAsequence, $start, $length);
				sequence($mRNA, \*RNA);

				### Creating protein sequence
				if ($type eq 'CDS'){
					translate($mRNA);
					sequence($protein, \*PROT);
					$protein = undef;
				}
			}

			### Forward, multiple exons
			elsif ($line =~ /FT   ($regex)\s+join\((.*)\)/){

				my $type = $1;
				my @array = split(',', $2);
				my @coordinates;

				while (my $segment = shift@array){
					chomp $segment;
					if ($segment =~ /(\d+)..(\d+)/){
						my $start = $1;
						my $stop = $2;
						push (@coordinates, --$start);
						push (@coordinates, --$stop);
					}
				}

				### Creating exons if requested
				if ($exon_feature){
					for (my $i = 0; $i < $#coordinates; $i += 2){
						$exon_counter++;
						my $start = $coordinates[$i];
						my $end = $coordinates[$i + 1];
						my $exon = substr($DNAsequence, $start, $end-$start + 1 );
						sequence($exon, \*EXON);
					}
				}

				### Creating introns if requested
				if ($intron_feature){
					for (my $i = 1; $i < $#coordinates; $i += 2){
						$intron_counter++;
						my $start = $coordinates[$i] + 1;
						my $end = $coordinates[$i + 1] - 1;
						my $intron = substr($DNAsequence, $start, $end-$start + 1 );
						sequence($intron, \*INTRON);
					}
				}

				### Creating mRNA sequence
				$mRNA = undef;
				for (my $i = 0; $i < $#coordinates; $i += 2){
					my $start = $coordinates[$i];
					my $end = $coordinates[$i + 1];
					$mRNA .= substr($DNAsequence, $start, $end-$start + 1 );
				}
				sequence($mRNA, \*RNA);

				### Creating protein sequence
				if ($type eq 'CDS'){
					translate($mRNA);
					sequence($protein, \*PROT);
					$protein = undef;
				}
			}

			### Reverse, single exon
			elsif ($line =~ /FT   ($regex)\s+complement\((\d+)..(\d+)\)/){

				my $type = $1;
				my $start = $3; $start--;
				my $stop = $2; $stop--;

				### Creating mRNA sequence
				my $mRNA =  substr($DNAsequence, ($stop), ($start-$stop+1));
				reverse_complement($mRNA);
				sequence($mRNA, \*RNA);

				### Creating protein sequence
				if ($type eq 'CDS'){
					translate($mRNA);
					sequence($protein, \*PROT);
					$protein = undef;
				}
			}

			### Reverse, mutiple exons
			elsif ($line =~ /FT   ($regex)\s+complement\(join\((.*)/){

				my $type = $1;
				my @array = split(',',$2);
				my @coordinates;

				while (my $segment = shift@array){
					chomp $segment;
					if ($segment =~ /(\d+)..(\d+)/){
						my $strt = $1;
						my $stp = $2;
						push (@coordinates, --$strt);
						push (@coordinates, --$stp);
					}
				}

				### Creating exons if requested
				if ($exon_feature){
					$exon_counter = scalar (@coordinates / 2);
					for (my $i = 0; $i < $#coordinates; $i += 2){
						my $start = $coordinates[$i];
						my $end = $coordinates[$i + 1];
						my $exon = substr($DNAsequence, $start, $end - $start + 1 );

						reverse_complement($exon);
						sequence($exon, \*EXON);
						$exon_counter--;
					}
				}

				### Creating introns if requested
				if ($intron_feature){
					$intron_counter = ((scalar(@coordinates))/2) - 1;
					for (my $i = 1; $i < $#coordinates; $i += 2){
						my $start = $coordinates[$i] + 1;
						my $end = $coordinates[$i + 1] - 1;
						my $intron = substr($DNAsequence, $start, $end-$start + 1 );

						reverse_complement($intron);
						sequence($intron, \*INTRON);
						$intron_counter--;
					}
				}

				### Creating mRNA sequence
				$mRNA = undef;

				for (my $i = 0; $i < $#coordinates; $i += 2){
					my $start = $coordinates[$i];
					my $end = $coordinates[$i + 1];
					$mRNA .= substr($DNAsequence, $start, ($end - $start + 1));
				}

				reverse_complement($mRNA);
				sequence($mRNA, \*RNA);

				### Creating protein sequence
				if ($type eq 'CDS'){
					translate($mRNA);
					sequence($protein, \*PROT);
					$protein = undef;
				}
			}
		}
	}
}
close IN;
close PROT;

## subroutines
sub reverse_complement{
	$_[0] = reverse($_[0]);
	$_[0] =~ tr/ATGCRYSWKMBDHVatgcryswkmbdhv/TACGYRWSMKVHDBtacgyrwsmkvhdb/;
}

sub translate{
	my $seq = $_[0];

	### Adjusting codon_start feature, if required
	my $cdn_start = undef;
	
	if (exists $codon_start{$locus}){
		$cdn_start = $codon_start{$locus};
	}
	else{
		$cdn_start = 1;
	}
	$cdn_start--;
	### End of adjustment

	my $phase = 0 + $cdn_start;

	for(my $i = $phase; $i < (length($seq) - 5); $i += 3){
		my $codon = substr($seq, $i, 3);
		if (exists  $gcodes{$gc}{$codon}){ $protein .=  $gcodes{$gc}{$codon}; }
		else { $protein .= 'X'; }
	}
}

sub sequence{
	my ($sequence, $fh) = @_;
	my $length = length $sequence;
	my $header = $locus;
	
	## Changing units if proteins
	my $unit;
	if ($fh eq \*PROT) { $unit = 'aa'; }
	else { $unit = 'bp' ;}

	## Changing header if exon or intron
	if ($fh eq \*EXON) {
		$header = "$locus".'_exon_#'."$exon_counter";
	}
	if ($fh eq \*INTRON) {
		$header = "$locus".'_intron_#'."$intron_counter";
	}

	print $fh ">$header \[length=$length $unit\]\n";
	
	my @SEQUENCE = unpack ("(A60)*", $sequence);
	while (my $seq = shift@SEQUENCE){
		print $fh "$seq\n";
	}
}

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
