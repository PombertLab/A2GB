#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'EMBLtoTBL.pl';
my $version = '1.7';
my $updated = '2024-01-15';

use strict;
use warnings;
use Bio::SeqIO;
use File::Basename;
use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Converts EMBL files to NCBI TBL format for TBL2ASN
REQUIREMENTS	BioPerl's Bio::SeqIO module
NOTE		The EMBL (*.embl) and FASTA (*.fsa) files must be in the same folder.
		Requires locus_tags to be defined in the EMBL files.
		
USAGE		${name} \\
		  -id IITBIO \\
		  -p product_list.txt \\
		  -embl *.embl \\
		  -c 1

OPTIONS:
-id			Desired institute ID [default: IITBIO]
-p (--prod)		Tab-delimited list of locus_tags and their products
-e (--embl)		EMBL files to convert
-c (--gcode)		NCBI genetic code [Default: 1]
			1  - The Standard Code
			2  - The Vertebrate Mitochondrial Code
			3  - The Yeast Mitochondrial Code
			4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
			11 - The Bacterial, Archaeal and Plant Plastid Code
			# For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
-o (--organelle)	Organelle mode; turns off protein_id / transcript_id
-t (--partial)		Add partials flags (< and >) to all gene/mRNA features [Default: off]
OPTIONS
die "\n$usage\n" unless @ARGV;

my $instID = 'IITBIO'; ## 
my $products; ## protein_list.txt
my @embl;
my $gcode = 1;
my $organelle;
my $partial;
GetOptions(
	'id=s' => \$instID,
	'p|prod=s' => \$products,
	'e|embl=s@{1,}' => \@embl,
	'c|gcode=i' => \$gcode,
	'o|organelle' => \$organelle,
	't|partial' => \$partial
);

## Checking for input files
unless (@embl){
	print "\nPlease specify at least one EMBL file with the -e option...\n\n";
	exit;
}

### Filling the products database
my %hash = ();
open HASH, "<", $products or die "Can't open $products: $!\n";
while (my $dbkey = <HASH>){
	chomp $dbkey;
	if ($dbkey =~ /^(\S+)\t(.*)$/){
		my $prot = $1;
		my $prod = $2;
		$hash{$prot} = $prod;
	}
}

### Working on EMBL files
my $locus_tag;
while (my $file = shift@embl){

	open IN, "<", $file or die "Can't open EMBL file file: $file\n";
	$file =~ s/.embl$//;
	my ($head, $dir) = fileparse($file);

	open DNA, "<", "$file.fsa" or die "Can't open FASTA file file: $file.fsa\n";
	open TBL, ">", "$file.tbl" or die "Can't create TBL output file: $file.tbl\n";
	print TBL ">Feature $head\n"; ## Generate TBL header
	
	### Creating a single DNA string for codon verification
	my $DNAseq = undef;
	while (my $dna = <DNA>){
		chomp $dna;
		if ($dna =~ /^>/){ next; }
		else { $DNAseq.= $dna; }
	}
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
		if ($line =~ /^FT\s+\/locus_tag="(\S+)"/){ $locus_tag = $1; }
		## If a pseudogene is present, it is assumed that the start and stop, as well as
		## pseudogene note, is provided in the "" of the feature. The start and stop are
		## seperated by ".." and the range is separated from the note by ";". If the
		## pseudogene is present on the complementary strand, reverse the start and stop.
		## e.g.: /pseudogene="21011..19202; similar to ECU_06_0080"
		elsif($line =~ /^FT\s+\/pseudogene/){
			my $start;
			my $stop;
			my $note="";
			if($line =~ /"complement(\d+)\.\.(\d+);(.*?)"/){
				$start = $2;
				$stop = $1;
				$note = $3;
			}
			elsif($line =~ /"(\d+)\.\.(\d+);(.*?)"/){
				$start = $1;
				$stop = $2;
				$note = $3;
			}
			print TBL $start."\t".$stop."\tgene\n";
			print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
			print TBL "\t\t\tpseudo\n";
			print TBL "\t\t\tnote\t$note\n";
		}
		### Working on tRNAs/rRNAs
		elsif ($line =~ /^FT\s+(tRNA|rRNA)\s+(\d+)..(\d+)/){ ## Forward, single exon
			my $type = $1;
			my $start = $2;
			my $stop = $3;
			print TBL $start."\t".$stop."\tgene\n";
			print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
			print TBL $start."\t".$stop."\t".$type."\n";
			RNA();
		}
		elsif ($line =~ /^FT\s+(tRNA|rRNA)\s+join\((.*)\)/){ ## Forward, multiple exons
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
			print TBL $start[0]."\t".$stop[$num]."\tgene\n";
			print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
			
			### Printing tRNA/rRNA info
			print TBL $start[0]."\t".$stop[0]."\t".$feat."\n"; ## printing the 1st exon
			if ($asize == 2){
				print TBL $start[$num]."\t".$stop[$num]."\n";
			}
			elsif ($asize >= 3){
				foreach my $subs (1..$dum){
					print TBL $start[$subs]."\t".$stop[$subs]."\n";
				}
				print TBL $start[$num]."\t".$stop[$num]."\n"; ## printing the last exon	
			}
			RNA();
		}
		elsif ($line =~ /^FT\s+(tRNA|rRNA)\s+complement\((\d+)..(\d+)\)/){ ## Reverse, single exon
			my $type = $1;
			my $start = $3;
			my $stop = $2;
			print TBL $start."\t".$stop."\tgene\n";
			print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
			print TBL $start."\t".$stop."\t".$type."\n";
			RNA();
		}
		elsif ($line =~ /^FT\s+(tRNA|rRNA)\s+complement\(join\((.*)/){ ## Reverse, mutiple exons
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
			print TBL $start[0]."\t".$stop[$num]."\tgene\n";
			print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
			
			### Printing tRNA/rRNA info
			print TBL $start[0]."\t".$stop[0]."\t".$feat."\n"; ## printing the 1st exon
			if ($asize == 2){
				print TBL $start[$num]."\t".$stop[$num]."\n";
			}
			elsif ($asize >= 3){
				foreach my $subs (1..$dum){
					print TBL $start[$subs]."\t".$stop[$subs]."\n";
				}
				print TBL $start[$num]."\t".$stop[$num]."\n"; ## printing the last exon
			}
			RNA();
		}
		
		### Working on CDS
		elsif ($line =~ /^FT\s+\/codon_start=(\d)/){ ## Looking for phased codons
			print TBL "\t\t\tcodon_start\t$1\n";
		}
		elsif ($line =~ /^FT\s+CDS\s+(\d+)..(\d+)/){ ## Forward, single exon

			my $start = $1;
			my $stop = $2;
			my $startcodon = substr($DNAsequence, $start-1, 3);
			my $stopcodon = substr($DNAsequence, $stop-3, 3);

			my $left_b = '';
			my $right_b = '';
			
			if (($startcodon ne 'atg') || ($start == 1)){
				$left_b = '<'
			}

			unless (($stopcodon eq 'taa') || ($stopcodon eq 'tag') || ($stopcodon eq 'tga')){
				$right_b = '>';
			}
			else {
				if ($stop == $contig_length){
					$right_b = '>';
				}
			}

			## Gene + mRNA tag
			if ($partial){
				print TBL '<'.$start."\t".'>'.$stop."\tgene\n";
				print TBL "\t\t\tlocus_tag\t$locus_tag\n";
				print TBL '<'.$start."\t".'>'.$stop."\tmRNA\n";
			}
			else{
				print TBL $left_b.$start."\t".$right_b.$stop."\tgene\n";
				print TBL "\t\t\tlocus_tag\t$locus_tag\n";
				print TBL $left_b.$start."\t".$right_b.$stop."\tmRNA\n";
			}
			mRNA();

			## CDS tag
			print TBL $left_b.$start."\t".$right_b.$stop."\tCDS\n";
			CDS();
		}

		## Reverse, single exon
		elsif ($line =~ /^FT\s+CDS\s+complement\((\d+)..(\d+)\)/){

			my $start = $2;
			my $stop = $1;
			my $startcodon = substr($DNAsequence, $start-3, 3);
			my $stopcodon = substr($DNAsequence, $stop-1, 3);

			my $left_b = '';
			my $right_b = '';
			
			if (($start == $contig_length) or ($startcodon ne 'cat')){
				$left_b = '<';
			}
			if ($stop == 1){
				$right_b = '>';
			}
			unless (($stopcodon eq 'tta') || ($stopcodon eq 'cta') || ($stopcodon eq 'tca')){
				$right_b = '>';
			}

			## Gene + mRNA tag
			if ($partial){
				print TBL '<'.$start."\t".'>'.$stop."\tgene\n";
				print TBL "\t\t\tlocus_tag\t$locus_tag\n";
				print TBL '<'.$start."\t".'>'.$stop."\tmRNA\n";
			}
			else{
				print TBL $left_b.$start."\t".$right_b.$stop."\tgene\n";
				print TBL "\t\t\tlocus_tag\t$locus_tag\n";
				print TBL $left_b.$start."\t".$right_b.$stop."\tmRNA\n";
			}
			mRNA();

			## CDS tag
			print TBL $left_b.$start."\t".$right_b.$stop."\tCDS\n";
			CDS();
		}

		## Forward, multiple exons
		elsif ($line =~ /^FT\s+CDS\s+join\((.*)\)/){
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

			my $startcodon = substr($DNAsequence, $start[0]-1, 3);
			my $stopcodon = substr($DNAsequence, $stop[$num]-3, 3);

			my $left_b = '';
			my $right_b = '';

			if (($startcodon ne 'atg') || ($start[0] == 1)){
				$left_b = '<';
			}
			if ($stop[$num] == $contig_length){
				$right_b = '>';
			}
			unless (($stopcodon eq 'taa') || ($stopcodon eq 'tag') || ($stopcodon eq 'tga')){
				$right_b = '>';
			}

			## Gene + mRNA info
			my $lb = '';
			my $rb = '';

			if ($partial){
				$lb = '<';
				$rb = '>';
			}
			else{
				$lb = $left_b;
				$rb = $right_b;
			}

			### Printing gene info
			print TBL $lb.$start[0]."\t".$rb.$stop[$num]."\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			
			### Printing mRNA info
			if ($asize == 1){ 
				print TBL $lb.$start[0]."\t".$rb.$stop[0]."\tmRNA\n";
			}
			elsif ($asize == 2){ 
				print TBL $lb.$start[0]."\t".$stop[0]."\tmRNA\n";
				print TBL $start[1]."\t".$rb.$stop[1]."\n";
			}
			elsif ($asize >= 3){
				print TBL $lb.$start[0]."\t".$stop[0]."\tmRNA\n";
				foreach my $subs (1..$dum){
					print TBL $start[$subs]."\t".$stop[$subs]."\n";
				}
				print TBL $start[$num]."\t".$rb.$stop[$num]."\n";
			}
			mRNA();

			### Printing CDS info
			print TBL $left_b.$start[0]."\t".$stop[0]."\tCDS\n";

			if ($asize == 2){
				print TBL $start[$num]."\t".$right_b.$stop[$num]."\n";
			}

			elsif ($asize >= 3){
				foreach my $subs (1..$dum){
					print TBL $start[$subs]."\t".$stop[$subs]."\n";
				}
				print TBL $start[$num]."\t".$right_b.$stop[$num]."\n";
			}
			CDS();
		}

		## Reverse, mutiple exons
		elsif ($line =~ /^FT\s+CDS\s+complement\(join\((.*)/){ 
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

			my $startcodon = substr($DNAsequence, $start[0]-3, 3);
			my $stopcodon = substr($DNAsequence, $stop[$num]-1, 3);

			my $left_b = '';
			my $right_b = '';

			if (($startcodon ne 'cat') || ($start[0] == $contig_length)){
				$left_b = '<';
			}
			if ($stop[$num] == 1){
				$right_b = '>';
			}
			unless (($stopcodon eq 'tta') || ($stopcodon eq 'cta') || ($stopcodon eq 'tca')){
				$right_b = '>';
			}

			## Gene + mRNA info
			my $lb = '';
			my $rb = '';

			if ($partial){
				$lb = '<';
				$rb = '>';
			}
			else{
				$lb = $left_b;
				$rb = $right_b;
			}

			### Printing gene info
			print TBL $lb.$start[0]."\t".$rb.$stop[$num]."\tgene\n";
			print TBL "\t\t\tlocus_tag\t$locus_tag\n";
			
			### Printing mRNA info
			if ($asize == 1){
				print TBL $lb.$start[0]."\t".$rb.$stop[0]."\tmRNA\n";
			}
			elsif ($asize == 2){
				print TBL $lb.$start[0]."\t".$stop[0]."\tmRNA\n";
				print TBL $start[1]."\t".$rb.$stop[1]."\n";
			}
			elsif ($asize >= 3){
				print TBL $lb.$start[0]."\t".$stop[0]."\tmRNA\n";
				foreach my $subs (1..$dum){
					print TBL $start[$subs]."\t".$stop[$subs]."\n";
				}
				print TBL $start[$num]."\t".$rb.$stop[$num]."\n";
			}
			mRNA();

			### Printing CDS info
			print TBL $left_b.$start[0]."\t".$stop[0]."\tCDS\n";
			if ($asize == 2){
				print TBL $start[$num]."\t".$right_b.$stop[$num]."\n";
			}
			elsif ($asize >= 3){
				foreach my $subs (1..$dum){
					print TBL $start[$subs]."\t".$stop[$subs]."\n";
				}
				print TBL $start[$num]."\t".$right_b.$stop[$num]."\n";
			}

			CDS();

		}
	}
}
close IN;
close TBL;

### Subroutines
sub numSort {
	if ($a < $b){ return -1; }
	elsif ($a == $b){ return 0; }
	elsif ($a > $b) { return 1; }
}

sub mRNA {
	print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
	if (exists $hash{$locus_tag}){
		print TBL "\t\t\tproduct\t".$hash{$locus_tag}."\n";
	}
	else{
		print TBL "\t\t\tproduct\thypothetical protein\n";
	}
	unless ($organelle){
		print TBL "\t\t\tprotein_id\tgnl|$instID|$locus_tag\n";
		print TBL "\t\t\ttranscript_id\tgnl|$instID|$locus_tag"."_mRNA\n";
	}
}

sub CDS {
	print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
	if (exists $hash{$locus_tag}){
		print TBL "\t\t\tproduct\t".$hash{$locus_tag}."\n";
	}
	else{
		print STDERR "Cannot find database entry for locus_tag: $locus_tag\n";
		print TBL "\t\t\tproduct\thypothetical protein\n";
	}
	unless ($gcode == 1){
		print TBL "\t\t\ttransl_table\t$gcode\n";
	}
	unless ($organelle){
		print TBL "\t\t\tprotein_id\tgnl|$instID|$locus_tag\n";
		print TBL "\t\t\ttranscript_id\tgnl|$instID|$locus_tag"."_mRNA\n";
	}
}

sub RNA {
	print TBL "\t\t\tlocus_tag\t".$locus_tag."\n";
	if (exists $hash{$locus_tag}){
		print TBL "\t\t\tproduct\t".$hash{$locus_tag}."\n";
	}
	else{
		print STDERR "Cannot find database entry for locus_tag: $locus_tag\n";
		print TBL "\t\t\tproduct\thypothetical RNA\n";
	}
}
