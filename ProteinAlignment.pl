#!/usr/bin/perl

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $name = 'ProteinAlignmentV2.pl';
my $version = '0.2a';
my $updated = '03/27/2021';
my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	The purpose of this script is to complete a blastp search on a protein of interest, and return a complete
			visual alignment of the query against the corresponding blastp results.

COMMAND		${name} \\
			-q protein.prot \\
			-i 6035 \\
			-t 16 \\
			-c 3 \\
			-v

OPTIONS

-q | --query			Input faa file
-i | --taxids			Taxonomic IDs to refine search (comma-separated list or new-line-separated file)
-t | --threads			Number of threads to run on [default = 8]
-c | --culling_limit	Number of results to return [default = 5]
-e | --evalue			Cutoff evalue for results [default = 1e-10]
-o | --output			Name of output file [default = aligned.log]
-q | --quiet			Prevent output to STDOUT as well as to output file [default = on]
EXIT
die("\n\n$usage\n\n") unless(@ARGV);

my @query;
my $taxid;
my $threads = 8;
my $cul_lim = 5;
my $eval = '1e-10';
my $out = 'aligned.log';
my $db = 'nr';
GetOptions(
	'q|query=s@{1,}' => \@query,
	'i|taxids=s' => \$taxid,
	't|threads=s' => \$threads,
	'c|culling_limit=s' => \$cul_lim,
	'e|evalue=s' => \$eval,
	'o|output=s' => \$out
);
die("[E] Query required to run blast") unless(@query);

## Create hashes of conserved amino acid substitutions
my %VHCS;
my %HCS;
my %CS;
Substitutions();

## Check to see if taxids are passed through a file or as comma seperated list
my $ids;
if(-f $taxid){
	open ID,"<","$taxid";
	while(my $line = <ID>){
		chomp($line);
		if($line =~ /,/){
			$ids .= $line;
		}
		else{
			$ids .= "$line,";
		}
	}
}
else{
	$ids = $taxid;
}

## Create output directory
unless(-d "align_files"){
	mkdir("align_files",0755);
}

## Iterate through all protein files
foreach my $file (@query){

	## blastp using taxids if provided
	if($taxid){
		system("blastp \\
				-task blastp-fast \\
				-query $file \\
				-db $db \\
				-taxids $ids \\
				-num_threads $threads \\
				-culling_limit $cul_lim \\
				-evalue $eval \\
				-outfmt 0 \\
				-out align_files/$file.blastp.0"
		);
	}
	## blastp without taxids
	else{
		system("blastp \\
				-task blastp-fast \\
				-query $file \\
				-db $db \\
				-num_threads $threads \\
				-culling_limit $cul_lim \\
				-evalue $eval \\
				-outfmt 0 \\
				-out align_files/$file.blastp.0"
		);
	}

	## Get the protein amino acid sequence 
	open QRY,"<","$file";
	my $qry_seq;
	while (my $line = <QRY>){
		chomp($line);
		unless($line =~ /^>/){
			$qry_seq .= $line;
		}
	}
	close QRY;

	## Parse blast file for the following details:
	## Match Data
	## [0] Accession number
	## [1] Protein annotation
	## [2] Organsim name
	## [3] Query start position
	## [4] Query stop position
	## [5] Reference start position
	## [6] Reference stop position
	## [7] Reference sequence
	open BLAST,"<","align_files/$file.blastp.0" or die"$!\n";
	my $accession;
	my $count = 0;
	my %match_data;
	my $prev_accession;
	my $grab_rest_of_line;
	while(my $line = <BLAST>){
		chomp($line);
		## Get the protein accession number, the protein annotation, and the organism name
		my $ref_seq;
		## If the organism name moves onto a second line, grab the rest of the organism line
		if($grab_rest_of_line){
			if($line =~ /(.+)/){
				@{$match_data{$accession}}[3] .= $1;
			}
			$grab_rest_of_line = undef;
		}
		if($line =~ /^>(\S+)\s+(.+)\s(\[.+)/){
			## $1 Accession number
			## $2 Protein Annotation
			## $3 Organism Name
			$accession = $1;
			if($prev_accession){
				if($prev_accession ne $accession){
					$count = 0;
				}
			}
			## Get the amino acid sequence for the blastp result
			system("blastdbcmd \\
				-entry $accession \\
				-db $db \\
				-out align_files/$accession.prot");
			## Populate reference protein database with amino acid sequence
			open P,"<","align_files/$accession.prot";
			$ref_seq = "";
			while(my $line = <P>){
				chomp($line);
				if($line =~ /^>(\S+)/){
					next;
				}
				$ref_seq .= $line;
			}
			close P;
			## [0] Reference Sequence
			push(@{$match_data{$accession}},$ref_seq);
			## [1] Accession number
			push(@{$match_data{$accession}},$1);
			## [2] Protein annotation
			push(@{$match_data{$accession}},$2);
			## [3] Organism
			push(@{$match_data{$accession}},$3);

			unless($3 =~ /\[.+?\]/){
				$grab_rest_of_line = 1;
			}
			else{
				$grab_rest_of_line = undef;
			}

		}
		if($line =~ /Score = (\d+)/){
			## [4] Bitscore
			push(@{$match_data{$accession}},$1);
		}
		## If we are on a query line, but this is the first time reading it, store the start point
		if($line =~ /^Query\s+(\d+)\s+(\S+)\s+(\d+)/){
			if($count < 2){
				## [5] Query start
				push(@{$match_data{$accession}},$1);
				## [6] Query end
				push(@{$match_data{$accession}},$3);
				## [7] Query data
				push(@{$match_data{$accession}},$2);
			}
			else{
				## [7] Query data
				@{$match_data{$accession}}[7] .= $2;
				## [6] Query end
				@{$match_data{$accession}}[6] = $3;
			}
			$count++;
		}
		if($line =~ /^Sbjct\s+(\d+)\s+(\S+)\s+(\d+)/){
			if($count < 2){
				## [8] Ref start
				push(@{$match_data{$accession}},$1);
				## [9] Ref stop
				push(@{$match_data{$accession}},$3);
				## [10] Reference data
				push(@{$match_data{$accession}},$2);
			}
			else{
				## [10] Ref data
				@{$match_data{$accession}}[10] .= $2;
				## [9] Ref stop
				@{$match_data{$accession}}[9] = $3;
			}
			$count++;
		}
		$prev_accession = $accession;
	}
	close BLAST;

	## Iterate through all the blastp results, sorted from best score to worst score
	foreach $accession (sort{$match_data{$b}[4] <=> $match_data{$a}[4]} keys %match_data){
		my @data = @{$match_data{$accession}};
		print("\n\n");
		Align(\@data,$qry_seq,$file);
		print("\n\n");
	}

}

sub Align{

	## Define input variables
	
	my @data = @{$_[0]};
	## Query sequence
	my $q_seq = $_[1];
	## Reference sequence
	my $r_seq = $data[0];
	## Accession number
	my $access = $data[1];
	## Protein annotation
	my $p_a = $data[2];
	## Organism name
	my $o_n = $data[3];
	## Query start
	my $q_s = $data[5];
	## Query end
	my $q_e = $data[6];
	## Query blast sequence
	my $q_b_seq = $data[7];
	## Reference start
	my $r_s = $data[8];
	## Reference end
	my $r_e = $data[9];
	## Reference blast sequence
	my $r_b_seq = $data[10];

	## Define variables needed for printing

	## Edited query sequence
	my $e_q_seq;
	## Edited reference sequence
	my $e_r_seq;
	## Edited alignment sequence
	my $e_a_seq;

	## Difference between reference start and query start
	my $q_dif = $r_s - $q_s;
	## Difference between query start and reference start
	my $r_dif = $q_s - $r_s;

	## If the query sequence starts before the reference, start printing at the beginning of the query
	if($q_dif > 0){
		$q_dif = 0;
	}
	## If the reference sequence starts before the query, start printing at the beginning of the reference
	if($r_dif > 0){
		$r_dif = 0;
	}

	## Index of the overall sequence
	my $i = 0;
	## Index of the blast alignmnet sequence
	my $b_i = 0;

	## The buffer offset for the query
	my $q_b = 0;
	## The buffer offset for the reference
	my $r_b = 0;

	## Run until the entirity of both sequences have been printed
	WHILE: while(0 == 0){
	
		## Stores whether looking into blast sequence or pure sequence
		my $blast;

		my $q_v = " ";
		my $r_v = " ";

		## If within the blast alignment, use the blast alignment values for both query and reference
		if($q_dif + $i >= $q_s - 1 && $q_dif + $i <= $q_e - 1){
			$e_q_seq .= substr($q_b_seq,$b_i,1);
			$q_v = substr($q_b_seq,$b_i,1);
			if(substr($q_b_seq,$b_i,1) eq "-"){
				$q_b++;
			}
			$e_r_seq .= substr($r_b_seq,$b_i,1);
			$r_v = substr($r_b_seq,$b_i,1);
			if(substr($r_b_seq,$b_i,1) eq "-"){
				$r_b++;
			}
			$blast = 1;
			$b_i++;
		}
		## If not in alignment, but in pure sequence, use the pure sequence value for the query
		elsif($q_dif + $i - $q_b > -1 && $q_dif + $i - $q_b < length($q_seq) - 1){
			$e_q_seq .= substr($q_seq,$q_dif+$i-$q_b,1);
			$q_v = substr($q_seq,$q_dif+$i-$q_b,1);
		}
		## If not in alignment, or pure sequence, use a space as the string value
		else{
			$e_q_seq .= " ";
		}

		unless($blast){
			## If not in alignment, but in pure sequence, use the pure sequence value for the reference
			if($r_dif + $i - $r_b > -1 && $r_dif + $i - $r_b < length($r_seq)){
				$e_r_seq .= substr($r_seq,$r_dif+$i-$r_b,1);
				$r_v = substr($r_seq,$r_dif+$i-$r_b,1);
			}
			## If not in alignment, or pure sequence, use a space as the string value
			else{
				$e_r_seq .= " ";
			}
		}

		## Create a string of possible conserved substitutions
		my $subs = $VHCS{$q_v}."/".$HCS{$q_v}."/".$CS{$q_v};
	
		## Determine the alignment pattern between the two sequences; can either be exact match, conserved match, or no
		## match
		if($q_v eq $r_v && $q_v ne " "){
			$e_a_seq .= "|";
		}
		elsif($subs =~ /$r_v/){
			$e_a_seq .= "+";
		}
		else{
			$e_a_seq .= " ";
		}

		## If the end of both sequences has been reaced, exit the loop
		if($q_v eq $r_v && $q_v eq " "){
			last WHILE;
		}

		## Add a space buffer between every ten amino acids
		if(($i+1)%10==0){
			$e_q_seq .= " ";
			$e_r_seq .= " ";
			$e_a_seq .= " ";
		}

		## Move to the next index point
		$i++;

	}

	## Split the three sequences into lines of 100 amino acids
	my @q_seq = unpack("(A110)*",$e_q_seq);
	my @r_seq = unpack("(A110)*",$e_r_seq);
	my @align = unpack("(A110)*",$e_a_seq);
	## Determine the padding size based on the longest sequence
	my $padding = length("$i");
	## Create a line buffer for a prettier print
	my $buffer = "-" x 120;
	## Print match data
	print("$access\t$p_a\t$o_n\n");
	for(my $i = 0;$i < scalar(@q_seq); $i++){
		my $locale = ($i*100)+1;
		my $location = sprintf("%0${padding}d",$locale);
		my $q_seq_line = $q_seq[$i];
		my $r_seq_line = $r_seq[$i];
		my $a_seq_line = $align[$i];
		print("$buffer\n");
		print("$location\n");
		print("$buffer\n");
		print("qry\t$q_seq_line\n");
		print("   \t$a_seq_line\n");
		print("ref\t$r_seq_line\n");
	}
	print("$buffer\n");

}

## Conserved substitution creation
sub Substitutions{
	%VHCS = ("A" => "S", "R" => "K", "N" => "Q/H", "D" => "E", "C" => "S", "Q" => "N", "E" => "D", "G" => "P",
			"H" => "N/Q", "I" => "L/V", "L" => "I/V", "K" => "R/Q/E", "M" => "L/I", "F" => "M/L/Y", "S" => "T",
			"T" => "S","W" => "Y", "Y" => "W/F", "V" => "I/L", "P" => "1", " " => "-1", "-" => "1"
	);
	%HCS = ("A" => "G/S/T", "R" => "Q/H/K", "N" => "D/Q/H/K/S/T", "D" => "N/Q", "C" => "", "Q" => "R/N/E/H/K/M", 
			"E" => "D/Q/K", "G" => "A",	"H" => "R/N/Q/Y", "I" => "L/M/V", "L" => "I,M,F,V", "K" => "R/D/Q/E", 
			"M" => "Q/I/L/V", "F" => "L/W/Y", "S" => "A/N/T", "T" => "A/N/S", "W" => "F/Y", "Y" => "H/F/W", 
			"V" => "I/L/M", "P" => "-1", " " => "1", "-" => "1"
	);
	%CS = ("A" => "C/G/S/T/V", "R" => "N/Q/E/H/K", "N" => "R/D/Q/E/H/K/S/T", "D" => "N/Q/E/S", "C" => "A", 
			"Q" => "R/N/D/E/H/K/M/S", "E" => "R/N/D/Q/H/K/S", "G" => "A/S", "H" => "R/N/Q/E/Y", "I" => "K/M/F/V", 
			"L" => "I/M/F/V", "K" => "R/N/Q/E/S", "M" => "Q/I/L/F/V", "F" => "I/L/M/W/Y", "S" => "A/N/D/Q/E/G/L/T",
			"T" => "A/N/S/V", "W" => "F/Y", "Y" => "H/F/W", "V" => "A/I/L/M/T", "P" => "-1", " " => "1", "-" => "1"
	);
}