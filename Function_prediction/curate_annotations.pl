#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'curate_annotations.pl';
my $version = '1.6a';
my $updated = '2021-04-09';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Displays lists of functions predicted per proteins. User can select or enter desired annotation.
		Creates a tab-delimited .curated list of annotations.

COMMAND		${name} \\
		  -i proteins.annotations \\
		  -r

OPTIONS
-i (--input)	Sequence homology based annotations (generated from parse_annotators.pl)
-r (--resume)		Resume annotation from last curated locus_tag
-c (--check)		Check loci marked with '?'
-d (--3D_annot)		3D structural homology based annotations (Generated with descriptive_GESAMT_matches.pl)
EXIT
die "\n\n$usage\n\n" unless @ARGV;

my $input;
my $input_3D;
my $continue;
my $review;

GetOptions(
	"i|input=s" => \$input,
	"r|resume" => \$continue,
	"c|check" => \$review,
	"d|3D_annot=s" => \$input_3D
);

## Make a temporary directory to store our temporary files.
unless (-d 'temp_files'){
	mkdir ('temp_files',0755) or die "Can't make directory temp_files: $!\n";
}

my ($filename,$dir) = fileparse($input);
open OUT, ">", "temp_files/$filename.temp" or die "Can't create temp_files/$filename.temp: $!\n";

## If reviewing annotations, load existing annotations into RAM
my $last_locus = undef;
my @to_review;
if ($review){
	open IN, "<", "$filename.curated" or die "Can't open $filename.curated: $!\n";
	while (my $line = <IN>){
		chomp $line;
		push (@to_review,$line);
	}
	close IN;
}
else{
	## Now, regardless of whether the "continue" flag was passed, check to see if annotations have taken place to not
	## inform the user that they are going to end up with a file with twice annotated locus tags.
	if (-f "$filename.curated"){
		unless ($continue){
			WHILE: while (0==0){
				print "\nWARNING: It appears that a curation file already exists for the provided annotation file.";
				print "If you proceed, you will overwrite annotated loci.\n\nWould you like to continue? \n\n[y]/[n]: ";
				chomp (my $proceed = lc(<STDIN>));
				if ( $proceed eq 'y' ) { last WHILE; }
				elsif ($proceed eq 'n') { die "\nScript terminating...\n\n"; }
				system "clear";
				system "mv $filename.curated $filename.backup.curated";
				print "\n $proceed is an invalid operator.\n";
			}
		}
		else{
			open IN, "<", "$filename.curated" or die "Can't read $filename.curated: $!\n";
			while (my $line = <IN>){
				chomp $line;
				print OUT "$line\n";
				my @data = split("\t",$line);
				$last_locus = $data[0];
			}
			close IN;
		}
	}
}

## If we are given a 3D file, parse through it to get the match data
my %three_d; my $three_d_query;
if ($input_3D){
	open DD, "<", "$input_3D" or die "\n[E] Can't open 3D matches file: $input_3D\n\n";
	while (my $line = <DD>){
		chomp $line;
		if ($line =~ /^### .*?(\w+)\-\w+\-\w+$/){ $three_d_query = $1; }
		else {
			if ($line =~ /^\S+\s+\d+\s+(\w+)\s+\w\s+(\S+).*pdb\w+.ent.gz\s(.*)$/){
				my $hit = $1;
				my $qscore = $2;
				my $match = $3;
				push (@{$three_d{$three_d_query}}, "$qscore\t$hit - $match");
			}
		}
	}
}

## Get total protein count. Start at -1 because, you know, file header
open IN, "<", "$input" or die "Can't read $input: $!\n";
my $protein_total = -1;
while (my $line = <IN>){
	$protein_total++;
}
my $padding = length($protein_total);
close IN;

## Do we start curating at the beginning? undef == no, 1 == yes. If we have a $last_locus, we aren't going to start at
## the beginning, so start == undef
my $current_protein = sprintf("%0${padding}d",0);
my @header_info;
my $start = 1;
my $string_length;
my $tab;
open IN, "<", $input or die "Can't read $input: $!\n";
if ($last_locus) { $start = undef; }
while (my $line = <IN>){
	chomp $line;
	## If we are at the header of the file, get the header info and skip to the next line
	if ($line =~ /^#/){
		$line =~ s/#//;
		@header_info = split("\t",$line);
		next;
	}
	$current_protein++;
	
	## If we are reviewing, skip the check start process because we are not going to be curating the same way.
	unless ($review){
		## Check to see if we are the $last_locus, and if we are, start the anntotation process
		unless ($start){
			if ($line =~ /$last_locus/){
				$start = 1;
			}
			next;
		}
	}
	else {
		my $val = scalar(@to_review);
		print "Values left to review $val\n";
		## If there are still more loci to review, check for review necessity
		if (@to_review){
			my $review_line = shift(@to_review);
			## Check to see if we are at an annotation that needs review, indicated by a '?'
			## If we are reading a line that does not require reviewing, print the line out, and move to the next one.
			## If we are reading a line that does require reviewing, print the normal annotation curation printouts.
			unless ($review_line =~ /\?/){
				print OUT "$review_line\n";
				next;
			}
		}
		## If there are no more loci to review, let user know and exit script.
		else{
			print "\nAnnotation review completed.\n\n";
			exit;
		}
	}
	## Locus	E-Value_1	Annot_1	E-Value_2	Annot_2	etc...
	## Odds are evalues, evens are predictions, except for the first even, which is the locus tag
	my @data = split("\t",$line);
	my $locus = $data[0];
	my @sources;
	my @predictions;
	my @evalues;
	my $NA_count = 0;
	for (my $i = 1; $i < scalar(@data); $i++){
		if ($i%2 == 0){
			push(@predictions,$data[$i]);
			push(@sources,$header_info[$i]);
		}
		else {
			push(@evalues,$data[$i]);
			if ("$data[$i]" eq "NA"){
				$NA_count++;
			}
		}
	}
	## Add 3D information if provided
	my $three_d_sources = 0;
	if (defined $input_3D){
		if (defined $three_d{$locus}){
			foreach my $struct (@{$three_d{$locus}}){
				push(@sources,"GESAMT");
				my @data = split ("\t", $struct);
				push(@predictions,$data[1]);
				push(@evalues,$data[0]);
				# print "3D.\tGESAMT:\t\t"."$_"."\n";
				$three_d_sources++;
			}
		}
		else {
			push(@sources,"GESAMT");
			push(@predictions,"no match found");
			push(@evalues,"NA");
			$NA_count++;
			$three_d_sources++;
			# print "3D.\tGESAMT:\t\t"."NA\t\tno match found"."\n";
		}
	}

	## Count the number of NAs present. If the number of NAs == number of predictions, the protein is hypothetical
	if ($NA_count == scalar(@evalues)){
		print OUT "${locus}\thypothetical protein\n";
		next;
	}

	## Make a lovely progress bar so our users know just how much suffering they have left to go
	my $progress = "|" x (int(($current_protein/$protein_total)*100));
	my $remaining = "." x (100-int(($current_protein/$protein_total)*100));
	my $status = "[".$progress.$remaining."]";
	my $options = scalar(@sources);
	my $choices = $options - $three_d_sources;
	WHILE: while (0==0){
		my $status_3D = 0;
		print "\n$status\t$current_protein/$protein_total\n";
		print "\nPutative annotation(s) found for protein $data[0]:\n";
		## Loop through our information arrays so we don't have a million conditionals, and we are more robust this way
		for (my $i = 1; $i <= $options; $i++){
			my $source = $sources[$i-1];
			my $prediction = $predictions[$i-1];
			my $evalue = $evalues[$i-1];
			$string_length = length("$source"); tab();
			if ($i <= $choices){
				print "${i}.\t${source}${tab}";
				$string_length = length("$evalue"); tab();
				print "${evalue}${tab}${prediction}\n";
			}
			else {
				if ($status_3D == 0){
					print "\n".'## 3D structural homologs (if any)'."\n";
					$status_3D = 1;
				}
				print "3D.\t${source}${tab}";
				$string_length = length("$evalue"); tab();
				print "${evalue}${tab}${prediction}\n";
			}
		}
		print "\nPlease enter [1-$choices] to assign annotation, [0] to annotate the locus as a 'hypothetical protein', ";
		print "[m] to manually annotate the locus, [?] to mark this annotation for review, or [x] to exit.\n";
		print "Selection: ";
		chomp (my $select = <STDIN>);
		if ($select eq 'x'){
			if ($review){
				print OUT "$locus\t?\n";
			}
			cleanup();
			die "Exiting annotation curation..." ;
		}
		elsif ($select eq 'm'){
			print "Enter desired annotation: ";
			chomp (my $manual = <STDIN>);
			print OUT "$locus\t$manual\n";
			system "clear";
			last WHILE;
		}
		elsif ($select eq '?'){
			print OUT "$locus\t?\n";
			system "clear";
			last WHILE;
		}
		elsif ($select eq '0'){
			print OUT "$locus\thypothetical protein\n";
			system "clear";
			last WHILE;
		}
		elsif ($select =~ /^\d+$/){
			my $selection = int($select)-1;
			if ($selection < $choices){
				print OUT "$locus\t$predictions[$selection]\n";
				system "clear";
				last WHILE;
			}
			else {
				system "clear";
				print "\nERROR: Invalid input value '$select'.\n\n";
			}
		}
		else {
			system "clear";
			print "\nERROR: Invalid input value '$select'.\n\n";
		}
	}
}
cleanup();

##subroutines

sub tab{
	if ($string_length >= 8){ $tab = "\t"; }
	else { $tab = "\t\t"; }
}

sub cleanup{
	if ($review){
		while (my $line = shift(@to_review)){
			print OUT "$line\n";
		}
	}
	system "mv temp_files/$filename.temp $filename.curated";
	system "rm -r temp_files";
}