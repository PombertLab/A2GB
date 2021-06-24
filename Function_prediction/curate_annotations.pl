#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'curate_annotations.pl';
my $version = '1.9a';
my $updated = '2021-06-24';

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	Displays lists of functions predicted per proteins. User can select or enter desired annotation.
		Creates a tab-delimited .curated list of annotations.

COMMAND		${name} \\
		  -sq proteins.annotations \\
		  -r

OPTIONS
-sq (--seq_hom)		Sequence homology based annotations (generated from parse_annotators.pl)
-rd (--rcsb_3d)		3D structural homology based annotations (Generated with descriptive_GESAMT_matches.pl)
-pd (--pfam_3d)		3D structural homology annotations based on predicted stuctures (Generated with descriptive)
-cx (--chimerax)	Path to ChimeraX pdb sessions
-r (--resume)		Resume annotation from last curated locus_tag
-c (--check)		Check loci marked with '?'
-v (--verify)		Check loci marked for 3D verification
EXIT
die "\n\n$usage\n\n" unless @ARGV;

my $input;
my $rcsb_input;
my $pfam_input;
my $chimerax;
my $continue;
my $review;
my $verify;
GetOptions(
	"sq|sequence=s" => \$input,
	"rd|rcsb_3d=s" => \$rcsb_input,
	"pd|pfam_3d=s" => \$pfam_input,
	"cx|chimerax=s" => \$chimerax,
	"r|resume" => \$continue,
	"v|verify" => \$verify,
	"c|check" => \$review,
);

my ($h_filename,$h_dir) = fileparse($0);
my $script = $h_dir."ChimeraX_helper_scripts/Restore_ChimeraX_Session.py";

## Make a temporary directory to store our temporary files.
unless (-d 'temp_files'){
	mkdir ('temp_files',0755) or die "Can't make directory temp_files: $!\n";
}

my ($filename,$dir) = fileparse($input);
open OUT, ">", "temp_files/$filename.temp" or die "Can't create temp_files/$filename.temp: $!\n";

## If reviewing annotations, load existing annotations into RAM
my $last_locus = undef;
my @to_review;
if ($review||$verify){
	print "verifying\n";
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
				print "\n\nWould you like to [c]ontinue, [r]estart, or [e]xit?: ";
				chomp (my $proceed = lc(<STDIN>));
				if ( $proceed eq 'r' ) {
					system "clear";
					system "mv $filename.curated $filename.backup.curated";
					last WHILE;
				}
				elsif ( $proceed eq 'e' ) {
					print "\nScript terminating...\n\n";
					exit;
					last WHILE;
				}
				elsif ($proceed eq 'c') {
					last_locus();
					print "continuing\n";
					print "$last_locus\n";
					last WHILE;
				}
				else {
					system "clear";
					print "\nERROR: Invalid input '$proceed'\n\n" 
				}
			}
		}
		else{
			last_locus();
		}
	}
}

## If we are given a rcsb 3D file, parse through it to get the match data
my %rcsb; my $rcsb_query;
if ($rcsb_input){

	open R3D, "<", "$rcsb_input" or die "\n[E] Can't open 3D matches file: $rcsb_input\n\n";
	
	while (my $line = <R3D>){
		chomp $line;
		if ($line =~ /^### .*?(\w+)\-\w+\-\w+$/){ $rcsb_query = $1; }
		else {
			my @data = split("\t",$line);
			my $hit = $data[2];
			my $chain = $data[3];
			my $qscore = $data[4];
			my $match = $data[$#data];
			push (@{$rcsb{$rcsb_query}}, "$qscore\t$hit - Chain: $chain - $match");
		}
	}
}

## If we are given a pfam 3D file, parse through it to get the match data
my %pfam; my $pfam_query;
if ($pfam_input){

	open P3D, "<", "$pfam_input" or die "\n[E] Can't open 3D matches file: $pfam_input\n\n";
	
	while (my $line = <P3D>){
		chomp $line;
		if ($line =~ /^### .*?(\w+)\-\w+\-\w+$/){ $pfam_query = $1; }
		else {
			my @data = split("\t",$line);
			my $hit = $data[$#data-1];
			my $qscore = $data[3];
			my $match = $data[$#data];
			$hit =~ s/.pdb$//;
			push (@{$pfam{$pfam_query}}, "$qscore\t$hit - $match");
			
		}
	}
}

my %cxs;
if ($chimerax){
	opendir(DIR,$chimerax);
	while (my $file = readdir(DIR)){
		if($file =~ /(\S+)\.cxs$/){
			$cxs{$1} = "$chimerax/$file";
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

my $last_annotation;
my $previous_annotation;
while (my $line = <IN>){
	
	system "clear";
	chomp $line;
	
	## If we are at the header of the file, get the header info and skip to the next line
	if ($line =~ /^#/){
		$line =~ s/#//;
		@header_info = split("\t",$line);
		next;
	}
	$current_protein++;
	
	## If we are reviewing, skip the check start process because we are not going to be curating the same way.
	my $annon_notes;
	unless ($review||$verify){
		## Check to see if we are the $last_locus, and if we are, start the anntotation process
		unless ($start){
			if ($line =~ /$last_locus/){
				$start = 1;
			}
			next;
		}
	}
	else {
		## If there are still more loci to review, check for review necessity
		if (@to_review){
			my $review_line = shift(@to_review);
			## Check to see if we are at an annotation that needs review, indicated by a '?'
			## If we are reading a line that does not require reviewing or verification, print the line out, and move to the next one.
			## If we are reading a line that does require reviewing, print the normal annotation curation printouts with annotaiton notes.
			## If we are reading a line that does require verification, print the normal annotation printout
			unless ($review_line =~ /\?|Verify 3D Structural Homology/){
				print OUT "$review_line\n";
				next;
			}
			elsif ($review_line =~ /\?/ && !$review) {
				print OUT "$review_line\n";
			}
			elsif ($review_line =~ /Verify 3D Structural Homology/ && !$verify){
				print OUT "$review_line\n";
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

	## Add rcsb information if provided
	my $rcsb_sources = 0;
	my $rcsb_predictions = 0;
	if (defined $rcsb_input){
		if (defined $rcsb{$locus}){
			foreach my $struct (@{$rcsb{$locus}}){
				push(@sources,"GESAMT");
				my @data = split ("\t", $struct);
				push(@predictions,$data[1]);
				push(@evalues,$data[0]);
				$rcsb_sources++;
				$rcsb_predictions++;
			}
		}
		else {
			push(@sources,"GESAMT");
			push(@predictions,"no match found");
			push(@evalues,"NA");
			$NA_count++;
			$rcsb_sources++;
		}
	}

	## Add pfam information if provided
	my $pfam_sources = 0;
	my $pfam_predictions = 0;
	if (defined $pfam_input){
		if (defined $pfam{$locus}){
			foreach my $struct (@{$pfam{$locus}}){
				push(@sources,"GESAMT");
				my @data = split ("\t", $struct);
				push(@predictions,$data[1]);
				push(@evalues,$data[0]);
				$pfam_sources++;
				$pfam_predictions++;
			}
		}
		else {
			push(@sources,"GESAMT");
			push(@predictions,"no match found");
			push(@evalues,"NA");
			$NA_count++;
			$pfam_sources++;
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
	my $choices = $options - $rcsb_sources - $pfam_sources;
	WHILE: while (0==0){
		
		print "\n$status\t$current_protein/$protein_total\n";
		
		if($review){
			print "\n## Annotation Notes:\n";
			if($annon_notes) {
				print "\t$annon_notes\n";
			}
			else{
				print "\tN/A\n";
			}
		}

		print "\n## Putative annotation(s) found for protein $data[0]:\n";
		
		## Loop through our information arrays so we don't have a million conditionals, and we are more robust this way
		my $rcsb_status = 0;
		my $pfam_status = 0;
		for (my $i = 1; $i <= $options; $i++){
			my $source = $sources[$i-1];
			my $prediction = $predictions[$i-1];
			my $evalue = $evalues[$i-1];
			$string_length = length("$source"); tab();
			if ($i <= $choices){
				print "\t${i}.\t${source}${tab}";
				$string_length = length("$evalue"); tab();
				print "${evalue}${tab}${prediction}\n";
			}
			elsif($rcsb_status < $rcsb_sources) {
				if ($rcsb_status == 0){
					print "\n".'## 3D structural homologs based on experimentally determined structures (if any):'."\n";
				}
				$rcsb_status ++;
				print "\tRCSB.\t${source}${tab}";
				$string_length = length("$evalue"); tab();
				print "${evalue}${tab}${prediction}\n";
			}
			else{
				if ($pfam_status == 0){
					print "\n".'## 3D structural homologs based on predicted structures (if any):'."\n";
				}
				$pfam_status++;
				print "\tPFAM.\t${source}${tab}";
				$string_length = length("$evalue"); tab();
				print "${evalue}${tab}${prediction}\n";
			}
		}

		print "\nPlease enter:\n\n";
		print "\t[1-$choices] to assign annotation\n";
		print "\t[0] to annotate the locus as a 'hypothetical protein'\n";
		print "\t[m] to manually annotate the locus, e.g. DUFxxx domain-containing protein\n";
		print "\t[n] to manually annotate the locus with annotation notes, e.g. structural homolog\n";
		if ($annon_notes) { print "\t[k] to keep annotation\n"; }
		if ($rcsb_predictions > 0) { print "\t[v] to mark this annotation for 3D structural verification\n"; }
		print "\t[?] to mark this annotation for review and add annotation notes (optional)\n\n";
		if ($rcsb_predictions > 0 || $pfam_predictions > 0){ print "\t[d] to display 3D homology for $locus\n\n"; }
		print "\t[x] to exit.\n";
		print "\nSelection: ";
		chomp (my $select = <STDIN>);

		## Select exit curation
		if ($select eq 'x'){
			if ($review){
				print OUT "$locus\t?\t$annon_notes\n";
			}
			elsif ($verify){
				print OUT "$locus\tVerify 3D Structural Homology\t";
				foreach my $struct (@{$rcsb{$locus}}){
					if ($struct =~ /^\S+\t(\S+)/) {
						print OUT "$1,";
					}
				}
				print OUT "\n";
			}
			cleanup();
			print "\nExiting annotation curation...\n\n";
			exit;
		}
		## Enter manual curation (no notes)
		elsif ($select eq 'm'){
			print "Enter desired annotation: ";
			chomp (my $manual = <STDIN>);
			if ($manual !~ /^\S+/ && length($manual) < 5){
				system "clear";
				print "\nERROR: Annotation cannot start with ' '.\n\n";
			}
			else{
				print OUT "$locus\t$manual\n";
				system "clear";
				last WHILE;
			}
		}
		## Enter manual curation (with notes)
		elsif ($select eq 'n'){
			print "Enter desired annotation: ";
			chomp (my $manual = <STDIN>);
			if ($manual !~ /^\S+/ && length($manual) < 5){
				system "clear";
				print "\nERROR: Annotation cannot start with ' '.\n\n";
			}
			else{
				print "Desired Annotation Note: ";
				chomp (my $note = <STDIN>);
				print OUT "$locus\t$manual\tNote: $note\n";
				system "clear";
				last WHILE;
			}
		}
		## Keep annotation
		elsif ($select eq 'k' && $annon_notes){
			print OUT "$locus\t? $annon_notes\n";
			last WHILE;
		}
		## Mark for 3D verification
		elsif ($select eq 'v' && ($rcsb_predictions > 0 || $pfam_predictions > 0)){
			print OUT "$locus\tVerify 3D Structural Homology\t";
			foreach my $struct (@{$rcsb{$locus}}){
				if ($struct =~ /^\S+\t(\S+)/) {
					print OUT "$1,";
				}
			}
			print OUT "\n";
			system "clear";
			last WHILE;
		}
		## Mark locus for review
		elsif ($select eq '?'){
			print "\nAnnotation note(s) (if any): ";
			chomp(my $note = <STDIN>);
			print OUT "$locus\t? $note\n";
			system "clear";
			last WHILE;
		}
		## Annotate locus as hypothetical protein
		elsif ($select eq '0'){
			print OUT "$locus\thypothetical protein\n";
			system "clear";
			last WHILE;
		}
		## Select annotation from sequence homology
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
		elsif ($select eq "d" && ($rcsb_predictions > 0 || $pfam_predictions > 0)){
			if (exists($cxs{$locus})){
				system "clear";
				system "chimerax $cxs{$locus} $script &";
			}
		}
		## If option isn't one of the following, throw an error and let the user know
		else {
			system "clear";
			print "\nERROR: Invalid input value '$select'.\n\n";
		}
	}
	system "pkill chimerax";
}
cleanup();
print "Annotation curation completed. Be sure to run with the -c or -v flag to check any non-annotated loci...\n";

##subroutines

sub tab {
	if ($string_length >= 8){ $tab = "\t"; }
	else { $tab = "\t\t"; }
}

sub cleanup {
	if ($review||$verify){
		while (my $line = shift(@to_review)){
			print OUT "$line\n";
		}
	}
	system "mv temp_files/$filename.temp $filename.curated";
	system "rm -r temp_files";
}

sub last_locus {
	open IN, "<", "$filename.curated" or die "Can't read $filename.curated: $!\n";
		while (my $line = <IN>){
			chomp $line;
			print OUT "$line\n";
			my @data = split("\t",$line);
			$last_locus = $data[0];
		}
	close IN;
}
