#!/usr/bin/perl
## Pombert Lab, IIT, 2020
my $name = 'parse_annotators.pl';
my $version = '1.6';
my $updated = '2021-06-28';

use strict; use warnings; use Getopt::Long qw(GetOptions);

my $usage = <<"OPTIONS";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	This script parses the output of annotators to help assign putative functions to predicted proteins.
		Annotators are:
		- BLASTP/DIAMOND searches against SwissProt/trEMBL databases
		- InterProScan 5 searches
		- BLASTP/DIAMOND searches against reference organism (optional)

USAGE	${name} \\
		  -q BEOM2.proteins.queries \\
		  -o BEOM2.annotations \\
		  -ip BEOM2.interpro.tsv \\
		  -sl sprot.list \\
		  -sb BEOM2.sprot.blastp.6 \\
		  -tl trembl.list \\
		  -tb BEOM2.trembl.blastp.6 \\
		  -rl reference.list \\
		  -rb reference.blastp.6 ## Searches against reference organism(s) (Optional)

OPTIONS:
-q	List of proteins queried against annotators
-o	Output file
-v	Verbose

## InterProScan5
-ip	TSV output from InterProScan

## BLAST/DIAMOND searches
-sl	SwissProt list of proteins and their products
-sb	SwissProt BLAST/DIAMOND outfmt 6 results
-tl	TREMBL list of proteins and their products
-tb	TREMBL BLAST/DIAMOND outfmt 6 results
-rl	Reference genome(s) list(s) of proteins and their products
-rb	Reference(s) BLAST/DIAMOND outfmt 6 results

## KEGG searches
-ko	KofamKOALA output file
-gk	GhostKOALA output file
-bk	BlastKOALA output file

OPTIONS
die "\n$usage\n" unless @ARGV;

## GetOptions
my $queries;
my $output;
my $verbose;
my $splist; my $spblast;
my $tblist; my $tbblast;
my @rblist; my @rbblast;
my $ipro;
my $kofamkoala;
my $ghostkoala;
my $blastkoala;
GetOptions(
	'q=s' => \$queries,
	'o=s' => \$output,
	'v' => \$verbose,
	'sl=s' => \$splist,
	'sb=s' => \$spblast,
	'tl=s' => \$tblist,
	'tb=s' => \$tbblast,
	'rl=s@{1,}' => \@rblist,
	'rb=s@{1,}' => \@rbblast,
	'ip=s' => \$ipro,
	'ko=s' => \$kofamkoala,
	'gk=s' => \$ghostkoala,
	'bk=s' => \$blastkoala,
);

## Connecting reference feature lists to its corresponding blast file, and storing it in a database to reference during
## parsing step
my %references;
unless (scalar(@rblist) == scalar(@rbblast)){
	die "[E] the number of reference feature lists do not equal the number of reference blast files\n";
}
else {
	for (my $i = 0; $i < scalar(@rblist); $i++){
		$references{$rblist[$i]} = $rbblast[$i];
	}
}

my $time = localtime(); 
my $tstart = time;
print "$time: Obtaining annotations for DIAMOND blastp hits in $spblast and $tbblast...\n";

## Parsing SwissProt blast.6
## Using a double pass for memory optimization and reduce the size of the hash
my %sprot; 
open SB, "<", "$spblast" or die "Can't open $spblast: $!\n"; 
my %sphits;
while (my $line = <SB>){
	chomp $line;
	my @cols = split("\t", $line);
	my $hit = $cols[1];
	$sprot{$hit} = 1;
}
close SB;

open SP, "<", "$splist" or die "Can't open $splist: $!\n";
while (my $line = <SP>){
	chomp $line;
	my @cols = split("\t", $line);
	my $locus = $cols[0];
	my $desc = $cols[1];
	if ( exists $sprot{$locus} ) { $sprot{$locus} = $desc; }
}
close SP; 

open SB, "<", "$spblast" or die "Can't open $spblast: $!\n";
while (my $line = <SB>){
	chomp $line;
	my @cols = split("\t", $line);
	my $query = $cols[0]; 
	my $hit = $cols[1]; 
	my $evalue = $cols[10];
	if ( exists $sphits{$query} ) { next; }
	elsif ( $sprot{$hit} =~ /uncharacterized/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
	elsif ( $sprot{$hit} =~ /hypothetical/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
	elsif ( $sprot{$hit} =~ /predicted protein/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
	else{
		$sphits{$query}[0] = $sprot{$hit};
		$sphits{$query}[1] = $evalue;
	}
}
close SB;

## Parsing TREMBL blast.6
## Using a double pass for memory optimization and reduce the size of the hash
my %trembl;
my %tbhits;
open TBB, "<", "$tbblast" or die "Can't open $tbblast: $!\n"; 
while(my $line = <TBB>){
	chomp $line;
	my @cols = split("\t", $line);
	my $hit = $cols[1];
	$trembl{$hit} = 1;
}
close TBB;

open TB, "<", "$tblist" or die "Can't open $tblist: $!\n";
while(my $line = <TB>){
	chomp $line;
	my @cols = split("\t", $line);
	my $locus = $cols[0];
	my $desc = $cols[1];
	if ( exists $trembl{$locus} ) { $trembl{$locus} = $desc; }
}
close TB;

open TBB, "<", "$tbblast" or die "Can't open $tbblast: $!\n";
while (my $line = <TBB>){
	chomp $line;
	my @cols = split("\t", $line);
	my $query = $cols[0];
	my $hit = $cols[1];
	my $evalue = $cols[10];
	if (exists $tbhits{$query}){ next; }
	elsif ( $trembl{$hit} =~ /uncharacterized/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
	elsif ( $trembl{$hit} =~ /hypothetical/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
	elsif ( $trembl{$hit} =~ /predicted protein/i ){ next; } ## Discarding uninformative BLAST/DIAMOND hits
	else {
		$tbhits{$query}[0] = $trembl{$hit};
		$tbhits{$query}[1] = $evalue;
	}
}
close TBB;

my $time_taken = time - $tstart; 
print "$time: Finished obtaining annotations for $splist and $tblist in $time_taken seconds.\n";

### Parsing InterProScan 5 output
$time = localtime(); 
print "$time: Parsing InterProScan5 $ipro...\n";
$time = time;
open IP, "<", "$ipro" or die "Can't open $ipro: $!\n";
my %pfam = ();
my %hamap = ();
my %tigr = ();
my %cdd = ();
while (my $line = <IP>){
	chomp $line;
	my @cols = split("\t", $line);
	## Columns info from https://github.com/ebi-pf-team/interproscan/wiki/OutputFormats
	## $cols[0] Protein Accession (e.g. P51587)
	## $cols[1] Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
	## $cols[2] Sequence Length (e.g. 3418)
	## $cols[3] Analysis (e.g. Pfam / PRINTS / Gene3D)
	## $cols[4] Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
	## $cols[5] Signature Description (e.g. BRCA2 repeat profile)
	## $cols[6] Start location
	## $cols[7] Stop location
	## $cols[8] Score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
	## $cols[9] Status - is the status of the match (T: true)
	## $cols[10] Date - is the date of the run
	## $cols[11] (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
	## $cols[12] (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
	## $cols[13] (GO annotations (e.g. GO:0005515) - optional column; only displayed if --goterms option is switched on)
	## $cols[14] (Pathways annotations (e.g. REACT_71) - optional column; only displayed if --pathways option is switched on)

	if ($cols[3] eq 'Pfam'){
		if (exists $pfam{$cols[0]}){
			if ($cols[8] < $pfam{$cols[0]}[1]){
				$pfam{$cols[0]}[0] = $cols[5];
				$pfam{$cols[0]}[1] = $cols[8];
			}
		}
		else {
			$pfam{$cols[0]}[0] = $cols[5];
			$pfam{$cols[0]}[1] = $cols[8];
		}
	}
	elsif ($cols[3] eq 'TIGRFAM'){
		if (exists $tigr{$cols[0]}){
			if ($cols[8] < $tigr{$cols[0]}[1]){
				$tigr{$cols[0]}[0] = $cols[5]; 
				$tigr{$cols[0]}[1] = $cols[8];
			}
		}
		else {
			$tigr{$cols[0]}[0] = $cols[5];
			$tigr{$cols[0]}[1] = $cols[8];
		}
	}
	elsif ($cols[3] eq 'Hamap'){
		if (exists $hamap{$cols[0]}){
			if ($cols[8] > $hamap{$cols[0]}[1]){
				$hamap{$cols[0]}[0] = $cols[5];
				$hamap{$cols[0]}[1] = $cols[8];
			}
		}
		else {
			$hamap{$cols[0]}[0] = $cols[5];
			$hamap{$cols[0]}[1] = $cols[8];
		}
	}
	elsif ($cols[3] eq 'CDD'){
		if (exists $cdd{$cols[0]}){
			if ($cols[8] > $cdd{$cols[0]}[1]){
				$cdd{$cols[0]}[0] = $cols[5];
				$cdd{$cols[0]}[1] = $cols[8];
			}
		}
		else {
			$cdd{$cols[0]}[0] = $cols[5]; 
			$cdd{$cols[0]}[1] = $cols[8];
		}
	}
}
close IP;
$time_taken = time - $time;
$time = localtime();
print "$time: Finished parsing InterProScan5 $ipro in $time_taken seconds...\n";

## KofamKOALA results, if any
my %kofam;
open KOFAM, "<", "$kofamkoala" or die "Can't open $kofamkoala: $!\n";
while (my $line = <KOFAM>){
	chomp $line;
	if ($line =~ /^#/) { next; }
	elsif ($line =~ /^\* /){
		my @columns = split (/\s+/, $line);
		my $query = $columns[1];
		my $ko = $columns[2];
		my $evalue = $columns[5];
		my $definition;
		for my $num (6..$#columns){
			$definition .= "$columns[$num] ";
		}
		$definition =~ s/ $//;

		## Adding to %kofam db
		$kofam{$query}{'ko'} = $ko;
		$kofam{$query}{'evalue'} = $evalue;
		$kofam{$query}{'def'} = $definition;
	}
}
close KOFAM;

## GhostKOALA results, if any
my %ghost;
open GHOST, "<", "$ghostkoala" or die "Can't open $ghostkoala: $!\n";
while (my $line = <GHOST>){
	chomp $line;
	if ($line =~ /^#/) { next; }
	else {
		my @columns = split ("\t", $line);
		my $query = $columns[0];
		my $ko = $columns[1];
		my $definition = $columns[2];
		my $score = $columns[3];

		## Adding to %ghost db
		if ($ko){
			$ghost{$query}{'ko'} = $ko;
			$ghost{$query}{'score'} = $score;
			$ghost{$query}{'def'} = $definition;
		}
	}
}
close GHOST;

## BlastKOALA results, if any
my %blastko;
open BLASTKO, "<", "$blastkoala" or die "Can't open $blastkoala: $!\n";
while (my $line = <BLASTKO>){
	chomp $line;
	if ($line =~ /^#/) { next; }
	else {
		my @columns = split ("\t", $line);
		my $query = $columns[0];
		my $ko = $columns[1];
		my $definition = $columns[2];
		my $score = $columns[3];

		## Adding to %blastko db
		if ($ko){
			$blastko{$query}{'ko'} = $ko;
			$blastko{$query}{'score'} = $score;
			$blastko{$query}{'def'} = $definition;
		}
	}
}
close BLASTKO;

### Reference organism, if any
my %refhits;
if ($rblist[0]){
	foreach my $ref (keys(%references)){
		## Parsing product list
		my $time = localtime();
		print "$time: Parsing the product list $ref...\n";
		my $tstart = time;
		my %reforg;
		open REFL, "<", "$ref" or die "Can't open $ref: $!\n"; 
		while (my $line = <REFL>){
			chomp $line; my @cols = split("\t", $line); 
			$reforg{$cols[0]} = $cols[1];
		}
		my $time_taken = time - $tstart;
		$time = localtime(); 
		print "$time: Finished parsing reference product list $ref in $time_taken seconds.\n";
		close REFL;
		## Parsing BLAST/DIAMOND outfmt 6 file from the reference organism
		$time = localtime(); 
		print "$time: Parsing DIAMOND blastp file $references{$ref}...\n";
		$tstart = time;
		open REFB, "<", "$references{$ref}" or die "Can't open $references{$ref}: $!\n";
		while (my $line = <REFB>){
			chomp $line;
			my @cols = split("\t", $line);
			my $query = $cols[0];
			my $hit = $cols[1];
			my $evalue = $cols[10];
			if ( exists $refhits{$query} ) { next; }
			elsif (!exists $reforg{$hit}) { ## Checking if hit is missing from reference file
				if ($verbose) { print "Description for $hit (BLAST hit) not found in $ref. Skipping...\n"; }
				next;
			}
			elsif ( $reforg{$hit} =~ /uncharacterized/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
			elsif ( $reforg{$hit} =~ /hypothetical/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
			elsif ( $reforg{$hit} =~ /predicted protein/i ) { next; } ## Discarding uninformative BLAST/DIAMOND hits
			else {
				$refhits{$ref}{$query}[0] = $reforg{$hit}; 
				$refhits{$ref}{$query}[1] = $evalue;
			}
		}
		close REFB;
		$time = localtime();
		$time_taken = time - $tstart;
		print "$time: Finished parsing DIAMOND blastp file $references{$ref} in $time_taken seconds...\n";
	}
}

## Creating the parsed annotations file
open QUE, "<", "$queries" or die "Can't open $queries: $!\n";
$time = localtime(); 
print "$time: Writing annotations to $output...\n";
$time = time;
open OUT, ">", "$output" or die "Can't create $output: $!\n";

### Writing header
print OUT '#'."Locus_tag"."\tEvalue\tSwissProt"."\tEvalue\ttrEMBL"."\tEvalue\tPFAM"."\tEvalue\tTIGR"."\tScore\tHAMAP"."\tEvalue\tCDD";
if ($kofamkoala){ print OUT "\tEvalue\tKofamKOALA"; }
if ($ghostkoala){ print OUT "\tEvalue\tGhostKOALA"; }
if ($blastkoala){ print OUT "\tEvalue\tBlastKOALA"; }
if ($rblist[0]){
	foreach my $ref (sort(keys(%references))){
		$ref =~ s/\.list//;
		print OUT "\tEvalue\t$ref";
	}
}
print OUT "\n";

### Writing data
while (my $line = <QUE>){
	chomp $line;
	print OUT "$line\t";
	if ( exists $sphits{$line} ) { print OUT "$sphits{$line}[1]\t$sphits{$line}[0]\t"; }
	else { print OUT "NA\thypothetical protein\t"; }
	if ( exists $tbhits{$line}) { print OUT "$tbhits{$line}[1]\t$tbhits{$line}[0]\t"; }
	else { print OUT "NA\thypothetical protein\t"; }
	if ( exists $pfam{$line}) { print OUT "$pfam{$line}[1]\t$pfam{$line}[0]\t"; }
	else { print OUT "NA\thypothetical protein\t"; }
	if ( exists $tigr{$line}) { print OUT "$tigr{$line}[1]\t$tigr{$line}[0]\t"; }
	else { print OUT "NA\thypothetical protein\t"; }
	if ( exists $hamap{$line}) { print OUT "$hamap{$line}[1]\t$hamap{$line}[0]\t"; }
	else { print OUT "NA\thypothetical protein\t"; }
	if ( exists $cdd{$line} ){ print OUT "$cdd{$line}[1]\t$cdd{$line}[0]"; }
	else { print OUT "NA\tno motif found"; }

	## Printing KofamKOALA results, if any
	if ($kofamkoala){
		if (exists $kofam{$line}){
			print OUT "\t$kofam{$line}{'evalue'}\t$kofam{$line}{'ko'}; $kofam{$line}{'def'}";
		}
		else { print OUT "\tNA\tno match found"; }
	}

	## Printing GhostKOALA results, if any
	if ($ghostkoala){
		if (exists $ghost{$line}){
			print OUT "\t$ghost{$line}{'score'}\t$ghost{$line}{'ko'}; $ghost{$line}{'def'}";
		}
		else { print OUT "\tNA\tno match found"; }
	}

	## Printing blastKOALA results, if any
	if ($blastkoala){
		if (exists $blastko{$line}){
			print OUT "\t$blastko{$line}{'score'}\t$blastko{$line}{'ko'}; $blastko{$line}{'def'}";
		}
		else { print OUT "\tNA\tno match found"; }
	}

	## Printing reference(s) results, if any
	if ($rblist[0]){
		foreach my $ref (sort(keys(%references))){
			if ( exists $refhits{$ref}{$line} ) { print OUT "\t$refhits{$ref}{$line}[1]\t$refhits{$ref}{$line}[0]"; }
			else { print OUT "\tNA\tno match found"; }
		}
	}

	print OUT "\n";
}
close QUE; 
close OUT;
$time_taken = time - $time;
my $ttime = time - $tstart;
$time = localtime(); 
print "$time: Done writing annotations in $time_taken seconds...\n";
$time = localtime();
print "$time: Task completed in $ttime seconds. Exiting...\n";
exit();