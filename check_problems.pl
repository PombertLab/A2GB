#!/usr/bin/perl
## Pombert Lab, IIT, 2021
use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

my $name = "check_problems.pl";
my $version = "0.4";
my $updated = '2021-04-05';

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}
SYNOPSIS	This script is used to notify the user of protein abnormalities, such as a missing 
		start methionines and internal stop codons. This script can also run EMBLtoFeatures.pl to 
		update the .prot files before checking for these abnormalities.

COMMAND		${name} \\
		  -p *.prot \\
		  -o ProteinCheck.log \\
		  -u \\
		  -v

OPTIONS
-p | --prot	FASTA files (.prot) to be checked for abnormalities
-o | --out	Print the output to a log file
-u | --update	Update .prot files using EMBLtoFeatures.pl
-v | --verb	Add verbosity
EXIT
die "\n$usage\n\n" unless @ARGV;

my @prot_files;
my $out;
my $update;
my $verb;
GetOptions(
	'p|prot=s@{1,}' => \@prot_files,
	'o|out=s' => \$out,
	'u|update' => \$update,
	'v|verb' => \$verb
);

if ($out){
	open OUT, ">", "$out" or die "Can't create $out: $!\n";
}

## Runs EMBLtoFeatures.pl if flag update is on
if ($update){
	foreach my $file (@prot_files){
		my ($filename,$dir) = fileparse($file);
		my $name = basename($file,".prot");
		if (-f "$name.embl"){
			if ($verb){
				system "EMBLtoFeatures.pl \\
				  -e $name.embl \\
				  -v";
			}
			else{
				system "EMBLtoFeatures.pl \\
				  -e $name.embl";
			}
		}
		else {
			print "[W] $file has no correlating .fsa file\n";
		}
	}
}

## Checking for problems in the .prot files 
for my $prot_file (@prot_files){

	open IN, "<", "$prot_file" or die "Can't open $prot_file: $!\n";
	my ($filename,$dir) = fileparse($prot_file);

	my %sequences;
	my $locus;
	while (my $line = <IN>){
		chomp $line;
		if($line =~ /^>(\S+)/){
			$locus = $1;
			next;
		}
		$sequences{$locus} .= $line;
	}

	if ($verb) { print "\nChecking for sequence errors in $prot_file located in $dir\n"; }
	my $count = undef;

	## Iterating through locus tags in database %sequences
	foreach my $locus_tag (sort(keys %sequences)){
		my $line = $sequences{$locus_tag};
		unless ($count){
			if (($line !~ /^M/) || ($line =~ /\W/)){
				if ($out){
					print OUT "\n\t\tInvalid Start Codon\tInternal Stop Codon\n\n";
				}
				if ($verb){
					print "\n\t\tInvalid Start Codon\tInternal Stop Codon\n\n";
				}
			}
		}

		## No start methionine (bad) + internal stop codon found (bad)
		if (($line !~ /^(M)/) && ($line =~ /\W/)){
			my ($aa) = $line =~ /^(\w)/;
			if ($out){
				print OUT "$locus_tag\t\t$aa\t\t\tX\n";
			}
			if ($verb){
				print "$locus_tag\t\t$aa\t\t\tX\n";
			}
			$count = 1;
		}

		## No start methionine (bad), but no internal stop codon (good)
		elsif ($line !~ /^(M)/){
			my ($aa) = $line =~ /^(\w)/;
			if ($out){
				print OUT "$locus_tag\t\t$aa\t\t\t.\n";
			}
			if ($verb){
				print "$locus_tag\t\t$aa\t\t\t.\n";
			}
			$count = 1;
		}

		## Start methionine found (good), internal stop codon found (bad)
		elsif ($line =~ /\W/){
			if ($out){
				print OUT "$locus_tag\t\t.\t\t\tX\n";
			}
			if ($verb){
				print "$locus_tag\t\t.\t\t\tX\n";
			}
			$count = 1;
		}
	}
	
	## No error found
	unless ($count){ 
		if ($out){
			print OUT "\nOK: No error found in $prot_file\n";
		}
		if ($verb){
			print "\nOK: No error found in $prot_file\n";
		}
	}
}
print "\n";
