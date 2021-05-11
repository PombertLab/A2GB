#!/usr/bin/perl

use strict; use warnings; use Getopt::Long qw(GetOptions); use File::Basename;

my $name = 'Create_ChimeraX_Session.pl';
my $version = '0.1a';
my $updated = '2021-05-11';

my $usage = <<"EXIT";
NAME		${name}
VERSION		${version}
UPDATED		${updated}

SYNOPSIS	The purpose of this script is to iterate thorugh all GESAMT match predictions,
		run Create_ChimeraX_Session.pl to generate ChimeraX sessions for easy 3D comparison
		during annotation creation.

COMMAND		${name} \
		-gr .../GESAMT/*.gesamt \
		-r .../RCSB/PDB/

OPTIONS
-i (--pred_pdb)		Path to predicted .pdb files
-gr (--gesamt_rcsb)	RCSB GESAMT match file
-gp (--gesamt_pfam)	PFAM GESAMT match file
-r (--rcsb_pdb)		Path to RCSB PDB directories
-p (--pfam_pdb)		Path to PFAM .pdb files
-o (--outdir)		Output directory [Default: ./CXS]
EXIT

die("\n$usage\n") unless(@ARGV);

my $pdb;
my $g_rcsb;
my $g_pfam;
my $rcsb;
my $pfam;
my $outdir = './CXS';
GetOptions(
	"i|pred_pdb=s" => \$pdb,
	"gr|gesamt_rcsb=s" => \$g_rcsb,
	"gp|gesamt_pfam=s" => \$g_pfam,
	"r|rcsb_pdb=s" => \$rcsb,
	"p|pfam_pdb=s" => \$pfam,
	"o|outdir=s" => \$outdir
);

unless(-d $outdir){
	mkdir($outdir,0775) or die("Can't create $outdir: $!\n");
}

## Load predicted pdb filenames into database
my %pred;
opendir(DIR,$pdb) or die("Can't open $pdb: $!\n");
while (my $file = readdir(DIR)){
	if ($file =~ /^(\w+)/){
		$pred{$1} = "$pdb/$file";
	}
}
## Load RCSB file locations into database
my %db;
if ($rcsb){
	opendir(DIR,$rcsb) or die("Can't open $rcsb: $!\n");
	while (my $dir = readdir(DIR)){
		opendir(DIR2,"$rcsb/$dir") or die("Can't open $rcsb/$dir: $!\n");
		while (my $file = readdir(DIR2)){
			if ($file =~ /^pdb/){
				$db{$file} = "$rcsb/$dir/$file";
			}
		}
		closedir DIR2;
	}
	closedir DIR;
}

## Load PFAM file locations into database
if ($pfam){
	opendir(DIR,$pfam) or die("Can't open $pfam: $!\n");
	while (my $file = readdir(DIR)){
		if($file =~ /^PF/){
			$db{$file} = "$pfam/$file";
		}
	}
	closedir DIR;
}

## For each match, get the filename of the best match for PFAM
my %sessions; ## ChimeraX session to be created
if($g_pfam){
	open IN, "<", "$g_pfam";
	my $locus_tag;
	my $match_found = 0;
	while (my $line = <IN>){
		if ($line =~ /^###\s(.*)+$/){
			$locus_tag = $1;
			$match_found = 1;
			next;
		}
		elsif (($match_found == 1) && ($line =~ /(\S+.pdb)\t.+?$/)){
			push(@{$sessions{$locus_tag}},$db{$1});
			$match_found = 0;
		}
	}
}

## For each match, get the filename of the best match for RCSB
if ($g_rcsb){
	open IN, "<", "$g_rcsb" or die("Can't open $g_rcsb: $!\n");
	my $locus_tag;
	my $match_found = 0;
	while (my $line = <IN>){
		chomp($line);
		if ($line =~ /^###\s(.*)$/){
			$locus_tag = $1;
			$match_found = 1;
			next;
		}
		elsif (($match_found == 1) && ($line =~ /(\S+\.ent\S*)\t.+?$/)){
			my $file = $1; $file =~ s/.gz//;
			push(@{$sessions{$locus_tag}},$db{$file});
			$match_found = 0;
		}
	}
}

my ($filename,$dir) = fileparse($0);
my $script = $dir."ChimeraX_helper_scripts/Create_ChimeraX_Session.py";
foreach my $locus (keys(%sessions)){
	if(scalar(@{$sessions{$locus}}) == 2){
		my $pred_pdb_file = $pred{$locus};
		my $pfam_pdb_file = $sessions{$locus}[0];
		my $rcsb_pdb_file = $sessions{$locus}[1];
		system("chimerax --nogui $script \\
			-p_f $pred_pdb_file \\
			-p_m $pfam_pdb_file \\
			-r_m $rcsb_pdb_file \\
			-o_d $outdir"
		);
	}
	elsif(scalar(@{$sessions{$locus}}) == 1){
		my $pred_pdb_file = $pred{$locus};
		my $ref_pdb_file = $sessions{$locus}[0];
		system("chimerax --nogui $script \\
			-p_f $pred_pdb_file \\
			-p_m $ref_pdb_file \\
			-o_d $outdir"
		);
	}
}