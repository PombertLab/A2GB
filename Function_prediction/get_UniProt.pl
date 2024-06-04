#!/usr/bin/env perl
## Pombert Lab, 2020

my $name = 'get_UniProt.pl';
my $version = '0.4';
my $updated = '2024-06-03';

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

###################################################################################################
## Command line options
###################################################################################################

my $usage = <<"OPTIONS";
NAME        ${name}
VERSION     ${version}
UPDATED     ${updated}
SYNOPSIS    Downloads the SwissProt and/or trEMBL databases from UniProt

RECOMMENDED     Aria2 - https://aria2.github.io/
                # sudo apt install aria2 (Ubuntu/Debian)
                # sudo dnf install aria2 (Fedora/RedHat)

EXAMPLE     ${name} \\
              -s \\
              -t \\
              -f ./ \\
              -n 20 \\
              -l download.log 

OPTIONS:
-s  (--swiss)         Download Swiss-Prot
-t  (--trembl)        Download trEMBL
-f  (--folder)        Download folder [Default: ./]
-l  (--log)           Print download information to log file
-dt (--dtool)         Specify download tool: aria2c, wget or curl ## Tries to autodetect otherwise
-x  (--connex)        Number of aria connections [Default: 10]
-d  (--decompress)    Decompress downloaded files with gunzip     ## trEMBL is huge, off by default
-v  (--version)       Show script version
OPTIONS

unless (@ARGV){
    print "\n$usage\n";
    exit(0);
};

my $swiss;
my $trembl;
my $dtool;
my $folder = './';
my $log;
my $decomp;
my $aria_connections = 10;
my $sc_version;
GetOptions(
    's|swiss' => \$swiss,
    't|trembl' => \$trembl,
    'f|folder=s' => \$folder,
    'l|log=s' => \$log,
    'd|decompress' => \$decomp,
    'x|connex=i' => \$aria_connections,
    'dt|dtool=s' => \$dtool,
    'v|version' => \$sc_version
);

###################################################################################################
## Version
###################################################################################################

if ($sc_version){
    print "\n";
    print "Script:     $name\n";
    print "Version:    $version\n";
    print "Updated:    $updated\n\n";
    exit(0);
}

###################################################################################################
## Aria2 / wget / curl precheck
###################################################################################################

unless ($dtool){

    my @tools = ('aria2c', 'wget', 'curl');

    for my $tool (@tools){

        my $tcheck = `echo \$(command -v $tool)`;
        chomp $tcheck;

        if ($tcheck ne ''){
            $dtool = $tcheck;
            last;
        }

    }

    if (!defined $dtool){

        print "[E] Neither aria2c, wget or curl were found in the \$PATH. ";
        print "Please install one of these.\n";

        print "[E] To install aria2 on Fedora, type: sudo dnf install aria2\n";
        print "Exiting...\n";

        exit(1);

    }
    else{
        print "\n"."Using $dtool as the downloading tool."."\n\n";
    }

}

###################################################################################################
## Precheck db to download
###################################################################################################

if ( (!defined $swiss) && (!defined $trembl) ){
    print "\n";
    print '[E] Neither SwissProt (-s/--swiss) nor trEMBL (-t/--trembl) was requested.'."\n";
    print '[E] Exiting...'."\n\n";
    exit(1);
}

###################################################################################################
## Output directory
###################################################################################################

unless (-d $folder){ 
    mkdir ($folder,0755) or die "Can't create $folder: $!\n";
}

print "\nOutput files will be located in : $folder\n";

if ($log){
    open LOG, ">", "${folder}/${log}";
}

###################################################################################################
## Downloading UniProt SwissProt/TrEMBL databases
###################################################################################################

my $url = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete";

##### SwissProt
if ($swiss){

    if ($log){
        my $date = `date`;
        print LOG "Downloading SwissProt on $date";
    }

    # aria2c
    if ($dtool =~ /aria2c/){

        print "\nDownloading SwissProt with aria2c\n";

        system ("
            aria2c \\
              -x$aria_connections \\
              $url/uniprot_sprot.fasta.gz \\
              --dir=$folder
        ") == 0 or checksig();

    }

    # wget
    elsif ($dtool =~ /wget/){

        print "\nDownloading SwissProt with wget\n";

        system ("
            wget \\
              -P $folder \\
              $url/uniprot_sprot.fasta.gz
        ");

    }

    # curl
    elsif ($dtool =~ /curl/){

        print "\nDownloading SwissProt with curl\n";

        system ("
            curl \\
              -o $folder/uniprot_sprot.fasta.gz \\
              -L $url/uniprot_sprot.fasta.gz
        ") == 0 or checksig();

    }

    if ($log){
        my $size = `du -h $folder/uniprot_sprot.fasta.gz`;
        print LOG "$size";
    }

}

##### TrEMBL
if ($trembl){

    if ($log){
        my $date = `date`;
        print LOG "Downloading trEMBL on $date";
    }

    # aria2c
    if ($dtool =~ /aria2c/){

        print "\nDownloading trEMBL with aria2c. This will take a while...\n\n";

        system ("
            aria2c \\
              -x$aria_connections \\
              $url/uniprot_trembl.fasta.gz \\
              --dir=$folder
        ") == 0 or checksig();

    }

    # wget
    elsif ($dtool =~ /wget/){

        print "\nDownloading trEMBL with wget. This will take a while...\n\n";

        system ("
            wget \\
              -P $folder \\
              $url/uniprot_trembl.fasta.gz
        ");

    }

    # curl
    elsif ($dtool =~ /curl/){

        print "\nDownloading trEMBL with curl. This will take a while...\n\n";

        system ("
            curl \\
              -o $folder/uniprot_trembl.fasta.gz \\
              -L $url/uniprot_trembl.fasta.gz
        ") == 0 or checksig();

    }

    if ($log){
        my $size = `du -h $folder/uniprot_trembl.fasta.gz`;
        print LOG "$size";
    }

}

###################################################################################################
## Decompressing the downloaded databases (if requested)
###################################################################################################

my $gunz = 'gunzip';

## Using unpigz if available
my $unpigz_check = `echo \$(command -v unpigz)`;
chomp $unpigz_check;

if ($unpigz_check ne ''){
    $gunz = 'unpigz';
}

if ($decomp){

    if ($swiss){
        print "\nDecompressing downloaded SwissProt database with $gunz...\n\n";
        system ("
            $gunz \\
            $folder/uniprot_sprot.fasta.gz
        ") == 0 or checksig();
    }

    if ($trembl){
        print "\nDecompressing downloaded trEMBL database with $gunz. Will take a while...\n\n";
        system ("
            $gunz \\
            $folder/trembl.fasta.gz
        ") == 0 or checksig();
    }

}

###################################################################################################
## Completion
###################################################################################################

if ($log){
    my $end = `date`;
    print LOG "Finished downloading database(s) on $end";
}

###################################################################################################
## Subroutine(s)
###################################################################################################

sub checksig {

    my $exit_code = $?;
    my $modulo = $exit_code % 255;

    print "\nExit code = $exit_code; modulo = $modulo \n";

    if ($modulo == 2) {
        print "\nSIGINT detected: Ctrl+C => exiting...\n";
        exit(2);
    }
    elsif ($modulo == 131) {
        print "\nSIGTERM detected: Ctrl+\\ => exiting...\n";
        exit(131);
    }

}