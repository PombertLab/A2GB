# A2GB – Annotations to GenBank

[![DOI](https://zenodo.org/badge/300303510.svg)](https://zenodo.org/doi/10.5281/zenodo.5532811)

[A2GB](https://github.com/PombertLab/A2GB/) is a eukaryote genome annotation pipeline that will transform the annotation files exported from [Apollo](https://genomearchitect.readthedocs.io/en/latest/) (fomerly known as [WebApollo](http://gmod.org/wiki/WebApollo)) in preparation for sequence submission to NCBI’s [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). While its primary goal is to submit annotations to GenBank for the generation of accession numbers, the conversion of file formats will permit the use of different tools in downstream analyses. 

In a stepwise approach, the A2GB pipeline converts annotation files from GFF3 -> EMBL -> TBL -> ASN. Each format will be useful for diagnostic quality checks of the annotations or become the springboard for other analyses, such as protein function prediction.  

Furthermore, A2GB acts as a guide to prepare sequence submissions according to [NCBI’s guidelines](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/), including project registration with [BioSample](https://www.ncbi.nlm.nih.gov/biosample) and [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) for the generation of [locus_tag](https://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf) prefixes. Upon acceptance from NCBI, these essential steps will ensure that sequence data will be made publicly available through GenBank and other member databases within the [International Nucleotide Sequence Database Collaboration](http://www.insdc.org/).

## Table of contents
* [Introduction](#introduction)
* [Requirements](#requirements)
* [A2GB workflow](#A2GB-workflow)
   * [Exporting annotations from Apollo](#Exporting-annotations-from-Apollo)
        *	[Splitting Apollo GFF3 files](#Splitting-Apollo-GFF3-files)
   *	[Converting GFF3 files to EMBL format](#Converting-GFF3-files-to-EMBL-format)
        *	[Checking for internal stop codons and missing start methionines](#Checking-for-internal-stop-codons-and-missing-start-methionines)
        *	[Creating tab-delimited lists of RNA locus tags and their products](#Creating-tab-delimited-lists-of-RNA-locus-tags-and-their-products)
   *	[Protein function prediction](#Protein-function-prediction)
        *	[Predicting functions with InterProScan5](#Predicting-functions-with-InterProScan5)
        *	[Performing homology searches against UniProt databases](#Performing-homology-searches-against-UniProt-databases)
	        *	[Downloading the SwissProt and TrEMBL databases](#Downloading-the-SwissProt-and-TrEMBL-databases)
       		*	[Creating tab-delimited product lists from UniProt databases](#Creating-tab-delimited-product-lists-from-UniProt-databases)
      		*	[Running DIAMOND or BLAST searches against UniProt databases](#Running-DIAMOND-or-BLAST-searches-against-UniProt-databases)
        *	[Performing homology searches against reference datasets](#Performing-homology-searches-against-reference-datasets)
        *	[Searching KEGG databases for orthologs](#Searching-KEGG-databases-for-orthologs)
        *	[Parsing the results from homology searches](#Parsing-the-results-from-homology-searches)
        *	[Curating the protein annotations](#Curating-the-protein-annotations)
   *	[Converting EMBL files to ASN format](#Converting-EMBL-files-to-ASN-format)
        *	[Converting EMBL files to TBL format](#Converting-EMBL-files-to-TBL-format)
        *	[Converting TBL files to ASN format](#Converting-TBL-files-to-ASN-format)
	        *	[Adding metadata to FASTA files](#Adding-metadata-to-FASTA-files)
      		*	[Generating a GenBank submission template](#Generating-a-GenBank-submission-template)
        	*	[Creating structured comments](#Creating-structured-comments)
        	*	[Using table2asn](#Using-table2asn)
        *	[Checking for errors and fixing them](#Checking-for-errors-and-fixing-them)
        	*	[Partial genes](#Partial-genes)
        	*	[Missing stop codons and GT-AG splice sites](#Missing-stop-codons-and-GT-AG-splice-sites)
        	*	[Fixing errors](#Fixing-errors)	
   *	[Submitting ASN files to GenBank](#Submitting-ASN-files-to-GenBank)
* [Funding and acknowledgments](#Funding-and-acknowledgments)
* [References](#References)

## Introduction

Various analyses contribute to assigning biological interpretation to a DNA sequence. The goal of genome annotation is to identify the location and function of a genome's encoded features, particularly of protein-coding and RNA-coding genes. Thus, generating the best descriptions of encoded features (*i.e.* annotations) is paramount for any genomic sequencing project that intends to identify these encoded features and attribute biological meaning to them.

[Apollo](https://genomearchitect.readthedocs.io/en/latest/) provides sequencing projects with the tools for gene prediction via evidence-based annotation. This is accomplished with the automatic generation of sequence features which can be refined through expert user curation. Exporting and utilizing these  annotations is the preliminary step to assigning user-curated biological function to predictions.

To extend the value of these efforts to the broader scientific community, annotations should be made available to the public in widely accessible biological databases. Upon successful submission, [GenBank](https://www.ncbi.nlm.nih.gov/genbank/) will assign accession numbers to submitted data to act as unique identifiers.  To achieve this, careful adherence to [NCBI guidelines for genome submission](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/) is essential. This pipeline is equipped with a series of checks to assess the quality of annotations and minimize errors.

The A2GB pipeline will:
1)	Reformat the GFF3 annotations from [Apollo](https://genomearchitect.readthedocs.io/en/latest/) for deposition into the [NCBI](https://www.ncbi.nlm.nih.gov/) databases. 
2)	Run sequence searches against [UniProt](https://www.uniprot.org/)'s SwissProt/TrEMBL and [InterPro](https://www.ebi.ac.uk/interpro/) databases for protein function prediction. 
3)	Run intermittent checks to assess the quality of annotations.

## Requirements
- A UNIX-like environment (UNIX/Linux, MacOS X/11, or Miscrosoft's [WSL2](https://docs.microsoft.com/en-us/windows/wsl/compare-versions))
- [Perl 5](https://www.perl.org/)
- [BioPerl](https://bioperl.org) ## Bio::SeqIO
- [Apollo](https://genomearchitect.readthedocs.io/en/latest/) (2.5.0+)
- [RNAmmer](https://services.healthtech.dtu.dk/software.php) (1.2+)
- [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) (2.0+)
- [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) (18.0.0+)
- [aria2](https://aria2.github.io/)
- [InterProScan 5](https://github.com/ebi-pf-team/interproscan) (latest version)
- [DIAMOND](https://github.com/bbuchfink/diamond) (2.0+) or [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (2.10+)
- [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/)
- [ChimeraX](https://www.cgl.ucsf.edu/chimerax/download.html)

## A2GB workflow
### Exporting annotations from Apollo
After genomic annotations are completed in Apollo, export the curated annotations.  Begin by selecting the 'Ref Sequence' tab. 
<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Apollo.png" alt="How to export Apollo annotations" width="1000"></p>

Then, select Export -> select GFF3; select ALL; select GFF3 with FASTA; click Export. The file created will be called Annotations.gff3.gz.
<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Apollo2.png" alt="How to export Apollo annotations" width="600"></p>

We recommend exporting only protein features (CDS) from Apollo. Altough rRNAs and tRNAs inferences (*e.g.* from [RNAmmer](https://services.healthtech.dtu.dk/software.php) and [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/), respectively) can be added to [Apollo](https://genomearchitect.readthedocs.io/en/latest/), the process of exporting those back is finicky and prone to errors (some rRNAs/tRNAs appear to be missing when exporting features from Apollo 2.5.0).

First, let's create a directory to store annotations:

```Bash
export ANNOT=/media/FatCat/user/raw_data     ## Replace /media/FatCat/user/raw_data by desired annotation directory
mkdir $ANNOT;                    ## Create directory $ANNOT
mv Annotations.gff3.gz $ANNOT/   ## move Annotations.gff3.gz into $ANNOT
cd $ANNOT/;
gunzip Annotations.gff3.gz       ## Decompress the GZIP file
```

Second, let's predict ribosomal RNAs with RNAmmer using [run_RNAmmer.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_RNAmmer.pl); then convert the annotations to GFF3 format with [RNAmmer_to_GFF3.pl](https://github.com/PombertLab/A2GB/blob/master/Apollo_tools/RNAmmer_to_GFF3.pl):

```Bash
mkdir $ANNOT/RNAmmer/

run_RNAmmer.pl \
   -f *.fasta \
   -d $ANNOT/RNAmmer/

RNAmmer_to_GFF3.pl \
   -g  $ANNOT/RNAmmer/*.gff2 \
   -d  $ANNOT/RNAmmer/
```


Third , let's predict transfer RNAs with tRNAscan-SE using [run_tRNAscan.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_tRNAscan.pl); then convert the annotations to GFF3 format with [tRNAscan_to_GFF3.pl](https://github.com/PombertLab/A2GB/blob/master/Apollo_tools/tRNAscan_to_GFF3.pl):

```Bash
mkdir $ANNOT/tRNAscan/

run_tRNAscan.pl \
   -f *.fasta \
   -t E \
   -d $ANNOT/tRNAscan/

tRNAscan_to_GFF3.pl \
   -t $ANNOT/tRNAscan/*.tRNAs \
   -d $ANNOT/tRNAscan/
```

Fourth, let's concatenate the tRNA, rRNA and CDS GFF3 annotations from RNAmmer, tRNAscan-SE, and Apollo:

```Bash
cat $ANNOT/RNAmmer/*.gff3 \
   $ANNOT/tRNAscan/*.gff3 \
   $ANNOT/Annotations.gff3 \
   > $ANNOT/all_annotations.gff3
## We concatenate Apollo's GFF3 file last as it includes sequence data
```

#### Splitting Apollo GFF3 files
Because debugging issues with annotations is easier when working with single files in [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/), let's split the concatenated Apollo GFF3 file into distinct GFF3 (.gff3) and FASTA (.fsa) files, one per contig/chromosome with [splitGFF3.pl](https://github.com/PombertLab/A2GB/blob/master/Apollo_tools/splitGFF3.pl).

```Bash
splitGFF3.pl \
   -g $ANNOT/all_annotations.gff3 \
   -d $ANNOT/splitGFF3
```

### Converting GFF3 files to EMBL format
This step requires a [locus_tag prefix](https://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf). If a locus_tag prefix has not been created, visit the [BioSample](https://www.ncbi.nlm.nih.gov/biosample) and [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) databases to submit all relevant sample metadata and project details. Once the sample has been accepted, the submitter will receive a BioSample accession number and a unique locus_tag prefix to be referenced during submission of corresponding experimental data to the [NCBI](https://www.ncbi.nlm.nih.gov/), [EMBL-EBI](https://www.ebi.ac.uk/) and [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) databases. Alternatively, to proceed without the locus_tag prefix, simply use a temporary prefix to be replaced later. 

Let's convert the GFF3 files to EMBL format with [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl). This script will generate locus tags automatically based on the provided prefix from NCBI. 

```Bash
ApolloGFF3toEMBL.pl \
   -p LOCUS_TAG_PREFIX \
   -g $ANNOT/splitGFF3/*.gff3 \
   -f $ANNOT/features.list \
   -c 1
```
Options for [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl) are:

```
-p (--prefix)	locus_tag prefix
-g (--gff)	GFF3 files generated by Apollo
-o (--outdir)	Output directory [Default: ./]
-f (--features)	Generate a tab-delimited list of features [Default: features.list]
-z (--zeroes)	Number of padding zeroes for locus tags [Default: 5]
-x (--exon)	Create exon features for genes with introns
-i (--intron)	Create intron features
-r (--rgb)	Change default colors of EMBL features for Artemis
-l (--lcolors)	Display a list of possible RGB colors
-c (--gcode)	NCBI genetic code [Default: 1]
		1  - The Standard Code
		2  - The Vertebrate Mitochondrial Code
		3  - The Yeast Mitochondrial Code
		4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		11 - The Bacterial, Archaeal and Plant Plastid Code
		# For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
```

After converting files to EMBL format, the content of the directory should look like this:

```Bash
ls -la $ANNOT/splitGFF3/*

-rw-rw-r--. 1 jpombert jpombert 2701562 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_01.embl
-rw-rw-r--. 1 jpombert jpombert 1886580 Oct 10 13:40 /media/FatCat/user/raw_data/splitGFF3/chromosome_01.fsa
-rw-rw-r--. 1 jpombert jpombert 1073070 Oct 10 13:40 /media/FatCat/user/raw_data/splitGFF3/chromosome_01.gff3
-rw-rw-r--. 1 jpombert jpombert  506424 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_01.prot
-rw-rw-r--. 1 jpombert jpombert 1536924 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_01.RNA
-rw-rw-r--. 1 jpombert jpombert 2573724 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_02.embl
-rw-rw-r--. 1 jpombert jpombert 1798505 Oct 10 13:40 /media/FatCat/user/raw_data/splitGFF3/chromosome_02.fsa
-rw-rw-r--. 1 jpombert jpombert 1021787 Oct 10 13:40 /media/FatCat/user/raw_data/splitGFF3/chromosome_02.gff3
-rw-rw-r--. 1 jpombert jpombert  471868 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_02.prot
-rw-rw-r--. 1 jpombert jpombert 1427737 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_02.RNA
-rw-rw-r--. 1 jpombert jpombert 2180866 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_03.embl
-rw-rw-r--. 1 jpombert jpombert 1528854 Oct 10 13:40 /media/FatCat/user/raw_data/splitGFF3/chromosome_03.fsa
-rw-rw-r--. 1 jpombert jpombert  808054 Oct 10 13:40 /media/FatCat/user/raw_data/splitGFF3/chromosome_03.gff3
-rw-rw-r--. 1 jpombert jpombert  408288 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_03.prot
-rw-rw-r--. 1 jpombert jpombert 1214029 Oct 10 14:02 /media/FatCat/user/raw_data/splitGFF3/chromosome_03.RNA
...
```

The EMBL files should resemble this:

```Bash
head -n 16 $ANNOT/splitGFF3/chromosome_01.embl

ID   chromosome_01; genomic DNA; 1855637 bp
XX
DE   Generated on Sat Mar 27 13:17:03 2021 by ApolloGFF3toEMBL.pl version 4.0
XX
FT   gene             2..148
FT                   /locus_tag="H0P50_01g00010"
FT                   /colour=255 255 255
FT   CDS             2..148
FT                   /locus_tag="H0P50_01g00010"
FT                   /colour=0 255 255
FT   gene             complement(239..3043)
FT                   /locus_tag="H0P50_01g00020"
FT                   /colour=255 255 255
FT   CDS             complement(239..3043)
FT                   /locus_tag="H0P50_01g00020"
FT                   /colour=0 255 255
```

If created properly, the files should be easy to open with [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/):

```Bash
art $ANNOT/splitGFF3/chromosome_01.embl
```

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Artemis.png" alt="Artemis opening en EMBl file generated with ApolloGFF3toEMBL.pl" width="1000"></p>

#### Checking for internal stop codons and missing start methionines
Note that [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl) will also create FASTA files of proteins and RNAs with the .prot and .RNA extensions, respectively, and which can be used for debugging issues with the corresponding annotations.

For example, we can use [check_problems.pl](https://github.com/PombertLab/A2GB/blob/master/check_problems.pl) to check for missing start methionines and for internal stop codons in proteins:

```Bash
check_problems.pl \
   -p $ANNOT/splitGFF3/*.prot \
   -v
```

If present, we should see error messages like this:
```
Checking for errors in chromosome_01.prot located in /media/FatCat/user/raw_data/splitGFF3/

		Invalid Start Codon	Internal Stop Codon

HOP50_01g00010		K			.
HOP50_01g00140		V			.
HOP50_01g07580		V			X

Checking for errors in chromosome_02.prot located in /media/FatCat/user/raw_data/splitGFF3/

		Invalid Start Codon	Internal Stop Codon

HOP50_02g14980		.			X
HOP50_02g18670		V			.

Checking for errors in chromosome_03.prot located in /media/FatCat/user/raw_data/splitGFF3/

		Invalid Start Codon	Internal Stop Codon

HOP50_03g26990		.			X
```

If present, internal stop codons and missing start methionines (wrong amino acids) can be corrected in [Apollo](https://genomearchitect.readthedocs.io/en/latest/), the GFF exported again, and the subsequent steps performed anew. Alternatively, the errors can be corrected directly on the EMBL files with [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/):

```Bash
art $ANNOT/splitGFF3/chromosome_01.embl
```

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Artemis_2.png" alt="Fixing an internal stop codon with Artemis" width="1000"></p>

We can check if the issues have been fixed by regenerating the .prot files with [check_problems.pl](https://github.com/PombertLab/A2GB/blob/master/check_problems.pl) again using the -u flag:
```Bash
check_problems.pl \
   -p $ANNOT/splitGFF3/*.prot \
   -u \
   -v
```
If fixed, the error message(s) should be gone:
```
Checking for errors in chromosome_01.prot located in /media/FatCat/user/raw_data/splitGFF3/
OK: No error found in chromosome_01.prot
```

Options for [check_problems.pl](https://github.com/PombertLab/A2GB/blob/master/check_problems.pl) are:
```
-p (--prot)	FASTA files (.prot) to be checked for abnormalities
-o (--out)	Print the output to a log file
-u (--update)	Update .prot files using EMBLtoFeatures.pl
-v (--verb)	Add verbosity
```

Additionally, features can be extracted from the EMBL files with [EMBLtoFeatures.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoFeatures.pl). Options for EMBLtoFeatures.pl are:
```
-e (--embl)	EMBL files
-o (--outdir)	Output directory [Default: ./]
-x (--exon)	Create fasta files of exons (.exons) for genes with introns
-i (--intron)	Create fasta files of exons (.introns) for genes with introns
-m (--mRNA)	Export messenger RNAs of CDS features to .mRNAs
-v (--verbose)	Add Verbosity
-c (--gcode)	NCBI genetic code [Default: 1]
		1  - The Standard Code
		2  - The Vertebrate Mitochondrial Code
		3  - The Yeast Mitochondrial Code
		4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		11 - The Bacterial, Archaeal and Plant Plastid Code
		NOTE - For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
```

#### Creating tab-delimited lists of RNA locus tags and their products
When running [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl), locus tags are generated automatically from the provided prefix. We can generate tab-delimited lists of tRNAs/rRNAs and their products from the features.list generated by ApolloGFF3toEMBL.pl and from the files located in $ANNOT/tRNAscan/ and $ANNOT/RNAmmer/ (see [section](https://github.com/PombertLab/A2GB/blob/master/README.md#Exporting-annotations-from-Apollo) above). Protein function prediction will be performed separately in the next section.

To generate tab-delimited lists of tRNAs/rRNAs, we can use [get_RNA_locus_tags.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_RNA_locus_tags.pl). This script will convert automatically the odd 8S_rRNA naming scheme produced by RNAmmer to '5S ribosomal RNA', which is recognized by NCBI. To use [get_RNA_locus_tags.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_RNA_locus_tags.pl), type:

```Bash
get_RNA_locus_tags.pl \
   -f $ANNOT/features.list \
   -ti $ANNOT/tRNAscan/*.tRNAs \
   -ri $ANNOT/RNAmmer/*.gff3 \
   -to $ANNOT/tRNA.annotations \
   -ro $ANNOT/rRNA.annotations
```

The lists created should look like this:

```
head -n 3 $ANNOT/tRNA.annotations $ANNOT/rRNA.annotations

==> /media/FatCat/user/raw_data/tRNA.annotations <==
HOP50_01g07820  tRNA-Leu(CAG)
HOP50_01g08050  tRNA-Ile(UAU)
HOP50_01g07840  tRNA-Leu(AAG)

==> /media/FatCat/user/raw_data/rRNA.annotations <==
HOP50_06g41910  28s ribosomal RNA
HOP50_06g42050  28s ribosomal RNA
HOP50_17g79280  28s ribosomal RNA

```

### Protein function prediction
In this step, individual protein sequences will be characterized using [InterProScan 5](https://github.com/ebi-pf-team/interproscan) searches, [BLASTP](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)/[DIAMOND](https://github.com/bbuchfink/diamond) searches against [UniProt](https://www.uniprot.org/)'s SwissProt/TrEMBL databases, and [BLASTP](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)/[DIAMOND](https://github.com/bbuchfink/diamond) searches against reference genome(s), if available. These annotators will help assign putative functions to predicted proteins.

First, let's generate a single multifasta file containing all of the predicted protein sequences. Ideally, internal stop codons and missing methionines should have been corrected prior to this point:

```Bash
cat $ANNOT/splitGFF3/*.prot > proteins.fasta
```

#### Predicting functions with InterProScan5
[InterPro](https://www.ebi.ac.uk/interpro/) is a free, widely used database which functionally characterizes unknown protein sequences by classifying them into families and predicts the presence of domains, repeats, and various functional sites. Unknown sequences are queried against predictive models built from identified domains and families. These models, or diagnostic signatures, are provided by InterPro’s diverse set of member databases. The result of pooling distinct signatures from member databases into a single searchable database makes InterPro a robust tool for protein functional prediction.

[InterProScan 5](https://github.com/ebi-pf-team/interproscan) can be run using the interproscan.sh script provided with its distribution or with the [run_InterProScan.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_InterProScan.pl) Perl wrapper. To run InterProScan 5 using [run_InterProScan.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_InterProScan.pl):

```Bash
run_InterProScan.pl \
   -c 10 \
   -ip \
   -go \
   -pa \
   -f $ANNOT/proteins.fasta \
   -d $ANNOT/Interproscan/ \
   -l interproscan.log
```

Options for [run_InterProScan.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_InterProScan.pl) are:

```
-c (--cpu)		Number of CPU cores to use [Default: 10]
-f (--fasta)		FASTA file(s) to query
-ip (--iprlookup)	Use InterPro's pre-calculated match lookup service
-go (--goterms)		Gene ontology search (requires --iprlookup)
-pa (--pathways)	KEGG pathways (requires --iprlookup)
-d (--dir)		Output directory (Optional)
-l (--log)		Log name [Default: interproscan.log]
```

Important, if any stop codon is present in the queries, InterProScan will throw an error message looking like this:
```
ERROR: uk.ac.ebi.interpro.scan.jms.worker.LocalJobQueueListener - 2. The exception is :java.lang.IllegalArgumentException: You have submitted a protein sequence which contains an asterix (*). This may be from an ORF prediction program. '*' is not a valid IUPAC amino acid character and amino acid sequences which go through our pipeline should not contain it. Please strip out all asterix characters from your sequence and resubmit your search.
```

If this happens, you can fix the issues with the stop codons by following the steps described in this [section](https://github.com/PombertLab/A2GB/blob/master/README.md#checking-for-internal-stop-codons-and-missing-start-methionines).

#### Performing homology searches against UniProt databases
The [UniProt](https://www.uniprot.org/) Knowledgebase (UniProtKB) is a wide-ranging database of extensively curated information of protein sequence and functional information. UniProtKB is comprised of UniProtKB/Swiss-Prot and UniProtKB/TrEMBL. Each of these offer a varying level of reliability and quality. The Swiss-Prot database contains proteins that have been tested experimentally and are manually annotated and reviewed. The TrEMBL database utilizes semi-automatic annotation, which is computationally analyzed and typically not reviewed.  Together, these databases provide a substantial collection of functional information on proteins.

Homology searches against the SwissProt and TrEMBL databases can be performed with [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or [DIAMOND](https://github.com/bbuchfink/diamond). We recommend using [DIAMOND](https://github.com/bbuchfink/diamond) due to its significantly decreased computation time.

#####  Downloading the SwissProt and TrEMBL databases
We can use [get_UniProt.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_UniProt.pl) to [download](https://www.uniprot.org/downloads) the SwissProt and/or TrEMBL databases from UniProt by leveraging [aria2](https://aria2.github.io/):

```Bash
get_UniProt.pl \
   -s \
   -t \
   -f $ANNOT/UNIPROT/ \
   -n 20 \
   -l download.log
```
Options for [get_UniProt.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_UniProt.pl) are:

```
-s  (--swiss)         Download Swiss-Prot
-t  (--trembl)        Download trEMBL
-f  (--folder)        Download folder [Default: ./]
-l  (--log)           Print download information to log file
-dt (--dtool)         Specify download tool: aria2c, wget or curl ## Tries to autodetect otherwise
-x  (--connex)        Number of aria connections [Default: 10]
-d  (--decompress)    Decompress downloaded files with gunzip     ## trEMBL is huge, off by default
-v  (--version)       Show script version
```

##### Creating tab-delimited product lists from UniProt databases
Homology searches against the [UniProt](https://www.uniprot.org/) databases will return positive matches against the corresponding accession numbers. However, these matches will not include product names. To facilitate downstream analyses, we can create tab-delimited lists of accession numbers and their products with [get_uniprot_products.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_uniprot_products.pl):

```Bash
get_uniprot_products.pl \
   -f $ANNOT/UNIPROT/uniprot_*.fasta.gz \
   -o $ANNOT/UNIPROT
```

These lists should be regenerated everytime the local UniProt databases are updated. Note that creating a product list from the TrEMBL database will take time due to its size. The tab-delimited lists should look like this:

```Bash
head -n 4 $ANNOT/UNIPROT/*.list
==> /media/FatCat/user/raw_data/UNIPROT/uniprot_sprot.list <==
sp|Q6GZX4|001R_FRG3G    Putative transcription factor 001R
sp|Q6GZX3|002L_FRG3G    Uncharacterized protein 002L
sp|Q197F8|002R_IIV3     Uncharacterized protein 002R
sp|Q197F7|003L_IIV3     Uncharacterized protein 003L

==> /media/FatCat/user/raw_data/UNIPROT/uniprot_trembl.list <==
tr|Q51723|Q51723_9EURY  Beta-galactosidase
tr|A2TI13|A2TI13_9EURY  Methyl-coenzyme M reductase (Fragment)
tr|A8USH6|A8USH6_9EURY  Methyl-coenzyme M reductase alpha (Fragment)
tr|C0LL04|C0LL04_9EURY  Methyl-coenzyme M reductase I alpha subunit (Fragment)
```

##### Running DIAMOND or BLAST searches against UniProt databases
We can use [DIAMOND](https://github.com/bbuchfink/diamond) to perform homology searches against the [UniProt](https://www.uniprot.org/) databases. Documentation on how to use DIAMOND can be found [here](https://github.com/bbuchfink/diamond/wiki).

First, let's create DIAMOND-formatted databases:

```Bash
mkdir $ANNOT/DIAMOND/; mkdir $ANNOT/DIAMOND/DB/; 

diamond makedb \
   --in $ANNOT/UNIPROT/uniprot_sprot.fasta.gz \
   -d $ANNOT/DIAMOND/DB/sprot

diamond makedb \
   --in $ANNOT/UNIPROT/uniprot_trembl.fasta.gz \
   -d $ANNOT/DIAMOND/DB/trembl
```

Second, let's perform protein-protein homology searches against the UniProt databases with a tabular output format (same as [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)'s -outfmt 6 format). Note that these searches will likely take a while depending the total number of proteins queried and the size of these databases:

```Bash
diamond blastp \
   -d $ANNOT/DIAMOND/DB/sprot \
   -q $ANNOT/proteins.fasta \
   -o $ANNOT/DIAMOND/sprot.diamond.6 \
   -f 6
   
diamond blastp \
   -d $ANNOT/DIAMOND/DB/trembl \
   -q $ANNOT/proteins.fasta \
   -o $ANNOT/DIAMOND/trembl.diamond.6 \
   -f 6
```

The result of the DIAMOND homology searches should look like this:

```Bash
head -n 4 $ANNOT/DIAMOND/*.diamond.6

==> /media/FatCat/user/raw_data/DIAMOND/sprot.diamond.6 <==
HOP50_01g00020  sp|Q54YZ9|DHKJ_DICDI    29.8    514     224     6       543     924     1340    1848    2.5e-50 202.2
HOP50_01g00020  sp|Q8D5Z6|LUXQ_VIBVU    33.3    381     234     6       542     916     476     842     6.3e-49 197.6
HOP50_01g00020  sp|Q7MD16|LUXQ_VIBVY    33.3    381     234     6       542     916     476     842     8.2e-49 197.2
HOP50_01g00020  sp|Q5A599|NIK1_CANAL    29.2    520     224     9       543     925     494     1006    2.4e-48 195.7

==> /media/FatCat/user/raw_data/DIAMOND/trembl.diamond.6 <==
HOP50_01g00010  tr|A0A5B8MBL8|A0A5B8MBL8_9CHLO  100.0   48      0       0       1       48      244     291     3.2e-19 102.8
HOP50_01g00010  tr|A0A5B8MD09|A0A5B8MD09_9CHLO  93.3    45      3       0       1       45      231     275     3.3e-16 92.8
HOP50_01g00020  tr|A0A5B8MDW7|A0A5B8MDW7_9CHLO  100.0   890     0       0       45      934     1       890     0.0e+00 1426.8
HOP50_01g00020  tr|A0A5B8MTR1|A0A5B8MTR1_9CHLO  65.7    895     299     3       42      932     26      916     2.5e-277        964.5
```

#### Performing homology searches against reference datasets
If desired, reference datasets (custom or downloaded from NCBI) can also be used as databases in homology searches to help with annotations. NCBI datasets can be accessed from the [NCBI genome](https://www.ncbi.nlm.nih.gov/genome) database or directly from their new [dataset](https://www.ncbi.nlm.nih.gov/datasets/genomes/) repository.

For example, using two datasets downloaded from NCBI:

```Bash
## Downloading data from NCBI
mkdir $ANNOT/REFERENCES/;
wget -O $ANNOT/REFERENCES/ref1.faa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/007/859/695/GCA_007859695.1_ASM785969v1/GCA_007859695.1_ASM785969v1_protein.faa.gz
wget -O $ANNOT/REFERENCES/ref2.faa.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/214/015/GCF_000214015.3_version_140606/GCF_000214015.3_version_140606_protein.faa.gz

## Creating DIAMOND database; using zcat to concatenate the output to STDOUT, then feed it to diamond as input
zcat $ANNOT/REFERENCES/*.gz | diamond makedb -d $ANNOT/DIAMOND/DB/reference

## Running DIAMOND 
diamond blastp \
   -d $ANNOT/DIAMOND/DB/reference \
   -q $ANNOT/proteins.fasta \
   -o $ANNOT/REFERENCES/reference.diamond.6 \
   -f 6
```

To create a tab-delimited list of accession numbers and their associated proteins from the downloaded NCBI .faa.gz files, we can use [get_reference_products.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_reference_products.pl).

```Bash
get_reference_products.pl \
   -f $ANNOT/REFERENCES/*.gz \
   -l $ANNOT/REFERENCES/reference.list
```

The list created should look like this:
```
head -n 5 $ANNOT/REFERENCES/reference.list

QDZ17483.1      hypothetical protein A3770_01p00010
QDZ17484.1      cytochrome P450
QDZ17485.1      hypothetical protein A3770_01p00030
QDZ17486.1      coiled-coil domain-containing protein
QDZ17487.1      putative transmembrane protein
```

#### Searching KEGG databases for orthologs
The [Kyoto Encyclopedia of Genes and Genomes](https://www.genome.jp/kegg/) (KEGG) database is a useful resource to identify which metabolic pathways are present in an organism. The KEGG databases can be queried for orthologs; proteins with matches against KEGG proteins will be assigned KO numbers (for KEGG orthologs). These can be useful during the annotation process. KEGG orthologs can be identified with [BlastKOALA](https://www.kegg.jp/blastkoala/), [GhostKOALA](https://www.kegg.jp/ghostkoala/) and/or [KofamKOALA](https://www.genome.jp/tools/kofamkoala/) using the KEGG web portal.

#### Parsing the results from homology searches
First, let's start by creating a simple list of all proteins queries, even those that returned no homology in InterProScan5, DIAMOND and/or KEGG searches. We will use [get_queries.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_queries.pl) for this: 
```Bash
get_queries.pl $ANNOT/proteins.fasta

head -n 4 $ANNOT/proteins.queries ## Looking at the list produced by get_queries.pl; a simple list with one entry per line
HOP50_01g00010
HOP50_01g00020
HOP50_01g00030
HOP50_01g00040
```

Then, let's parse the output of the InterProScan 5 and DIAMOND searches using the list of queries produced by [get_queries.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_queries.pl) and the lists of accession numbers/products created with [get_uniprot_products.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_uniprot_products.pl). We will use [parse_annotators.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/parse_annotators.pl) to do this:

```Bash
parse_annotators.pl \
   -q $ANNOT/proteins.queries \
   -o $ANNOT/proteins.annotations \
   -ip $ANNOT/Interproscan/proteins.fasta.interpro.tsv \
   -sl $ANNOT/UNIPROT/uniprot_sprot.list \
   -tl $ANNOT/UNIPROT/uniprot_trembl.list \
   -sb $ANNOT/DIAMOND/sprot.diamond.6 \
   -tb $ANNOT/DIAMOND/trembl.diamond.6
```

The parsed output should look like this:

```Bash
head -n 4 $ANNOT/proteins.annotations

#Locus_tag      Evalue  SwissProt       Evalue  trEMBL  Evalue  PFAM    Evalue  TIGR    Score   HAMAP   Evalue  CDD
HOP50_01g00010  NA      hypothetical protein    NA      hypothetical protein    NA      hypothetical protein    NA      hypothetical protein    NA      hypothetical protein    NA      no motif found
HOP50_01g00020  2.5e-50 Hybrid signal transduction histidine kinase J   0.0e+00 Signal transduction histidine kinase    5.5E-30 Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase NA      hypothetical protein    NA      hypothetical protein    1.8907E-11      HisKA
HOP50_01g00030  NA      hypothetical protein    1.0e-07 Insulin-like growth factor binding, N-terminal  1.0E-6  Putative ephrin-receptor like   NA      hypothetical protein    NA      hypothetical protein    6.31891E-8      TNFRSF
```

To use [parse_annotators.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/parse_annotators.pl) with a reference dataset, type:

```Bash
parse_annotators.pl \
   -q $ANNOT/proteins.queries \
   -o $ANNOT/proteins.annotations \
   -ip $ANNOT/Interproscan/proteins.fasta.interpro.tsv \
   -sl $ANNOT/UNIPROT/uniprot_sprot.list \
   -tl $ANNOT/UNIPROT/uniprot_trembl.list \
   -sb $ANNOT/DIAMOND/sprot.diamond.6 \
   -tb $ANNOT/DIAMOND/trembl.diamond.6 \
   -rl $ANNOT/REFERENCES/reference.list \
   -rb $ANNOT/REFERENCES/reference.diamond.6
```

We can also use multiple reference datasets if desired (file prefixes should match between the corresponding .list and .diamond.6 files). To use [parse_annotators.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/parse_annotators.pl) with multiple reference datasets, type:

```Bash
parse_annotators.pl \
   -q $ANNOT/proteins.queries \
   -o $ANNOT/proteins.annotations \
   -ip $ANNOT/Interproscan/proteins.fasta.interpro.tsv \
   -sl $ANNOT/UNIPROT/uniprot_sprot.list \
   -tl $ANNOT/UNIPROT/uniprot_trembl.list \
   -sb $ANNOT/DIAMOND/sprot.diamond.6 \
   -tb $ANNOT/DIAMOND/trembl.diamond.6 \
   -rl $ANNOT/REFERENCES/*.list \
   -rb $ANNOT/REFERENCES/*.diamond.6
```

If desired, results from [KEGG](https://www.genome.jp/kegg/) and [dbCAN2](http://bcb.unl.edu/dbCAN2/) searches can also be parsed accordingly by invoking the corresponding command line options. Current options for [parse_annotators.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/parse_annotators.pl) are:
```
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

## dbCAN2 CAZy searches: http://bcb.unl.edu/dbCAN2/
-ca	dbCAN2 output file
-cl	CAZy families list ## http://bcb.unl.edu/dbCAN2/download/Databases/CAZyDB.07302020.fam-activities.txt
```

#### Curating the protein annotations
The script [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl) was designed to faciliate comparisons between function prediction tools. It requires as input the file generated with [parse_annotators.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/parse_annotators.pl). At minimum, this file should include the results from DIAMOND (or NCBI BLAST+) BLASTP homology searches against the SwissProt/TrEMBL databases and the InterProScan5 results. [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl) will generate a user-curated tab-delimited list of locus_tags and their predicted functions. This list will be stored in a file with the .curated file extension. 

To start curating annotations with [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl), simply type:

```
curate_annotations.pl -sq $ANNOT/proteins.annotations

## Putative annotation(s) found for protein HOP50_01g00020:
	1.	SwissProt	2.5e-50		Hybrid signal transduction histidine kinase J
	2.	trEMBL		0.0e+00		Signal transduction histidine kinase
	3.	PFAM		5.5E-30		Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
	4.	TIGR		NA		hypothetical protein
	5.	HAMAP		NA		hypothetical protein
	6.	CDD		1.8907E-11	HisKA

Please enter:

	[1-6] to assign annotation
	[0] to annotate the locus as a 'hypothetical protein'
	[m] to manually annotate the locus, e.g. DUFxxx domain-containing protein
	[n] to manually annotate the locus with annotation notes, e.g. structural homolog
	[?] to mark this annotation for review and add annotation notes

	[x] to exit
```

To speed up the manual annotation process, proteins without any homology/significant hit in any of the predictors used will be annotated automatically as 'hypothetical protein'. Proteins with one or more matches identified by the predictors will show a menu like the one above. Users can enter the desired selection from the menu to annotate the  proteins accordingly. The option [m] for manual annotation will likely be useful to fix typos and/or lower/uppercase character issues in the corresponding matches. The option [n] allows the user to enter a note in addition to the manual annotation. The option [?] allows the user to flag the entry for later review. The option [x] will enable the user to quit and resume at a later stage. If the option entered is not recognized, the script will exit automatically to prevent potential problems.

To resume annotations from the last annotated proteins, simply add -r (resume) to the command line:

```
curate_annotations.pl -r -sq $ANNOT/proteins.annotations

[....................................................................................................]	7/8631

## Putative annotation(s) found for protein HOP50_01g00070:
	1.	SwissProt	4.3e-186	Eukaryotic translation initiation factor 3 subunit A
	2.	trEMBL		0.0e+00		Eukaryotic translation initiation factor 3 subunit A
	3.	PFAM		4.4E-9		PCI domain
	4.	TIGR		NA		hypothetical protein
	5.	HAMAP		24.645		Eukaryotic translation initiation factor 3 subunit A [EIF3A].
	6.	CDD		NA		no motif found

Please enter:

	[1-6] to assign annotation
	[0] to annotate the locus as a 'hypothetical protein'
	[m] to manually annotate the locus, e.g. DUFxxx domain-containing protein
	[n] to manually annotate the locus with annotation notes, e.g. structural homolog
	[?] to mark this annotation for review and add annotation notes

	[x] to exit.

Selection: 

```

The  tab-delimited list of locus_tags and their predicted functions generated by [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl) should look like this:
```Bash
head -n 7 $ANNOT/proteins.annotations.curated

HOP50_01g00010  hypothetical protein
HOP50_01g00010  hypothetical protein
HOP50_01g00010  hypothetical protein
HOP50_01g00020  signal transduction histidine kinase
HOP50_01g00030  hypothetical protein
HOP50_01g00040  hypothetical protein
HOP50_01g00050  glutamine-dependent NAD(+) synthetase
```

If a reference dataset was used, the menu from [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl) should look like this:
```
curate_annotations.pl -sq $ANNOT/proteins.annotations

[....................................................................................................]	0002/8631

## Putative annotation(s) found for protein HOP50_01g00020:
	1.	SwissProt	2.5e-50		Hybrid signal transduction histidine kinase J
	2.	trEMBL		0.0e+00		Signal transduction histidine kinase
	3.	PFAM		5.5E-30		Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
	4.	TIGR		NA		hypothetical protein
	5.	HAMAP		NA		hypothetical protein
	6.	CDD		1.8907E-11	HisKA
	7.	Ref_organism	0.0e+00		signal transduction histidine kinase

Please enter:

	[1-7] to assign annotation
	[0] to annotate the locus as a 'hypothetical protein'
	[m] to manually annotate the locus, e.g. DUFxxx domain-containing protein
	[n] to manually annotate the locus with annotation notes, e.g. structural homolog
	[?] to mark this annotation for review and add annotation notes (optional)

	[x] to exit.

Selection: m
Enter desired annotation: signal transduction histidine kinase
```

Options for [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl) are:
```
-sq (--seq_hom)		Sequence homology based annotations (generated from parse_annotators.pl)
-rd (--rcsb_3d)		3D structural homology based annotations (Generated with descriptive_GESAMT_matches.pl)
-pd (--pfam_3d)		3D structural homology annotations based on predicted stuctures (Generated with descriptive_GESAMT_matches.pl)
-cx (--chimerax)	Path to ChimeraX pdb sessions
-r (--resume)		Resume annotation from last curated locus_tag
-c (--check)		Check loci marked with '?'
-v (--verify)		Check loci marked for 3D verification
```

### Converting EMBL files to ASN format
The conversion of EMBL files to TBL format in [A2GB](https://github.com/PombertLab/A2GB) is a two step process. EMBL files are first converted to TBL format with [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl), then NCBI's [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) converts the later format to ASN.

#### Converting EMBL files to TBL format
[EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl) converts EMBL files to TBL format. This script requires a single tab-limited list of the locus tags and their predicted annotations. We can create this list by concatenating the tRNAs.annotations and rRNAs.annotations files generated [previously](https://github.com/PombertLab/A2GB#Creating-tab-delimited-lists-of-RNA-locus-tags-and-their-products) together with the curated list of proteins annotations (see [above](https://github.com/PombertLab/A2GB#curating-the-protein-annotations)). Alternatively, any tab-delimited list of locus_tags and their products can be used. 

Concatenating the annotations can be quickly performed with:
```
cat \
   $ANNOT/tRNA.annotations \
   $ANNOT/rRNA.annotations \
   $ANNOT/proteins.annotations.curated \
   > $ANNOT/verified_annotations.tsv
```

The conversion from EMBL to TBL can then be performed with [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl):
```
EMBLtoTBL.pl \
   -id ITTBIO \
   -p $ANNOT/verified_annotations.tsv \
   -embl $ANNOT/splitGFF3/*.embl \
   1> $ANNOT/STD.log \
   2> $ANNOT/ERROR.log
```

Options for [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl) are:
```
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
```

Locus tag entries missing from the tab-delimited list of annotations will be reported in the $ANNOT/ERROR.log file. Missing entries will be annotated automatically as 'hypothetical protein', 'hypothetical tRNA' or as 'hypothetical RNA' for CDS, tRNA and rRNA features, respectively. If any entry is missing, the $ANNOT/ERROR.log file will look like this:
```
head -n 6 $ANNOT/ERROR.log

Cannot find database entry for locus_tag: HOP50_01g00020
Cannot find database entry for locus_tag: HOP50_01g00020
Cannot find database entry for locus_tag: HOP50_01g00030
Cannot find database entry for locus_tag: HOP50_01g00030
Cannot find database entry for locus_tag: HOP50_01g00040
Cannot find database entry for locus_tag: HOP50_01g00040
```

The TBL files created should look like below. In this example, the first gene is incomplete in 5'; this is indicated in the table by the presence of smaller than signs (<). The codon start entry is also added to ensure proper translation of the open reading frame. 
```
head -n 25 `ls $ANNOT/splitGFF3/*.tbl | head -n 1`

>Feature chromosome_01
<1      148     gene
                        locus_tag       HOP50_01g00010
<1      148     mRNA
                        locus_tag       HOP50_01g00010
                        product hypothetical protein
                        protein_id      gnl|IITBIO|HOP50_01g00010
                        transcript_id   gnl|IITBIO|HOP50_01g00010_mRNA
<1      148     CDS
                        locus_tag       HOP50_01g00010
                        product hypothetical protein
                        protein_id      gnl|IITBIO|HOP50_01g00010
                        transcript_id   gnl|IITBIO|HOP50_01g00010_mRNA
                        codon_start     2
3043    239     gene
                        locus_tag       HOP50_01g00020
3043    239     mRNA
                        locus_tag       HOP50_01g00020
                        product signal transduction histidine kinase
                        protein_id      gnl|IITBIO|HOP50_01g00020
                        transcript_id   gnl|IITBIO|HOP50_01g00020_mRNA
3043    239     CDS
                        locus_tag       HOP50_01g00020
                        product signal transduction histidine kinase
                        protein_id      gnl|IITBIO|HOP50_01g00020
```

#### Converting TBL files to ASN format
Metadata must be included together with genome sequences during the submission process to NCBI. Although some of this metadata can be entered from the online submission form(s), it is often easier to add it beforehand while generating the ASN files. Metadata for genome submission includes taxonomic information about the source of the data being submitted, details about the sequencing experiments/computational analyses performed, and general information about the author(s) and institution(s) submitting the genomes.

##### Adding metadata to FASTA files
Taxonomic metadata can be added directly to the FASTA files. The list of modifiers that can be added directly to the FASTA definition lines can be found [here](https://www.ncbi.nlm.nih.gov/genbank/mods_fastadefline/). Mandatory modifiers include the organism name [organism=XXX].

To add metadata with [add_metadata_to_fasta.pl](https://github.com/PombertLab/A2GB/blob/master/add_metadata_to_fasta.pl) using single metadata keys, type:

```Bash
add_metadata_to_fasta.pl \
   -f $ANNOT/splitGFF3/*.fsa \
   -o 'Chloropicon primus RCC138' \
   -s RCC138 \
   -g 1
```

Alternatively, to add metadata with [add_metadata_to_fasta.pl](https://github.com/PombertLab/A2GB/blob/master/add_metadata_to_fasta.pl) and tab-delimited metadata files, type:

```Bash
add_metadata_to_fasta.pl \
   -f $ANNOT/splitGFF3/*.fsa \
   -k metakeys_NCBI.tsv \
   -c chromosomes.tsv
```

Options for [add_metadata_to_fasta.pl](https://github.com/PombertLab/A2GB/blob/master/add_metadata_to_fasta.pl) are:
```
-f (--fasta)		Specifies which FASTA files to add metadata to

## Single metadata keys
-o (--organism)		Full organism name; e.g. 'Chloropicon primus RCC138'
-s (--strain)		Strain definition; e.g. RCC138
-i (--isolate)		Isolate name; e.g. 'Pacific Isolate'
-g (--gcode)		NCBI genetic code ## https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
-m (--moltype)		NCBI moltype descriptor (e.g. genomic)

## Metadata files
-k (--keys)		Tab-delimited NCBI metadata key -> value file
-c (--chromosome)	Tab-delimited contig -> chromosome assignment file
```

The tab-delimited NCBI metadata key -> value file (e.g. [metakeys_NCBI.tsv](https://github.com/PombertLab/A2GB/blob/master/Example_files/metakeys_NCBI.tsv)) should look like this:
```
### metadata key	metadata value
organism	Chloropicon primus
strain	RCC138
moltype	genomic
gcode	1
```

The tab-delimited contig -> chromosome assignment file (e.g. [chromosomes.tsv](https://github.com/PombertLab/A2GB/blob/master/Example_files/chromosomes.tsv)) should look like this:
```
#contig_name	chromosome
chromosome_01	I
chromosome_02	II
chromosome_03	III
chromosome_04	IV
```

Once modified, the FASTA definition lines should look like this:
```
head -n 1 $ANNOT/splitGFF3/*.fsa
==> /media/FatCat/user/raw_data/splitGFF3/chromosome_01.fsa <==
>chromosome_01 [organism=Chloropicon primus][moltype=genomic][gcode=1][strain=RCC138][location=chromosome][chromosome=I]

==> /media/FatCat/user/raw_data/splitGFF3/chromosome_02.fsa <==
>chromosome_02 [organism=Chloropicon primus][moltype=genomic][gcode=1][strain=RCC138][location=chromosome][chromosome=II]

==> /media/FatCat/user/raw_data/splitGFF3/chromosome_03.fsa <==
>chromosome_03 [organism=Chloropicon primus][moltype=genomic][gcode=1][strain=RCC138][location=chromosome][chromosome=III]
```

##### Generating a GenBank submission template
NCBI provides a simple web-based tool to generate the GenBank submission template file (template.sbt) required by table2asn. To generate a template.sbt file, visit: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/.

##### Creating structured comments
More information about NCBI's stuctured comments can be found [here](https://www.ncbi.nlm.nih.gov/genbank/structuredcomment/). Files containing stuctured comments are tab-delimited; assembly-data structured comments for genomes usually looks like this:
```
StructuredCommentPrefix	##Genome-Assembly-Data-START##
Assembly Method	SPAdes v. 3.13.0; Canu v. 1.8
Assembly Name	Version 1
Long Assembly Name	RCC138 version 1
Genome Coverage	345x
Sequencing Technology	Illumina MiSeq; Oxford Nanopore
```

##### Using table2asn
[table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) is a command-line program created by NCBI to automate the creation of sequence records for submission to [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). The table2asn tool will generate the .sqn file to be used for the submission. More information about the structure and content of the TBL files generated by [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl) and required as input for table2asn can be found in NCBI's [Eukaryotic Genome Annotation Guide](https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/).

To prepare a genome for submission with [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/), we must use:
1. a single multifasta file containing all sequences to be deposited
2. a single concatenated TBL file containing all corresponding annotations

While we can run [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) on separate files to generate a single SQN output, this SQN output file would result in X separate submissions rather than a single one, <i>e.g.</i>:

```bash
## A single SQN file with 20 sequences considered separate submissions
cat test_1.sqn | grep -c 'Seq-submit ::='
20

## A single SQN file with the 20 same sequences, but this time considered a single submission 
cat test_2.sqn | grep -c 'Seq-submit ::='
1
```

To generate a single submission from a set of sequences with [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/), we can use the following:


```Bash
## Concatenating the FASTA and TBL files:
ANNOT_DIR=$ANNOT/splitGFF3/
cat $ANNOT_DIR/*.tbl > genome.tbl
cat $ANNOT_DIR/*.fsa > genome.fasta

## Generating a single submission with table2asn
table2asn \
   -t template.sbt \
   -w genome.cmt \
   -i genome.fasta \
   -f genome.tbl \
   -o output_file.sqn \
   -euk \
   -J \
   -M n \
   -V vb \
   -Z \
   -H 12/31/2024

## Checking for the presence of a single submission in the output file
cat output_file.sqn | grep -c 'Seq-submit ::='
1     ## Expected number for a single submission
```

Options for [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) used above are:
```
-t	Template file (.sbt).
-w	File (.cmt) containing Genome Assembly structured comments.
-i	Creates a single submission from indicated .fasta file
-f	Concatenated TBL file
-o	Desired output .sqn file. Sets the basename for all output files.
-euk	Asserts eukaryotic lineage for the discrepancy report tests.
-J	Delayed Genomic Product Set
-M n	Master Genome Flags: n: Normal. Combines flags for genomes submissions (invokes FATAL calls when -Z discrep is included).
-V vb	Verification: v Validates the data records; b Generates GenBank flatfiles with a .gbf suffix.
-Z	Runs the Discrepancy Report.
-H	Desired date for data release.
```

#### Checking for errors and fixing them
[table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/) generates two distinct types of error reports. The first consists of validation files with the file extension .val (generated if potential errors are found). The second report (.dr => discrepancy report) contains a detailed summary including metrics, warnings and fatal errors.

Ideally, the .val files should be absent, indicating that no error has been found. We can check if errors were potentially found with:
```Bash
ls -lh *.val

-rw-rw-r--. 1 jpombert jpombert 5.6K Dec  8 14:12 /media/FatCat/user/raw_data/splitGFF3/chromosome_01.val
-rw-rw-r--. 1 jpombert jpombert 4.7K Dec  8 14:12 /media/FatCat/user/raw_data/splitGFF3/chromosome_02.val
-rw-rw-r--. 1 jpombert jpombert 1.8K Dec  8 14:12 /media/FatCat/user/raw_data/splitGFF3/chromosome_03.val

...
```

In the above example, a few .val files exist (and their sizes are not zero), which means that errors have been detected. Errors will vary per file, obviously, but the content of a .val file should look like:

```
cat $ANNOT/splitGFF3/chromosome_01.val

WARNING: valid [SEQ_FEAT.PartialProblem] PartialLocation: 5' partial is not at start AND is not at consensus splice site FEATURE: CDS: hypothetical protein [lcl|chromosome_01:<2-148] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g00010]
WARNING: valid [SEQ_FEAT.PartialProblem] PartialLocation: 5' partial is not at start AND is not at consensus splice site FEATURE: CDS: hypothetical protein [lcl|chromosome_01:c>27956-27786] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g00140]
WARNING: valid [SEQ_FEAT.PartialProblem] PartialLocation: 3' partial is not at stop AND is not at consensus splice site FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:194078-194080, 194546->194695)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g01170]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusDonor] Splice donor consensus (GT) not found after exon ending at position 194695 of lcl|chromosome_01 FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:194078-194080, 194546->194695)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g01170]
ERROR: valid [SEQ_FEAT.SeqLocOrder] Location: Intervals out of order in SeqLoc [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g07580]
WARNING: valid [SEQ_FEAT.PartialProblem] PartialLocation: 5' partial is not at start AND is not at consensus splice site FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g07580]
WARNING: valid [SEQ_FEAT.ShortExon] Internal coding region exon is too short FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g07580]
ERROR: valid [SEQ_FEAT.InternalStop] 1 internal stops. Genetic code [1] FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g07580]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusDonor] Splice donor consensus (GT) not found after exon ending at position 1360062 of lcl|chromosome_01 FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g07580]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusAcceptor] Splice acceptor consensus (AG) not found before exon starting at position 1361096 of lcl|chromosome_01 FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g07580]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusAcceptor] Splice acceptor consensus (AG) not found before exon starting at position 1361878 of lcl|chromosome_01 FEATURE: CDS: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-1360062)] [lcl|chromosome_01: raw, dna len= 1855637] -> [gnl|ITTBIO|HOP50_01g07580]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusDonor] Splice donor consensus (GT) not found after exon ending at position 194695 of lcl|chromosome_01 FEATURE: mRNA: hypothetical protein [(lcl|chromosome_01:<194078-194080, 194546->194695)] [lcl|chromosome_01: raw, dna len= 1855637]
WARNING: valid [SEQ_FEAT.PartialProblem] PartialLocation: Start does not include first/last residue of sequence FEATURE: Gene: HOP50_01g07580 [lcl|chromosome_01:c>1361096-<1360062] [lcl|chromosome_01: raw, dna len= 1855637]
WARNING: valid [SEQ_FEAT.PartialProblem] PartialLocation: Stop does not include first/last residue of sequence FEATURE: Gene: HOP50_01g07580 [lcl|chromosome_01:c>1361096-<1360062] [lcl|chromosome_01: raw, dna len= 1855637]
ERROR: valid [SEQ_FEAT.SeqLocOrder] Location: Intervals out of order in SeqLoc [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-<1360062)] FEATURE: mRNA: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-<1360062)] [lcl|chromosome_01: raw, dna len= 1855637]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusDonor] Splice donor consensus (GT) not found after exon ending at position 1360062 of lcl|chromosome_01 FEATURE: mRNA: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-<1360062)] [lcl|chromosome_01: raw, dna len= 1855637]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusAcceptor] Splice acceptor consensus (AG) not found before exon starting at position 1361096 of lcl|chromosome_01 FEATURE: mRNA: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-<1360062)] [lcl|chromosome_01: raw, dna len= 1855637]
WARNING: valid [SEQ_FEAT.NotSpliceConsensusAcceptor] Splice acceptor consensus (AG) not found before exon starting at position 1361878 of lcl|chromosome_01 FEATURE: mRNA: hypothetical protein [(lcl|chromosome_01:c>1361096-1360062, c1361878-1361864, c1361711-1361185, c1361097-<1360062)] [lcl|chromosome_01: raw, dna len= 1855637]
ERROR: valid [SEQ_INST.StopInProtein] [1] termination symbols in protein sequence (HOP50_01g07580 - hypothetical protein) BIOSEQ: gnl|ITTBIO|HOP50_01g07580: raw, aa len= 870
```

The discrepancy report generated by [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/)'s -Z command line switch is highly verbose, and includes many checks, warnings and error messages. It is composed of two sections; a summary and a detailed section. FATAL errors will prevent submissions to NCBI and should be corrected. Howewer, FATAL: BACTERIAL_ errors, like in the example below, can be safely ignored for eukaryote genomes.

```
head -n 40 $OUTDIR/*.dr

Discrepancy Report Results

Summary
COUNT_NUCLEOTIDES: 20 nucleotide Bioseqs are present
SOURCE_QUALS: taxname (all present, all same)
SOURCE_QUALS: 20 sources have taxname = Chloropicon primus
SOURCE_QUALS: location (all present, all same)
SOURCE_QUALS: 20 sources have location = chromosome
SOURCE_QUALS: chromosome (all present, all unique)
SOURCE_QUALS: strain (all present, all same)
SOURCE_QUALS: 20 sources have strain = RCC138
FEATURE_COUNT: CDS: 8627 present
FEATURE_COUNT: gene: 8683 present
FEATURE_COUNT: mRNA: 8627 present
FEATURE_COUNT: rRNA: 10 present
FEATURE_COUNT: tRNA: 46 present
SUSPECT_PRODUCT_NAMES: 15 product_names contain suspect phrases or characters
        Putative Typo
                1 feature May contain plural
        May contain database identifier more appropriate in note; remove from product name
                14 features contains three or more numbers together that may be identifiers more appropriate in note
FATAL: BACTERIAL_PARTIAL_NONEXTENDABLE_PROBLEMS: 2 features have partial ends that do not abut the end of the sequence or a gap, and cannot be extended by 3 or fewer nucleotides to do so
FATAL: BACTERIAL_JOINED_FEATURES_NO_EXCEPTION: 3219 coding regions with joined locations have no exceptions
JOINED_FEATURES: 6470 features have joined locations.
CHROMOSOME_PRESENT: one or more chromosomes are present
FEATURE_LIST: Feature List
INCONSISTENT_BIOSOURCE: 20 inconsistent contig sources (subsource qualifiers differ)
SUSPICIOUS_SEQUENCE_ID: 20 sequences have suspicious identifiers

Detailed Report

COUNT_NUCLEOTIDES: 20 nucleotide Bioseqs are present
output_file.sqn:chromosome_01 (length 1855637)
output_file.sqn:chromosome_02 (length 1769006)
output_file.sqn:chromosome_03 (length 1503776)
output_file.sqn:chromosome_04 (length 1460551)
output_file.sqn:chromosome_05 (length 1118318)
output_file.sqn:chromosome_06 (length 1154373)
output_file.sqn:chromosome_07 (length 832941)
output_file.sqn:chromosome_08 (length 860321)
```

##### Partial genes
A common issue, especially with fragmented assemblies, is the presence of partial genes that abut the edges of contigs or chromosomes. To fix this, we must extend the feature to the edge of the contig and then, for protein-coding genes, add a tag codon_start with the proper frame (e.g. /codon_start=2). This can be done easily with [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/). [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl) will recognize these tags automatically, and adjust the TBL files accordingly.

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Partial_1.png" alt="Issue with partial gene at the start of a contig" width="1000"></p>
<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Partial_2.png" alt="Fixing the issue with Artemis" width="1000"></p>

##### Missing stop codons and GT-AG splice sites
Another issue with gene predictors is that they sometimes do not include proper stop codons for predicted protein-coding genes. This can be fixed easily with [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) by selecting the features to modify (gene + CDS), then extending them by dragging the mouse to the proper stop codon. Note that this error is often mislabelled as a GT-AG rule issue by [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/). Real issues with improper GT-AG intron/exon junctions can also be adjusted easily by drag and drop.

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Missing_sc_1.png" alt="Missing stop codon in a protein gene" width="1000"></p>
<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/MIssing_sc_2.png" alt="Fixing the issue with Artemis" width="1000"></p>

##### Fixing errors
Errors fixed with [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) can be saved easily from the graphical interface by selecting the 'File > Save All Entries' option. Then, we can simply rerun [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl) followed by [table2asn](https://www.ncbi.nlm.nih.gov/genbank/table2asn/). As the errors are getting fixed, the .val files will gradually become empty.

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Saveas_artemis.png" alt="Missing stop codon in a protein gene" width="1000"></p>

```Bash
## Running EMBLtoTBL.pl and table2asn again on the updated EMBL files
EMBLtoTBL.pl \
   -id ITTBIO \
   -p $ANNOT/verified_annotations.tsv \
   -embl $ANNOT/splitGFF3/*.embl \
   1> $ANNOT/STD.log \
   2> $ANNOT/ERROR.log

## Concatenating the FASTA and TBL files:
ANNOT_DIR=$ANNOT/splitGFF3/
cat $ANNOT_DIR/*.tbl > genome.tbl
cat $ANNOT_DIR/*.fsa > genome.fasta

## Generating a single submission with table2asn
table2asn \
   -t template.sbt \
   -w genome.cmt \
   -i genome.fasta \
   -f genome.tbl \
   -o output_file.sqn \
   -euk \
   -J \
   -M n \
   -V vb \
   -Z \
   -H 12/31/2024

## Looking at the .val files again; no such file should remain.
ls *.val
ls: cannot access '*.val': No such file or directory
```

### Submitting ASN files to GenBank
Once the errors have been corrected, the file(s) should be ready to deposit in the <b><i>Submit assembled eukaryotic and prokaryotic genomes (WGS or Complete)</b></i> segment of the [NCBI's genome submission portal](https://submit.ncbi.nlm.nih.gov/). This step will require the corresponding [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) and [BioSample](https://www.ncbi.nlm.nih.gov/biosample) information as part of the submission process.

#### <b><i>Note</b></i>
NCBI now requests the submission of haplotigs generated as part of diploid/polyploid genome assemblies in addition to the consensus genome sequence. Instructions on how to do so are available at 
[Submitting Multiple Haplotype Assemblies](https://www.ncbi.nlm.nih.gov/genbank/diploid_haps/). These pseudohaplotypes can be submitted via the "Pseudohaplotypes of a diploid/polyploid assembly" option in the genome submission portal (each require their own BioProject however). Each file should also be named distinctively (<i>e.g.</i> genome.principal_haplotype.sqn / genome.alternate_haplotype.sqn) based on their relationships for greater clarity.

## Funding and acknowledgments
This work was supported in part by the National Institute of Allergy and Infectious Diseases of the National Institutes of Health (award number R15AI128627) to Jean-Francois Pombert. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## References
Buchfink B, Xie C, Huson DH. **Fast and sensitive protein alignment using DIAMOND.** *Nat Methods.* 2015 Jan;12(1):59-60. Epub 2014 Nov 17. PMID: 25402007 DOI: [10.1038/nmeth.3176](https://doi.org/10.1038/nmeth.3176).

Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. **Basic local alignment search tool.** *J Mol Biol.* 1990 Oct 5;215(3):403-10. PMID: 2231712 DOI: [10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/s0022-2836(05)80360-2).

Barrett T, Clark K, Gevorgyan R, Gorelenkov V, Gribov E, Karsch-Mizrachi I, Kimelman M, Pruitt KD, Resenchuk S, Tatusova T, Yaschenko E, Ostell J. **BioProject and BioSample databases at NCBI: facilitating capture and organization of metadata.** *Nucleic Acids Res.* 2012 Jan;40(Database issue):D57-63. Epub 2011 Dec 1. PMID: 22139929 PMCID: [PMC3245069](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc3245069/) DOI: [10.1093/nar/gkr1163](https://doi.org/10.1093/nar/gkr1163).

Dunn NA, Unni DR, Diesh C, Munoz-Torres M, Harris NL, Yao E, Rasche H, Holmes IH, Elsik CG, Lewis SE. **Apollo: Democratizing genome annotation.** *PLoS Comput Biol.* 2019 Feb 6;15(2):e1006790. PMID: 30726205 PMCID: [PMC6380598](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc6380598/) DOI: [10.1371/journal.pcbi.1006790](https://doi.org/10.1371/journal.pcbi.1006790). 

Carver T, Harris SR, Berriman M, Parkhill J, McQuillan JA. **Artemis: an integrated platform for visualization and analysis of high-throughput sequence-based experimental data.** *Bioinformatics.* 2012 Feb 15;28(4):464-9. Epub 2011 Dec 22. PMID: 22199388; PMCID: [PMC3278759](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc3278759/). DOI: [10.1093/bioinformatics/btr703](https://doi.org/10.1093/bioinformatics/btr703).

Chan PP, Lowe TM. **tRNAscan-SE: Searching for tRNA Genes in Genomic Sequences.** *Methods Mol Biol.* 2019;1962:1-14. PMID: 31020551 PMCID: [PMC6768409](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc6768409/) DOI: [10.1007/978-1-4939-9173-0_1](https://doi.org/10.1007/978-1-4939-9173-0_1). 

Lagesen K, Hallin P, Rødland EA, Staerfeldt HH, Rognes T, Ussery DW. **RNAmmer: consistent and rapid annotation of ribosomal RNA genes.** *Nucleic Acids Res.* 2007;35(9):3100-8. Epub 2007 Apr 22. PMID: 17452365 PMCID: [PMC1888812](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc1888812/) DOI: [10.1093/nar/gkm160](https://doi.org/10.1093/nar/gkm160).

Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, McWilliam H, Maslen J, Mitchell A, Nuka G, Pesseat S, Quinn AF, Sangrador-Vegas A, Scheremetjew M, Yong SY, Lopez R, Hunter S. **InterProScan 5: genome-scale protein function classification.** *Bioinformatics.* 2014 May 1;30(9):1236-40. Epub 2014 Jan 21. PMID: 24451626 PMCID: [PMC3998142](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc3998142/) DOI: [10.1093/bioinformatics/btu031](https://doi.org/10.1093/bioinformatics/btu031). 

Mitchell AL, Attwood TK, Babbitt PC, Blum M, Bork P, Bridge A, Brown SD, Chang HY, El-Gebali S, Fraser MI, Gough J, Haft DR, Huang H, Letunic I, Lopez R, Luciani A, Madeira F, Marchler-Bauer A, Mi H, Natale DA, Necci M, Nuka G, Orengo C, Pandurangan AP, Paysan-Lafosse T, Pesseat S, Potter SC, Qureshi MA, Rawlings ND, Redaschi N, Richardson LJ, Rivoire C, Salazar GA, Sangrador-Vegas A, Sigrist CJA, Sillitoe I, Sutton GG, Thanki N, Thomas PD, Tosatto SCE, Yong SY, Finn RD. **InterPro in 2019: improving coverage, classification and access to protein sequence annotations.** *Nucleic Acids Res.* 2019 Jan 8;47(D1):D351-D360. PMID: 30398656 PMCID: [PMC6323941](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc6323941/) DOI: [10.1093/nar/gky1100](https://doi.org/10.1093/nar/gky1100).

UniProt Consortium. **UniProt: a worldwide hub of protein knowledge.** Nucleic Acids Res. 2019 Jan 8;47(D1):D506-D515. PMID: 30395287; PMCID: [PMC6323992](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc6323992/). DOI: [10.1093/nar/gky1049](https://doi.org/10.1093/nar/gky1049). 

Sayers EW, Cavanaugh M, Clark K, Ostell J, Pruitt KD, Karsch-Mizrachi I. **GenBank.** *Nucleic Acids Res.* 2020 Jan 8;48(D1):D84-D86. PMID: 31665464 PMCID: [PMC7145611](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc7145611/) DOI: [10.1093/nar/gkz956](https://doi.org/10.1093/nar/gkz956)

Stein L. **Genome annotation: from sequence to biology.** *Nat Rev Genet.* 2001 Jul;2(7):493-503. DOI: [10.1038/35080529](https://doi.org/10.1038/35080529). PMID: 11433356. DOI: [10.1038/35080529](https://doi.org/10.1038/35080529)

Kanehisa M, Sato Y, Morishima K. **BlastKOALA and GhostKOALA: KEGG tools for functional characterization of genome and metagenome sequences.** *J Mol Biol.* 2016;428:726–731. PMID: 31423653 PMCID: [PMC6933857](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc6933857/) DOI: [10.1002/pro.3711](https://doi.org/10.1002/pro.3711)

Suzuki S, Kakuta M, Ishida T, Akiyama Y. **GHOSTX: An improved sequence homology search algorithm using a query suffix array and a database suffix array.** *PLoS One.* 2014;9(8):e103833 PMID: 25099887 PMCID: [PMC4123905](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc4123905/) DOI: [10.1371/journal.pone.0103833](https://doi.org/10.1371/journal.pone.0103833)

Aramaki T, Blanc‐Mathieu R, Endo H, et al. **KofamKOALA: KEGG ortholog assignment based on profile HMM and adaptive score threshold.** *Bioinformatics.* 2020 Apr 1;36(7):2251-2252. PMID: 31742321 PMCID: [PMC7141845](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc7141845/) DOI: [10.1093/bioinformatics/btz859](https://doi.org/10.1093/bioinformatics/btz859)

Zhang H, Yohe T, Huang L, Entwistle S, Wu P, Yang Z, Busk PK, Xu Y, Yin Y. **dbCAN2: a meta server for automated carbohydrate-active enzyme annotation.** *Nucleic Acids Res.* 2018 Jul 2;46(W1):W95-W101. PMID: 29771380 PMCID: [PMC6031026](http://www.ncbi.nlm.nih.gov/pmc/articles/pmc6031026/) DOI: [10.1093/nar/gky418](https://doi.org/10.1093/nar/gky418)
