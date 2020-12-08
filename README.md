# A2GB – Annotations to GenBank

[A2GB](https://github.com/PombertLab/A2GB/) is a pipeline that will transform the genome annotation files exported from [Apollo](https://genomearchitect.readthedocs.io/en/latest/) (fomerly known as [WebApollo](http://gmod.org/wiki/WebApollo)) in preparation for sequence submission to NCBI’s [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). While its primary goal is to submit annotations to GenBank for the generation of accession numbers, the conversion of file formats will permit the use of different tools in downstream analyses. 

In a stepwise approach, the A2GB pipeline converts annotation files from GFF3 -> EMBL -> TBL -> ASN. Each format will be useful for diagnostic quality checks of the annotations or become the springboard for other analyses, such as protein function prediction.  

Furthermore, A2GB acts as a guide to prepare sequence submissions according to [NCBI’s guidelines](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit/), including project registration with [BioSample](https://www.ncbi.nlm.nih.gov/biosample) and [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) for the generation of locus_tag prefixes. Upon acceptance from NCBI, these essential steps will ensure that sequence data will be made publicly available through GenBank and other member databases within the [International Nucleotide Sequence Database Collaboration](http://www.insdc.org/).

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
        *	[Predicting functions with InterProScan 5](#Predicting-functions-with-InterProScan-5)
        *	[Performing homology searches against UniProt databases](#Performing-homology-searches-against-UniProt-databases)
	        *	[Downloading the SwissProt and TrEMBL databases](#Downloading-the-SwissProt-and-TrEMBL-databases)
       		*	[Creating tab-delimited product lists from UniProt databases](#Creating-tab-delimited-product-lists-from-UniProt-databases)
      		*	[Running DIAMOND or BLAST searches against UniProt databases](#Running-DIAMOND-or-BLAST-searches-against-UniProt-databases)
        *	[Performing homology searches against reference datasets](#Performing-homology-searches-against-reference-datasets)
        *	[Parsing the result of InterProScan 5 and DIAMOND searches](#Parsing-the-result-of-InterProScan-5-and-DIAMOND-searches)
        *	[Curating the protein annotations](#Curating-the-protein-annotations)
   *	[Converting EMBL files to ASN format](#Converting-EMBL-files-to-ASN-format)
        *	[Converting EMBL files to TBL format](#Converting-EMBL-files-to-TBL-format)
        *	[Converting TBL files to ASN format](#Converting-TBL-files-to-ASN-format)
	        *	[Adding metadata to FASTA files](#Adding-metadata-to-FASTA-files)
      		*	[Generating a GenBank submission template](#Generating-a-GenBank-submission-template)
        	*	[Creating structured comments](#Creating-structured-comments)
        	*	[Using TBL2ASN](#Using-TBL2ASN)
        *	[Checking for errors](#Checking-for-errors)
   *	[Submitting ASN file to GenBank]
   *	[Miscellaneous] 
*	[References]

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
- [Apollo](https://genomearchitect.readthedocs.io/en/latest/) (2.5.0+)
- [RNAmmer](https://services.healthtech.dtu.dk/software.php) (1.2+)
- [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) (2.0+)
- [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) (18.0.0+)
- [InterProScan 5](https://github.com/ebi-pf-team/interproscan) (latest version)
- [DIAMOND](https://github.com/bbuchfink/diamond) (2.0+) or [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (2.10+)
- [TBL2ASN](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/)

### A2GB workflow
#### Exporting annotations from Apollo
After genomic annotations are completed in Apollo, export the curated annotations.  Begin by selecting the 'Ref Sequence' tab. 
<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Apollo.png" alt="How to export Apollo annotations" width="1000"></p>

Then, select Export -> select GFF3; select ALL; select GFF3 with FASTA; click Export. The file created will be called Annotations.gff3.gz.
<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Apollo2.png" alt="How to export Apollo annotations" width="600"></p>

We recommend exporting only protein features (CDS) from Apollo. Altough rRNAs and tRNAs inferences (*e.g.* from [RNAmmer](https://services.healthtech.dtu.dk/software.php) and [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/), respectively) can be added to [Apollo](https://genomearchitect.readthedocs.io/en/latest/), the process of exporting those back is finicky and prone to errors (some rRNAs/tRNAs appear to be missing when exporting features from Apollo 2.5.0).

First, let's create a directory to store annotations:

```Bash
export ANNOT=/media/FatCat/user/raw_data     ## Replace /media/FatCat/user/raw_data by desired annotation directory
mkdir $ANNOT; mv Annotations.gff3.gz $ANNOT/ ## Create directory; then move Annotations.gff3.gz into it
cd $ANNOT/; gunzip Annotations.gff3.gz       ## Decompress the GZIP file
```

Second, let's predict ribosomal RNAs with RNAmmer using [run_RNAmmer.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_RNAmmer.pl); then convert the annotations to GFF3 format with [RNAmmer_to_GFF3.pl](https://github.com/PombertLab/A2GB/blob/master/Apollo_tools/RNAmmer_to_GFF3.pl):

```Bash
mkdir $ANNOT/RNAmmer/
run_RNAmmer.pl -f *.fasta -d $ANNOT/RNAmmer/
RNAmmer_to_GFF3.pl -g  $ANNOT/RNAmmer/*.gff2 -d  $ANNOT/RNAmmer/
```


Third , let's predict transfer RNAs with tRNAscan-SE using [run_tRNAscan.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_tRNAscan.pl); then convert the annotations to GFF3 format with [tRNAscan_to_GFF3.pl](https://github.com/PombertLab/A2GB/blob/master/Apollo_tools/tRNAscan_to_GFF3.pl):

```Bash
mkdir $ANNOT/tRNAscan/
run_tRNAscan.pl -f *.fasta -t E -d $ANNOT/tRNAscan/
tRNAscan_to_GFF3.pl -t $ANNOT/tRNAscan/*.tRNAs -d $ANNOT/tRNAscan/
```

Fourth, let's concatenate the tRNA, rRNA and CDS GFF3 annotations from RNAmmer, tRNAscan-SE, and Apollo:

```Bash
cat $ANNOT/RNAmmer/*.gff3 $ANNOT/tRNAscan/*.gff3 $ANNOT/Annotations.gff3 > $ANNOT/all_annotations.gff3
## We concatenate Apollo's GFF3 file last as it includes sequence data
```

##### Splitting Apollo GFF3 files
Because debugging issues with annotations is easier when working with single files, let's split the concatenated Apollo GFF3 file into distinct GFF3 (.gff3) and FASTA (.fsa) files, one per contig/chromosome with [splitGFF3.pl](https://github.com/PombertLab/A2GB/blob/master/Apollo_tools/splitGFF3.pl).

```Bash
splitGFF3.pl -g $ANNOT/all_annotations.gff3 -d $ANNOT/splitGFF3
```

#### Converting GFF3 files to EMBL format
This step requires a [locus_tag prefix](https://www.ncbi.nlm.nih.gov/genomes/locustag/Proposal.pdf). If a locus_tag prefix has not been created, visit the [BioSample](https://www.ncbi.nlm.nih.gov/biosample) and [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) databases to submit all relevant sample metadata and project details. Once the sample has been accepted, the submitter will receive a BioSample accession number and a unique locus_tag prefix to be referenced during submission of corresponding experimental data to the [NCBI](https://www.ncbi.nlm.nih.gov/), [EBI](https://www.ebi.ac.uk/) and [DDBJ](https://www.ddbj.nig.ac.jp/index-e.html) databases. Alternatively, to proceed without the locus_tag prefix, simply use a temporary prefix to be replaced later. 

Let's convert the GFF3 files to EMBL format with [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl). This script will generate locus tags automatically based on the provided prefix from NCBI. 

```Bash
ApolloGFF3toEMBL.pl -p LOCUS_TAG_PREFIX -g $ANNOT/splitGFF3/*.gff3 -f  $ANNOT/features.list -c 1
```
Options for [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl) are:

```
-p (--prefix)	## locus_tag prefix
-g (--gff)	## GFF3 files generated by Apollo
-f (--features) ## Generates a tab-delimited list of features [Default: features.list]
-c (--gcode)	## NCBI genetic code [Default: 1], e.g:
		1  - The Standard Code
		2  - The Vertebrate Mitochondrial Code
		3  - The Yeast Mitochondrial Code
		4  - The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma/Spiroplasma Code
		11 - The Bacterial, Archaeal and Plant Plastid Code
		NOTE - For complete list; see https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
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

FT   gene             2..148
FT                   /locus_tag="HOP50_01g00010"
FT   CDS             2..148
FT                   /locus_tag="HOP50_01g00010"
FT   gene             complement(239..3043)
FT                   /locus_tag="HOP50_01g00020"
FT   CDS             complement(239..3043)
FT                   /locus_tag="HOP50_01g00020"
FT   gene             3823..7741
FT                   /locus_tag="HOP50_01g00030"
FT   CDS             join(3823..4103,4213..7741)
FT                   /locus_tag="HOP50_01g00030"
FT   gene             complement(7814..11050)
FT                   /locus_tag="HOP50_01g00040"
FT   CDS             complement(join(7814..8530,8615..11050)
FT                   /locus_tag="HOP50_01g00040"
```

If created properly, the files should be easy to open with [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/):

```Bash
art $ANNOT/splitGFF3/chromosome_01.embl
```

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Artemis.png" alt="Artemis opening en EMBl file generated with ApolloGFF3toEMBL.pl" width="1000"></p>

##### Checking for internal stop codons and missing start methionines
Note that [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl) will also create FASTA files of proteins and RNAs with the .prot and .RNA extensions, respectively, and which can be used for debugging issues with the corresponding annotations.

For example, we can use [check_problems.pl](https://github.com/PombertLab/A2GB/blob/master/check_problems.pl) to check for missing start methionines and for internal stop codons in proteins:

```Bash
check_problems.pl -s -m -f $ANNOT/splitGFF3/*.prot
```

If present, we should see error messages like this:
```
Checking for internal stop codons in chromosome_01.prot located in /media/FatCat/user/raw_data/splitGFF3/
ERROR: Protein HOP50_01g07580 contains one or more stop codon(s) in /media/FatCat/user/raw_data/splitGFF3/chromosome_01.prot

Checking for internal stop codons in chromosome_02.prot located in /media/FatCat/user/raw_data/splitGFF3/
ERROR: Protein HOP50_02g14900 contains one or more stop codon(s) in /media/FatCat/user/raw_data/splitGFF3/chromosome_02.prot

Checking for internal stop codons in chromosome_03.prot located in /media/FatCat/user/raw_data/splitGFF3/
OK: No internal stop codon found

...

Checking for missing start methionines in chromosome_01.prot located in /media/FatCat/user/raw_data/splitGFF3/
ERROR: Protein HOP50_01g07580 starts with V in /media/FatCat/user/raw_data/splitGFF3/chromosome_01.prot
ERROR: Protein HOP50_01g00140 starts with V in /media/FatCat/user/raw_data/splitGFF3/chromosome_01.prot
ERROR: Protein HOP50_01g00010 starts with K in /media/FatCat/user/raw_data/splitGFF3/chromosome_01.prot

Checking for missing start methionines in chromosome_02.prot located in /media/FatCat/user/raw_data/splitGFF3/
ERROR: Protein HOP50_02g18530 starts with V in /media/FatCat/user/raw_data/splitGFF3/chromosome_02.prot

Checking for missing start methionines in chromosome_03.prot located in /media/FatCat/user/raw_data/splitGFF3/
OK: All proteins start with methionines...
```

If present, internal stop codons and missing start methionines can be corrected in [Apollo](https://genomearchitect.readthedocs.io/en/latest/), the GFF exported again, and the subsequent steps performed anew. Alternatively, the errors can be corrected directly on the EMBL files with [Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/), then the .prot files regenerated with [EMBLtoPROT.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoPROT.pl):

```Bash
art $ANNOT/splitGFF3/chromosome_01.embl
```

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Artemis_2.png" alt="Fixing an internal stop codon with Artemis" width="1000"></p>

We can check if the issues have been fixed by regenerating the .prot files with [EMBLtoPROT.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoPROT.pl), then by running [check_problems.pl](https://github.com/PombertLab/A2GB/blob/master/check_problems.pl) again:
```Bash
EMBLtoPROT.pl -e $ANNOT/splitGFF3/*.embl -c 1
check_problems.pl -s -m -f $ANNOT/splitGFF3/*.prot
```
If fixed, the error message(s) should be gone:
```
Checking for internal stop codons in chromosome_01.prot located in /media/FatCat/user/raw_data/splitGFF3/
OK: No internal stop codon found
```

##### Creating tab-delimited lists of RNA locus tags and their products
When running [ApolloGFF3toEMBL.pl](https://github.com/PombertLab/A2GB/blob/master/ApolloGFF3toEMBL.pl), locus tags are generated automatically from the provided prefix. We can generate tab-delimited lists of tRNAs/rRNAs and their products from the features.list generated by ApolloGFF3toEMBL.pl and from the files located in $ANNOT/tRNAscan/ and $ANNOT/RNAmmer/ (see [section](https://github.com/PombertLab/A2GB/blob/master/README.md#Exporting-annotations-from-Apollo) above). Protein function prediction will be performed separately in the next section.

To generate tab-delimited lists of tRNAs/rRNAs, simply type:

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

#### Protein function prediction
In this step, individual protein sequences will be characterized using [InterProScan 5](https://github.com/ebi-pf-team/interproscan) searches, [BLASTP](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)/[DIAMOND](https://github.com/bbuchfink/diamond) searches against [UnitProt](https://www.uniprot.org/)'s SwissProt/TrEMBL databases, and [BLASTP](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)/[DIAMOND](https://github.com/bbuchfink/diamond) searches against reference genome(s), if available. These annotators will help assign putative functions to predicted proteins.

First, let's generate a single multifasta file containing all of the predicted protein sequences. Ideally, internal stop codons and missing methionines should have been corrected prior to this point:

```Bash
cat $ANNOT/splitGFF3/*.prot > proteins.fasta
```

##### Predicting functions with InterProScan 5
[InterPro](https://www.ebi.ac.uk/interpro/) is a free, widely used database which functionally characterizes unknown protein sequences by classifying them into families and predicts the presence of domains, repeats, and various functional sites. Unknown sequences are queried against predictive models built from identified domains and families. These models, or diagnostic signatures, are provided by InterPro’s diverse set of member databases. The result of pooling distinct signatures from member databases into a single searchable database makes InterPro a robust tool for protein functional prediction.

[InterProScan 5](https://github.com/ebi-pf-team/interproscan) can be run using the interproscan.sh script provided with its distribution or with the [run_InterProScan.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_InterProScan.pl) Perl wrapper. To run InterProScan 5 using [run_InterProScan.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/run_InterProScan.pl):

```Bash
run_InterProScan.pl -c 10 -ip -go -pa -f $ANNOT/proteins.fasta -d $ANNOT/Interproscan/ -l interproscan.log
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

##### Performing homology searches against UniProt databases
The [UniProt](https://www.uniprot.org/) Knowledgebase (UniProtKB) is a wide-ranging database of extensively curated information of protein sequence and functional information. UniProtKB is comprised of UniProtKB/Swiss-Prot and UniProtKB/TrEMBL. Each of these offer a varying level of reliability and quality. The Swiss-Prot database contains proteins that have been tested experimentally and are manually annotated and reviewed. The TrEMBL database utilizes semi-automatic annotation, which is computationally analyzed and typically not reviewed.  Together, these databases provide a substantial collection of functional information on proteins.

Homology searches against the SwissProt and TrEMBL databases can be performed with [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or [DIAMOND](https://github.com/bbuchfink/diamond). We recommend using [DIAMOND](https://github.com/bbuchfink/diamond) due to its significantly decreased computation time.

######  Downloading the SwissProt and TrEMBL databases
We can use [get_UniProt.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_UniProt.pl) to [download](https://www.uniprot.org/downloads) the SwissProt and/or TrEMBL databases from UniProt:

```Bash
get_UniProt.pl -s -t -f $ANNOT/UNIPROT/ -n 20 -l download.log
```
Options for [get_UniProt.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_UniProt.pl) are:

```
-s (--swiss)		Download Swiss-Prot
-t (--trembl)		Download trEMBL
-f (--folder)		Download folder [Default: ./]
-n (--nice)		Linux Process Priority [Default: 20] ## Runs downloads in the background
-l (--log)		Print download information to log file
-d (--decompress)	Decompresss downloaded files with gunzip ## trEMBL files will be huge, off by default
```

###### Creating tab-delimited product lists from UniProt databases
Homology searches against the [UniProt](https://www.uniprot.org/) databases will return positive matches against the corresponding accession numbers. However, these matches will not include product names. To facilitate downstream analyses, we can create tab-delimited lists of accession numbers and their products with [get_uniprot_products.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_uniprot_products.pl):

```Bash
get_uniprot_products.pl $ANNOT/UNIPROT/uniprot_*.fasta.gz
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

###### Running DIAMOND or BLAST searches against UniProt databases
We can use [DIAMOND](https://github.com/bbuchfink/diamond) to perform homology searches against the [UniProt](https://www.uniprot.org/) databases. Documentation on how to use DIAMOND can be found [here](http://www.diamondsearch.org/index.php).

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
   -o $ANNOT/DIAMOND/diamond.sprot.6 \
   -f 6
   
diamond blastp \
   -d $ANNOT/DIAMOND/DB/trembl \
   -q $ANNOT/proteins.fasta \
   -o $ANNOT/DIAMOND/diamond.trembl.6 \
   -f 6
```

The result of the DIAMOND homology searches should look like this:

```Bash
head -n 4 $ANNOT/DIAMOND/diamond.*.6

==> /media/FatCat/user/raw_data/DIAMOND/diamond.sprot.6 <==
HOP50_01g00020  sp|Q54YZ9|DHKJ_DICDI    29.8    514     224     6       543     924     1340    1848    2.5e-50 202.2
HOP50_01g00020  sp|Q8D5Z6|LUXQ_VIBVU    33.3    381     234     6       542     916     476     842     6.3e-49 197.6
HOP50_01g00020  sp|Q7MD16|LUXQ_VIBVY    33.3    381     234     6       542     916     476     842     8.2e-49 197.2
HOP50_01g00020  sp|Q5A599|NIK1_CANAL    29.2    520     224     9       543     925     494     1006    2.4e-48 195.7

==> /media/FatCat/user/raw_data/DIAMOND/diamond.trembl.6 <==
HOP50_01g00010  tr|A0A5B8MBL8|A0A5B8MBL8_9CHLO  100.0   48      0       0       1       48      244     291     3.2e-19 102.8
HOP50_01g00010  tr|A0A5B8MD09|A0A5B8MD09_9CHLO  93.3    45      3       0       1       45      231     275     3.3e-16 92.8
HOP50_01g00020  tr|A0A5B8MDW7|A0A5B8MDW7_9CHLO  100.0   890     0       0       45      934     1       890     0.0e+00 1426.8
HOP50_01g00020  tr|A0A5B8MTR1|A0A5B8MTR1_9CHLO  65.7    895     299     3       42      932     26      916     2.5e-277        964.5
```

##### Performing homology searches against reference datasets

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
   -o $ANNOT/DIAMOND/diamond.reference.6 \
   -f 6
```

To create a tab-delimited list of accession numbers and their associated proteins from the downloaded NCBI .faa.gz files, we can use [get_reference_products.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_reference_products.pl).

```Bash
get_reference_products.pl -f $ANNOT/REFERENCES/*.gz -l $ANNOT/REFERENCES/reference.list
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

##### Parsing the result of InterProScan 5 and DIAMOND searches

First, let's start by creating a simple list of all proteins queries, even those that returned no homology in InterProScan 5 and/or DIAMOND searches. We will use [get_queries.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_queries.pl) for this: 
```Bash
get_queries.pl $ANNOT/proteins.fasta

head -n 4 $ANNOT/proteins.queries ## Looking at the list produced by get_queries.pl; a simple list with one entry per line
HOP50_01g00010
HOP50_01g00020
HOP50_01g00030
HOP50_01g00040
```

Then, let's parse the output of the InterProScan 5 and DIAMOND searches using the list of queries produced by [get_queries.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_queries.pl) and the lists of accession numbers/products created with [get_uniprot_products.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/get_uniprot_products.pl). We will use [parse_annotators.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/parse_annotators.pl) to do this (Note that the uniprot_trembl.list file will be large and will eat up at least 5 Gb of RAM):

```Bash
parse_annotators.pl \
   -q $ANNOT/proteins.queries \
   -o $ANNOT/proteins.annotations \
   -sl $ANNOT/UNIPROT/uniprot_sprot.list \
   -tl $ANNOT/UNIPROT/uniprot_trembl.list \
   -sb $ANNOT/DIAMOND/diamond.sprot.6 \
   -tb $ANNOT/DIAMOND/diamond.trembl.6 \
   -ip $ANNOT/Interproscan/proteins.fasta.interpro.tsv
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
   -sl $ANNOT/UNIPROT/uniprot_sprot.list \
   -tl $ANNOT/UNIPROT/uniprot_trembl.list \
   -sb $ANNOT/DIAMOND/diamond.sprot.6 \
   -tb $ANNOT/DIAMOND/diamond.trembl.6 \
   -ip $ANNOT/Interproscan/proteins.fasta.interpro.tsv \
   -rl $ANNOT/REFERENCES/reference.list \
   -rb $ANNOT/DIAMOND/diamond.reference.6
```

##### Curating the protein annotations
The script [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl) was designed to faciliate comparisons between function prediction tools. It requires as input the file generated with [parse_annotators.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/parse_annotators.pl). At minimum, this file should include the result of DIAMOND (or NCBI BLAST+) BLASTP homology searches against the SwissProt/TrEMBL databases and the result of InterProScan 5 searches. [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl) will generate a user-curated tab-delimited list of locus_tags and their predicted functions. This list will be stored in a file with the .curated file extension. 

To start curating annotations with [curate_annotations.pl](https://github.com/PombertLab/A2GB/blob/master/Function_prediction/curate_annotations.pl), simply type:

```
curate_annotations.pl -i $ANNOT/proteins.annotations

Putative annotation(s) found for protein #0002: HOP50_01g00020:
1.      SWISSPROT:      2.5e-50         Hybrid signal transduction histidine kinase J
2.      TREMBL:         0.0e+00         Signal transduction histidine kinase
3.      Pfam:           5.5E-30         Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
4.      TIGRFAM:        NA              hypothetical protein
5.      HAMAP:          NA              hypothetical protein
6.      CDD:            1.8907E-11      HisKA

Please enter selection [1-6] to assign annotation, [0] to annotate as 'hypothetical protein', [m] for manual annotation, or [x] to exit
m
Enter desired annotation: signal transduction histidine kinase

Putative annotation(s) found for protein #0003: HOP50_01g00030:
1.      SWISSPROT:      NA              hypothetical protein
2.      TREMBL:         1.0e-07         Insulin-like growth factor binding, N-terminal
3.      Pfam:           1.0E-6          Putative ephrin-receptor like
4.      TIGRFAM:        NA              hypothetical protein
5.      HAMAP:          NA              hypothetical protein
6.      CDD:            6.31891E-8      TNFRSF

Please enter selection [1-6] to assign annotation, [0] to annotate as 'hypothetical protein', [m] for manual annotation, or [x] to exit
x
```

To speed up the manual annotation process, proteins without any homology/significant hit in any of the predictors used will be annotated automatically as 'hypothetical protein'. Proteins with one or more matches identified by the predictors will show a menu like the one above. Users can enter the desired selection from the menu to annotate the  proteins accordingly. The option [m] for manual annotation will likely be useful to fix typos and/or lower/uppercase character issues in the corresponding matches. To option [x] will enable the user to quit and resume at a later stage. If the option entered is not recognized, the script will exit automatically to prevent potential problems.

To resume annotations from the last annotated proteins, simply add -r (resume) to the command line:

```
curate_annotations.pl -r -i $ANNOT/proteins.annotations

Putative annotation(s) found for protein #0003: HOP50_01g00030:
1.      SWISSPROT:      NA              hypothetical protein
2.      TREMBL:         1.0e-07         Insulin-like growth factor binding, N-terminal
3.      Pfam:           1.0E-6          Putative ephrin-receptor like
4.      TIGRFAM:        NA              hypothetical protein
5.      HAMAP:          NA              hypothetical protein
6.      CDD:            6.31891E-8      TNFRSF

Please enter selection [1-6] to assign annotation, [0] to annotate as 'hypothetical protein', [m] for manual annotation, or [x] to exit
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
curate_annotations.pl -i $ANNOT/proteins.annotations

Putative annotation(s) found for protein #0002: HOP50_01g00020:
1.      SWISSPROT:      2.5e-50         Hybrid signal transduction histidine kinase J
2.      TREMBL:         0.0e+00         Signal transduction histidine kinase
3.      Pfam:           5.5E-30         Histidine kinase-, DNA gyrase B-, and HSP90-like ATPase
4.      TIGRFAM:        NA              hypothetical protein
5.      HAMAP:          NA              hypothetical protein
6.      CDD:            1.8907E-11      HisKA
7.      Reference:      0.0e+00         signal transduction histidine kinase

Please enter selection [1-7] to assign annotation, [0] to annotate as 'hypothetical protein', [m] for manual annotation, or [x] to exit
```

#### Converting EMBL files to ASN format
The conversion of EMBL files to TBL format in [A2GB](https://github.com/PombertLab/A2GB) is a two step process. EMBL files are first converted to TBL format with [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl), then NCBI's [TBL2ASN](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) converts the later format to ASN.

##### Converting EMBL files to TBL format
[EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl) converts EMBL files to TBL format. This script requires a single tab-limited list of the locus tags and their predicted annotations. We can create this list by concatenating the tRNAs.annotations and rRNAs.annotations files generated [previously](https://github.com/PombertLab/A2GB#Creating-tab-delimited-lists-of-RNA-locus-tags-and-their-products) together with the curated list of proteins annotations (see [above](https://github.com/PombertLab/A2GB#curating-the-protein-annotations)). Alternatively, any tab-delimited list of locus_tags and their products can be used. 

Concatenating the annotations can be quickly performed with:
```
cat \
$ANNOT/tRNA.annotations \
$ANNOT/rRNA.annotations \
$ANNOT/proteins.annotations.curated \
> $ANNOT/verified_annotations.tsv
```

The conversion from EMBL to TBL can then be performed with:
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
-id		Desired institute ID [default: IITBIO]
-p		Tab-delimited list of locus_tags and their products
-embl		EMBL files to convert
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

The TBL files created should look like this:
```
head -n 25 `ls $ANNOT/splitGFF3/*.tbl | head -n 1`

>Feature chromosome_01
<2      >148    gene
                        locus_tag       HOP50_01g00010
<2      >148    mRNA
                        locus_tag       HOP50_01g00010
                        product hypothetical protein
                        protein_id      gnl|ITTBIO|HOP50_01g00010
                        transcript_id   gnl|ITTBIO|HOP50_01g00010_mRNA
<2      148     CDS
                        locus_tag       HOP50_01g00010
                        product hypothetical protein
                        protein_id      gnl|ITTBIO|HOP50_01g00010
                        transcript_id   gnl|ITTBIO|HOP50_01g00010_mRNA
<3043   >239    gene
                        locus_tag       HOP50_01g00020
<3043   >239    mRNA
                        locus_tag       HOP50_01g00020
                        product signal transduction histidine kinase
                        protein_id      gnl|ITTBIO|HOP50_01g00020
                        transcript_id   gnl|ITTBIO|HOP50_01g00020_mRNA
3043    239     CDS
                        locus_tag       HOP50_01g00020
                        product signal transduction histidine kinase
                        protein_id      gnl|ITTBIO|HOP50_01g00020
                        transcript_id   gnl|ITTBIO|HOP50_01g00020_mRNA
```

##### Converting TBL files to ASN format
Metadata must be included together with genome sequences during the submission process to NCBI. Although some of this metadata can be entered from the online submission form(s), it is often easier to add it beforehand while generating the ASN files. Metadata for genome submission includes taxonomic information about the source of the data being submitted, details about the sequencing experiments/computational analyses performed, and general information about the author(s) and institution(s) submitting the genomes.

###### Adding metadata to FASTA files
Taxonomic metadata can be added directly to the FASTA files. The list of modifiers that can be added directly to the FASTA definition lines can be found [here](https://www.ncbi.nlm.nih.gov/Sequin/modifiers.html). Mandatory modifiers include the organism name [organism=XXX] and its taxonomic lineage [lineage=XXX]. The latter can be found from the NCBI taxonomy database (see figure below).

<p align="center"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/lineage.png" alt="Lineage information from the NCBI Taxonomy database" width="1000"></p>

The script [add_metadata_to_fasta.pl](https://github.com/PombertLab/A2GB/blob/master/add_metadata_to_fasta.pl) can be used to add some of the most common modifiers to the FASTA definition lines. Alternatively, source qualifiers can also be entered directly from the [TBL2ASN](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) command line when creating the ASN files with the -j switch; e.g. -j '[organism=Chloropicon primus RCC138][strain=RCC138]'

To add metadata with [add_metadata_to_fasta.pl](https://github.com/PombertLab/A2GB/blob/master/add_metadata_to_fasta.pl), type:

```Bash
add_metadata_to_fasta.pl \
   -f $ANNOT/splitGFF3/*.fsa \
   -o 'Chloropicon primus RCC138' \
   -s RCC138 \
   -l 'cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;' \
   -g 1 \
   -c \
   -w 2
```

Options for [add_metadata_to_fasta.pl](https://github.com/PombertLab/A2GB/blob/master/add_metadata_to_fasta.pl) are:
```
-f (--fasta)		Specifies which FASTA files to add metadata to
-o (--organism)		Full organism name; e.g. 'Chloropicon primus RCC138'
-s (--strain)		Strain definition; e.g. RCC138
-l (--lineage)		NCBI taxonomic lineage; e.g. 'cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;'
-g (--gcode)		NCBI genetic code [Default: 1]
-m (--moltype)		NCBI moltype descriptor [Default: genomic]
-c (--chromosome)	Annotate contigs as chromosomes
-w (--width)		Character width for chromosome numbers [Default: 2] ## Adds padding zeroes if below threshold
```

Once modified, the FASTA definition lines should look like this:
```
head -n 1 $ANNOT/splitGFF3/*.fsa
==> /media/FatCat/user/raw_data/splitGFF3/chromosome_01.fsa <==
>chromosome_01 [organism=Chloropicon primus RCC138][strain=RCC138][lineage=cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;][gcode=1][moltype=genomic][chromosome=01]

==> /media/FatCat/user/raw_data/splitGFF3/chromosome_02.fsa <==
>chromosome_02 [organism=Chloropicon primus RCC138][strain=RCC138][lineage=cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;][gcode=1][moltype=genomic][chromosome=02]

==> /media/FatCat/user/raw_data/splitGFF3/chromosome_03.fsa <==
>chromosome_03 [organism=Chloropicon primus RCC138][strain=RCC138][lineage=cellular organisms; Eukaryota; Viridiplantae; Chlorophyta;][gcode=1][moltype=genomic][chromosome=03]
```

###### Generating a GenBank submission template
NCBI provides a simple web-based tool to generate the GenBank submission template file (template.sbt) required by TBL2ASN. To generate a template.sbt file, visit: https://submit.ncbi.nlm.nih.gov/genbank/template/submission/.

###### Creating structured comments
More information about NCBI's stuctured comments can be found [here](https://www.ncbi.nlm.nih.gov/genbank/structuredcomment/). Files containing stuctured comments are tab-delimited; assembly-data structured comments for genomes usually looks like this:
```
StructuredCommentPrefix	##Genome-Assembly-Data-START##
Assembly Method	SPAdes v. 3.13.0; Canu v. 1.8
Assembly Name	Version 1
Long Assembly Name	RCC138 version 1
Genome Coverage	345x
Sequencing Technology	Illumina MiSeq; Oxford Nanopore
```

###### Using TBL2ASN
[TBL2ASN](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) is a command-line program created by NCBI to automate the creation of sequence records for submission to [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). TBL2ASN will generate .sqn files to be used for the submission. More information about the structure and content of the TBL files generated by [EMBLtoTBL.pl](https://github.com/PombertLab/A2GB/blob/master/EMBLtoTBL.pl) and required as input for TBL2ASN can be found in NCBI's [Eukaryotic Genome Annotation Guide](https://www.ncbi.nlm.nih.gov/genbank/eukaryotic_genome_submission_annotation/).

To run [TBL2ASN](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) on the TBL files generated by [A2GB](https://github.com/PombertLab/A2GB) using a submisssion template and structured comments, type:
```Bash
tbl2asn \
   -t $ANNOT/template.sbt \
   -w $ANNOT/genome.cmt \
   -p $ANNOT/splitGFF3/ \
   -g \
   -M n \
   -Z $ANNOT/discrepancy.report \
   -H 12/31/2021
```

Options for [TBL2ASN](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) used above are:
```
-t	Template file (.sbt).
-w	File (.cmt) containing Genome Assembly structured comments.
-p	Path to the directory. If files are in the current directory -p ./ should be used.
-g	Genomic Product Set [T/F] ## Should be used for eukaryote genomes
-M n 	Master Genome Flags: n: Normal. Combines flags for genomes submissions (invokes FATAL calls when -Z discrep is included).
-Z	Runs the Discrepancy Report. Must supply an output file name. Only for annotated genome submissions.
-H	Desired date for data release
```

##### Checking for errors
[TBL2ASN](https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) generates two distinct types of error reports. The first consists of validation files with the file extension .val; one per FASTA file plus a summary titled errorsummary.val. The second report generated with -Z will be stored in the corresponding filename (discrepancy.report in the above command line).

Ideally, the .val files should be empty, indicating that no error has been found. We can check the size of our files easily with:
```Bash
ls -lh $ANNOT/splitGFF3/*.val

-rw-rw-r--. 1 jpombert jpombert 5.6K Dec  8 14:12 /media/FatCat/user/raw_data/splitGFF3/chromosome_01.val
-rw-rw-r--. 1 jpombert jpombert 4.7K Dec  8 14:12 /media/FatCat/user/raw_data/splitGFF3/chromosome_02.val
-rw-rw-r--. 1 jpombert jpombert 1.8K Dec  8 14:12 /media/FatCat/user/raw_data/splitGFF3/chromosome_03.val
...
-rw-rw-r--. 1 jpombert jpombert  409 Dec  8 14:12 /media/FatCat/user/raw_data/splitGFF3/errorsummary.val

```

In the above example, the file sizes are not zero, which means that errors have been detected. In similar situations, the errorsummary.val file should look like this:
```Bash
cat $ANNOT/splitGFF3/errorsummary.val

    34 ERROR:   SEQ_FEAT.GeneXrefStrandProblem
    35 ERROR:   SEQ_FEAT.InternalStop
     1 ERROR:   SEQ_FEAT.NoStop
    40 ERROR:   SEQ_FEAT.SeqLocOrder
    35 ERROR:   SEQ_INST.StopInProtein
    70 WARNING: SEQ_FEAT.NotSpliceConsensusAcceptor
    74 WARNING: SEQ_FEAT.NotSpliceConsensusDonor
    54 WARNING: SEQ_FEAT.PartialProblem
    15 WARNING: SEQ_FEAT.ShortExon
     2 INFO:    SEQ_FEAT.PartialProblem
```

Errors will vary per file, obviously, but the content of a .val file should look like:
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

