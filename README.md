# A2GB – Annotations to GenBank

[A2GB](https://github.com/PombertLab/A2GB/) is a pipeline that will transform the genome annotation files exported from [Apollo](https://genomearchitect.readthedocs.io/en/latest/) (fomerly known as [WebApollo](http://gmod.org/wiki/WebApollo)) in preparation for sequence submission to NCBI’s [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). While its primary goal is to submit annotations to GenBank for the generation of accession numbers, the conversion of file formats will permit the use of different tools in downstream analyses. 

In a stepwise approach, the A2GB pipeline converts annotation files from GFF3 -> EMBL -> TBL -. ASN. Each format will be useful for diagnostic quality checks of the annotations or become the springboard for other analyses, such as protein function prediction.  

Furthermore, A2GB acts as a guide to prepare sequence submissions according to NCBI’s guidelines, including project registration with [BioSample](https://www.ncbi.nlm.nih.gov/biosample) and [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) for the generation of locus_tag prefixes. Upon acceptance from NCBI, these essential steps will ensure that sequence data will be made publicly available through GenBank and other member databases within the [International Nucleotide Sequence Database Collaboration](http://www.insdc.org/).

## Table of contents
* [Introduction](#introduction)
* [Requirements](#requirements)
* [A2GB workflow](#A2GB-workflow)
   * [Exporting annotations from Apollo](#Exporting-annotations-from-Apollo)
   *	[Splitting Apollo GFF3 files](#Splitting-Apollo-GFF3-files)
   *	[Converting the GFF3 files to EMBL format]
   *	[Function prediction]
        *	[Predicting functions with InterProScan 5]
        *	[Downloading the SwissProt/UniProt databases]
        *	[Creating tab-delimited lists of sequences in the SwissProt/UniProt databases]
        *	[Running BLAST searches against SwissProt/UniProt]
        *	[Generating a list of all proteins queried]
        *	[Parsing the result of InterProScan 5 and SwissProt/UniProt searches]
        *	[Curating the annotations]
   *	[Adding taxonomic info to FASTA files]
   *	[Converting EMBL files to TBL format]
   *	[Generating a template.sbt file per genome]
   *	[Creating a structure comments file (genome.cmt)]
   *	[Converting TBL files to ASN format]
   *	[Checking for errors]
   *	[Submitting ASN file to GenBank]
*	[Miscellaneous] 
*	[References]

## Introduction

Various analyses contribute to assigning biological interpretation to a DNA sequence. The goal of genome annotation is to identify the location and function of a genome's encoded features, particularly of protein-coding and RNA-coding genes. Thus, generating the best descriptions of encoded features (i.e. annotations) is paramount for any genomic sequencing project that intends to identify these encoded features and attribute biological meaning to them.

[Apollo](https://genomearchitect.readthedocs.io/en/latest/) provides sequencing projects with the tools for gene prediction via evidence-based annotation. This is accomplished with the automatic generation of sequence features which can be refined through expert user curation. Exporting and utilizing these  annotations is the preliminary step to assigning user-curated biological function to predictions.

To extend the value of these efforts to the broader scientific community, annotation should be made available to the public in widely accessible biological databases. Upon successful submission, GenBank will assign accession numbers to submitted data to act as unique identifiers.  To achieve this, careful adherence to NCBI guidelines for genome submission is essential. This pipeline is equipped with a series of checks to assess the quality of annotations and minimize errors.

The A2GB pipeline will:
1)	Reformat the annotations from Apollo for deposition into the NCBI database. 
2)	Run sequence searches against [UniProt](https://www.uniprot.org/)'s SwissProt/TrEMBL and [InterPro](https://www.ebi.ac.uk/interpro/) databases for protein function prediction. 
3)	Run intermittent checks to assess the quality of annotations.

## Requirements
- Unix/Linux, MacOS X, or Miscrosoft's [WSL2](https://docs.microsoft.com/en-us/windows/wsl/compare-versions)
-	[Perl](https://www.perl.org/) 5
-	[InterProScan](https://github.com/ebi-pf-team/interproscan) (latest version)
-	[Artemis](http://sanger-pathogens.github.io/Artemis/Artemis/) (18.0.0+)
- [Apollo](https://genomearchitect.readthedocs.io/en/latest/) (2.5.0+)
- [RNAmmer](https://services.healthtech.dtu.dk/software.php) (1.2+)
- [tRNAscan-SE](http://lowelab.ucsc.edu/tRNAscan-SE/) (2.0+)

#### A2GB workflow
##### Exporting annotations from Apollo
After genomic annotations are completed in Apollo, export the curated annotations.  Begin by selecting the 'Ref Sequence' tab. 
<p align="left"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Apollo.png" alt="How to export Apollo annotations" width="600"></p>

Then, select Export -> select GFF3; select ALL; select GFF3 with FASTA; click Export. The file created will be called Annotations.gff3.gz.
<p align="left"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Apollo2.png" alt="How to export Apollo annotations" width="400"></p>

We recommend exporting only protein features (CDS) from Apollo. Altough rRNAs and tRNAs inferences (e.g. from RNAmmer and tRNAscan-SE, respectively) can be added to Apollo, the process of exporting those back is finicky and prone to errors (some rRNAs/tRNAs appear to be missing when exporting features from Apollo 2.5.0).

First, let's create a directory to store annotations:

```Bash
export ANNOT=/media/FatCat/user/raw_data     ## Replace /media/FatCat/user/raw_data by desired annotation directory
mkdir $ANNOT; mv Annotations.gff3.gz $ANNOT/ ## Create directory; then move Annotations.gff3.gz into it
cd $ANNOT/; gunzip Annotations.gff3.gz       ## Decompress the GZIP file
```

Second, let's predict ribosomal RNAs with RNAmmer; then convert the annotations to GFF3 format:

```Bash
mkdir $ANNOT/RNAmmer/
run_RNAmmer.pl -f *.fasta -d $ANNOT/RNAmmer/
RNAmmer_to_GFF3.pl -g RNAmmer/*.gff2 -d RNAmmer/
```


Third , let's predict transfer RNAs with tRNAscan-SE; then convert the annotations to GFF3 format:

```Bash
mkdir $ANNOT/RNAmmer/tRNAscan

```

Fourth, let's concatenate the tRNA, rRNA and CDS GFF3 annotations from RNAmmer, tRNAscan-SE, and Apollo:

```Bash
cat RCC* Annotations.gff3 > all_genes.gff3  ## We concatenate Apollo's GFF3 file last as it includes sequence data
```

##### Splitting Apollo GFF3 files
