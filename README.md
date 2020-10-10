# A2GB – Annotations to GenBank

[A2GB](https://github.com/PombertLab/A2GB/) is a pipeline that will transform the genome annotation files exported from [Apollo](https://genomearchitect.readthedocs.io/en/latest/) (fomerly known as [WebApollo](http://gmod.org/wiki/WebApollo)) in preparation for sequence submission to NCBI’s [GenBank](https://www.ncbi.nlm.nih.gov/genbank/). While its primary goal is to submit annotations to GenBank for the generation of accession numbers, the conversion of file formats will permit the use of different tools in downstream analyses. 

In a stepwise approach, the A2GB pipeline converts annotation files from GFF3 -> EMBL -> TBL -. ASN. Each format will be useful for diagnostic quality checks of the annotations or become the springboard for other analyses, such as protein function prediction.  

Furthermore, A2GB acts as a guide to prepare your sequence submission according to NCBI’s guidelines, including registration of your project with [BioSample](https://www.ncbi.nlm.nih.gov/biosample) and [BioProject](https://www.ncbi.nlm.nih.gov/bioproject) for the generation of locus_tag prefixes. Upon acceptance from NCBI, these essential steps will ensure that your sequence data will be made publicly available through GenBank and other member databases within the [International Nucleotide Sequence Database Collaboration](http://www.insdc.org/).

## Table of contents
* [Introduction](#introduction)
* [Requirements](#requirements)
  * [A2GB workflow](#A2GB-workflow)
    * [Exporting annotations from Apollo](#Exporting-annotations-from-Apollo)
    *	[Splitting Apollo GFF3 files]
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
-	[Artemis](https://www.sanger.ac.uk/tool/artemis/) (latest version)

#### A2GB workflow
##### Exporting annotations from Apollo
After genomic annotations are completed in Apollo, export the curated annotations.  Begin by selecting the 'Ref Sequence' tab.
<p align="left"><img src="https://github.com/PombertLab/A2GB/blob/master/Misc/Apollo.png" alt="How to export Apollo annotations" width="600"></p>


