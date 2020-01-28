# synthia
Generate synthetic FASTQ from a reference including 1 variant per exon for QA
purposes.

This is a two-step process:
1. *generate vcf*: a bedfile is created from the UCSC knownCanonical genes
   list, optionally filtered by a specified bedfile containin a list of RefSeq
   IDs. From this bedfile, a VCF is generated contianing 1 random SNP in each
   exon.
2. *synthesise fastq*: a FASTQ file is generated from the reference genome used
    in Step 1, along with the VCF of desired variants which are introduced to
    the FASTQ files.

## generate_vcf.py
```
usage: generate_vcf.py [-h] [-b BEDFILE] {hg19,hg38}

A script which takes a reference genome (hg19 or hg38) name as input as well
as (optionally) a bed file which has a 'name' column describing the RefSeq ID
of the bed file regions. If no bed file is given, then the whole exome will be
used as default. Only canonical transcripts are recognised as genes when
generating exons. These transcripts and exons are pulled from the UCSC public
MySQL databases 'knownCanonical' to fetch canonical genes and associated
RefSeq ID, and

positional arguments:
  {hg19,hg38}           The reference genome used to generate synthetic
                        variants.

optional arguments:
  -h, --help            show this help message and exit
  -b BEDFILE, --bedfile BEDFILE
                        A bedfile which will be used to specify desired
                        regions for exonic variants.

```
