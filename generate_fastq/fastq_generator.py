#!/usr/bin/env python3
__author__      = "Theo Cole"
__email__       = "theo.cole@nhs.net"
__description__ = """
    Class which takes a VCF-like dataframe, then creates a number of files
    based on the number of processes the synthia program is set to run with.
    These VCFs are then passed to the DWGSIM program along with the relevant
    reference file to generate Illumina-like FASTQ files.
"""


class FASTQGenerator():
    def __init__(self, df_vcf, process_count, genome, genome_file, output_dir,
                 output_name)
