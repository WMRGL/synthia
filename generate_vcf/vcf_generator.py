#!/usr/bin/env python3
__author__      = "Theo Cole"
__email__       = "theo.cole@nhs.net"
__description__ = """
    Class which takes a BED-like dataframe, then proceeds to create a VCF by
    looking up coordinates in a reference genome, creating a random variant
    within the returned sequence.
"""
import multiprocessing
import os
import random
import sys

import pysam


VARIANT_OPTIONS = {
    'A': ['T', 'C', 'G'], 'a': ['t', 'c', 'g'],
    'T': ['A', 'C', 'G'], 't': ['a', 'c', 'g'],
    'C': ['T', 'A', 'G'], 'c': ['t', 'a', 'g'],
    'G': ['T', 'C', 'A'], 'g': ['t', 'c', 'a']
}


class VCFGenerator():
    def __init__(self, df_bed, process_count, genome, genome_file, output_dir,
                 output_name):
        self.df_bed = df_bed
        self.process_count = process_count
        self.genome = genome
        self.genome_file = genome_file
        self.outfile = os.path.join(
            output_dir, 
            f"{output_name}.{genome}.bed"
        )

    def run(self):
        # get sequence for each bed record in parallel
        process_pool = multiprocessing.Pool(self.process_count)

        regions = self.df_bed.to_dict('records')

        sys.stdout.write('Creating variants for bed file regions...')
        sys.stdout.flush()
        variants = process_pool.map(self.create_variants, regions)
        sys.stdout.write('OK!')
        sys.stdout.flush()
        print(variants)

    def create_variants(self, region):
        sequence = self.lookup_sequence(region)
        variant = self.create_random_variant(region, sequence)


    def lookup_sequence(self, region):
        # look up the DNA sequence from the ref given in config file
        dna_lookup = pysam.faidx(
            self.genome_file,
            f"{region['chrom']}:{region['chromStart']}-{region['chromEnd']}"
        )
        sequence = ''.join(dna_lookup.split('\n')[1:])
        return sequence

    def create_random_variant(self, region, sequence):
        # don't want to pull an unknown base from ref
        vcf_ref = 'N'
        while vcf_ref != 'N':
            exon_pos = random.randrange(len(sequence))
            vcf_pos = region['chromStart'] + exon_pos
            vcf_ref = sequence[exon_pos]

        # select random base that is not the same; account for lowercase in ref
        vcf_alt = random.choice(VARIANT_OPTIONS[vcf_ref])

        variant = {
            '#CHROM': region[0],
            'POS': vcf_pos,
            'REF': vcf_ref,
            'ALT': vcf_alt,
            'QUAL': '.',
            'FILTER': 'PASS',
            'INFO': '.',
            'FORMAT': 'GT',
            'SAMP001': '0|1'
        }
        return variant
