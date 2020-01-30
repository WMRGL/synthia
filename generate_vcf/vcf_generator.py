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
import sys

import pysam

from generate_vcf.utils import MultiprocessCounter


seq_counter = MultiprocessCounter(0)
var_counter = MultiprocessCounter(0)


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
        # set later
        self.total_regions = None  

    def run(self):
        # get sequence for each bed record in parallel
        process_pool = multiprocessing.Pool(self.process_count)

        regions = self.df_bed.to_dict('records')
        self.total_regions = len(regions)

        count = f"{seq_counter.value}/{self.total_regions}"

        variants = process_pool.map(self.create_variants, regions)

    def create_variants(self, region):
        sys.stdout.write(f'Fetching sequence for all regions...{count}\r')
        sys.stdout.flush()
        sequence = lookup_sequence(region)
        sys.stdout.write('\n')
        sys.stdout.flush()

        sys.stdout.write(f'Creating variant for all regions...{count}\r')
        sys.stdout.flush()
        variant = create_random_variant(region, sequence)
        sys.stdout.write('\n')
        sys.stdout.flush()

    def create_random_variant(self, region, sequence):
        raise NotImplementedError()

    def lookup_sequence(self, region):
        # look up the DNA sequence from the ref given in config file
        dna_lookup = pysam.faidx(
            self.genome_file,
            f"{region['chrom']}:{region['chromStart']}-{region['chromEnd']}"
        )
        sequence = ''.join(dna_lookup.split('\n')[1:])
        seq_counter.increment()
        count = f"{seq_counter.value}/{self.total_regions}"
        sys.stdout.write(f'Fetching sequence for all regions...{count}\r')
        return sequence
