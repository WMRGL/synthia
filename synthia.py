#!/usr/bin/env python3
__author__      = "Theo Cole"
__email__       = "theo.cole@nhs.net"
__description__ = """
    A program which takes a reference genome (hg19 or hg38) name as input as
    well as (optionally) a bed file which has a 'name' column describing the
    RefSeq ID of the bed file regions. If no bed file is given, then the whole
    exome will be used as default. A VCF will then be generated containing one
    random variant per exon. This VCF is then used in combination with a ref
    genome to generate synthetic FASTQs which can be used to test if a pipeline
    can pick up all the known synthetic variants (exome or in the desired
    panel).
"""
import argparse
import configparser
import os
import sys

from generate_vcf import bed_generator
from generate_vcf import vcf_generator


def main(args, config):
    sys.stdout.write(f'Creating bed file...\n')
    sys.stdout.flush()
    bed_gen = bed_generator.BEDGenerator(
        genome=args.genome,
        bedfile=args.bedfile,
        output_dir=args.output_dir,
        output_name=args.name
    )
    df_bedfile = bed_gen.run()
    sys.stdout.write(f'Saved bed file to {bed_gen.outfile}\n')

    sys.stdout.write(f'Creating VCF file...\n')
    sys.stdout.flush()

    vcf_gen = vcf_generator.VCFGenerator(
        df_bed=df_bedfile,
        process_count=args.process_count,
        genome=args.genome,
        genome_file=config['GENOME_RESOURCES'][args.genome],
        output_dir=args.output_dir,
        output_name=args.name
    )
    df_vcf = vcf_gen.run()
    sys.stdout.write(f'Saved vcf file to {vcf_gen.outfile}\n')
    sys.stdout.flush()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__description__)
    config = configparser.ConfigParser()
    config.read('config.ini')

    parser.add_argument(
        '-g', '--genome',
        required=True,
        choices=['hg19', 'hg38'],
        type=str,
        help="The reference genome used to generate synthetic variants."
    )
    parser.add_argument(
        '-o', '--output_dir',
        required=True,
        type=str,
        help="The output directory for BED, VCF, and FASTQ files."
    )
    parser.add_argument(
        '-n', '--name',
        type=str,
        default="allExons",
        help="""
            Name of panel (optional; defaults to 'allExons') for output
            prefix.
        """
    )
    parser.add_argument(
        '-b', '--bedfile',
        type=str,
        help="""
            A bedfile which will be used to specify desired regions for exonic
            variants.
        """
    )
    parser.add_argument(
        '-p', '--process_count',
        type=int,
        default=1,
        help="""
            The number of processes which should be used to generate the VCF
            and FASTQ files (optional; defaults to 1).
        """
    )

    args = parser.parse_args()
    # check genome file is indexed
    if not os.path.exists(config['GENOME_RESOURCES'][args.genome] + '.fai'):
        raise Exception(f"Cannot find index for {args.genome}.")

    main(args, config)
