#!/usr/bin/env python3
'''
this script generates a bedfile which contains each exon for every known
canonical exon in refseq, it requires two files generated using the UCSC table
browser:
1) UCSC knownCanonical 
    - track: UCSC Genes
    - table: knownCanonical
    - select: chrom, refseq, kgXref, 
2) RefSeq all
    - track: NCBI RefSeq
    - table: RefSeq All
    - select: name, chrom, strand, exonStarts, exonEnds
'''
import argparse
import numpy as np
import pandas as pd


class VCFGenerator():
    def __init__(self):
        pass

    def run(self):
        df_canonical = pd.read_csv(self.canonical_transcripts, sep='\t')
        df_canonical = df_canonical.set_index('hg19.kgXref.refseq')

        df_exons = pd.read_csv(self.canonical_exons, sep='\t')
        df_exons = df_exons.set_index('refSeqID')
        df_exons['refSeqID'] = df_exons['#name'].str.split('.', expand=True)[0]


    def create_canonical_exons_df():
        df_canonical_exons = df_canonical.join(df_exons)
        df_canonical_exons = df_canonical_exons.dropna()
        df_exons_for_bed = pd.concat(
            [
                pd.DataFrame(
                    data={
                        'chrom': row['#hg19.knownCanonical.chrom'], 
                        'chromStart': map(int, row.exonStarts.strip(',').split(',',)),
                        'chromEnd': map(int, row.exonEnds.strip(',').split(',')),
                        'name': map(
                            lambda x: idx + '_exon' + str(x + 1), 
                            range(len(row.exonStarts.strip(',').split(',')))
                        ), 
                        'strand': row.strand}) 
                for idx, row in df_canonical_exons.iterrows()
            ]
        ).sort_values(by=["chrom", "chromStart"]).reset_index(drop=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Read in files for I/O.")
    vcf_generator = VCFGenerator()
    vcf_generator.run()
