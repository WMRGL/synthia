#!/usr/bin/env python3
__author__      = "Theo Cole"
__email__       = "theo.cole@nhs.net"
__description__ = """
    A script which takes a reference genome (hg19 or hg38) name as input as
    well as (optionally) a bed file which has a 'name' column describing the
    RefSeq ID of the bed file regions. If no bed file is given, then the whole
    exome will be used as default.

    Only canonical transcripts are recognised as genes when generating exons.
    These transcripts and exons are pulled from the UCSC public MySQL databases
    'knownCanonical' to fetch canonical genes and associated RefSeq ID, and 
"""
import argparse
import sys

import sqlalchemy
import pandas as pd


class VCFGenerator():
    def __init__(self, genome):
        self.genome = genome

    def run(self):
        df_bedfile = self.create_df_bedfile(
            df_canonical_exons=self.create_canonical_exons_df(
                df_canonical=self.get_canonical_df(),
                df_exons=self.get_exons_df()
            )
        )
        print(df_bedfile)

    def get_canonical_df(self):
        sys.stdout.write(f"Connecting to UCSC {self.genome} knownCanonical...")
        sys.stdout.flush()

        db_connection = self.set_up_db_conn()
        df = pd.read_sql(
            'SELECT * FROM knownCanonical'
            ' LEFT JOIN kgXref ON knownCanonical.transcript = kgXref.kgID', 
            con=db_connection
        )
        df = df.set_index('refseq').dropna()

        sys.stdout.write('\tOK!\n')
        return df

    def get_exons_df(self):
        sys.stdout.write(f"Connecting to UCSC {self.genome} ncbiRefSeq...")
        sys.stdout.flush()

        db_connection = self.set_up_db_conn()
        df = pd.read_sql(
            'SELECT * FROM ncbiRefSeq',
            con=db_connection
        )
        df['refseq'] = df['name'].str.split('.', expand=True)[0]
        df= df.set_index('refseq')

        # bytes to string
        df.exonStarts = df.exonStarts.str.decode('utf-8')
        df.exonEnds = df.exonEnds.str.decode('utf-8')

        sys.stdout.write('\t\tOK!\n')
        return df

    def set_up_db_conn(self):
        protocol = 'mysql+pymysql'
        user = 'genome'
        host = 'genome-mysql.soe.ucsc.edu'
        db_connection_addr = f"{protocol}://{user}@{host}/{self.genome}"
        db_connection = sqlalchemy.create_engine(db_connection_addr)
        return db_connection

    def create_canonical_exons_df(self, df_canonical, df_exons):
        sys.stdout.write(f"Merging canonical and exon dataframes...")
        sys.stdout.flush()

        df_canonical_exons = df_canonical.join(df_exons, rsuffix='ncbi')
        df_canonical_exons = df_canonical_exons.dropna()

        sys.stdout.write('\tOK!\n')
        return df_canonical_exons

    def create_df_bedfile(self, df_canonical_exons):
        sys.stdout.write(f"Fitting exon coordinates to transcripts...")
        sys.stdout.flush()
        df_bedfile = pd.concat(
            [
                pd.DataFrame(
                    data={
                        'chrom': row['chrom'], 
                        'chromStart': map(
                            int, 
                            row.exonStarts.strip(',').split(',',)
                        ),
                        'chromEnd': map(
                            int,
                            row.exonEnds.strip(',').split(',')
                        ),
                        'name': map(
                            lambda x: row['name'] + '_exon' + str(x + 1), 
                            range(len(row.exonStarts.strip(',').split(',')))
                        ), 
                        'strand': row.strand
                    },
                ) 
                for idx, row in df_canonical_exons.iterrows()
            ]
        ).sort_values(by=["chrom", "chromStart"]).reset_index(drop=True)
        sys.stdout.write('\tOK!\n')
        return df_bedfile

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument(
        'genome',
        choices=['hg19', 'hg38'],
        type=str,
        help="The reference genome used to generate synthetic variants."
    )
    parser.add_argument(
        '-b', '--bedfile',
        type=str,
        help="""
            A bedfile which will be used to specify desired regions for exonic
            variants.
        """
    )
    args = parser.parse_args()

    vcf_generator = VCFGenerator(
        genome=args.genome
    )
    vcf_generator.run()
