#!/usr/bin/env python3
__author__      = "Theo Cole"
__email__       = "theo.cole@nhs.net"
__description__ = """
    Class to manage generation of a BED file containing coordinates for all
    exons included in a user-provided bed file, or if no bed file provided then
    all exons for all known canonical RefSeq transcripts (using the UCSC table
    API).
"""
import os
import re
import sys

import sqlalchemy
import pybedtools
import pandas as pd


CHROMOSOMES = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
               'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
               'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22',
               'chrX', 'chrY', 'chrMT']


class BEDGenerator():
    def __init__(self, genome, bedfile, output_dir, output_name):
        self.genome = genome
        self.bedfile = bedfile
        self.output_dir = output_dir
        self.output_name = output_name
        self.outfile = os.path.join(
            self.output_dir, 
            f"{self.output_name}.{self.genome}.bed"
        )

    def run(self):
        df_bedfile = self.create_df_bedfile(
            df_canonical_exons=self.create_canonical_exons_df(
                df_canonical=self.get_canonical_df(),
                df_exons=self.get_exons_df()
            )
        )
        df_bedfile.to_csv(
            path_or_buf=self.outfile,
            sep='\t',
            header=False,
            index=False
        )
        return df_bedfile

    def get_canonical_df(self):
        # only use canonical genes if user doesn't specify bed
        # maybe function name is a misnomer....
        if self.bedfile:
            table = 'knownGene'
            fkid = 'name'
        else:
            table = 'knownCanonical'
            fkid = 'transcript'

        sys.stdout.write(f"Connecting to UCSC {self.genome} {table}...")
        sys.stdout.flush()

        # sql connection, requires join to fetch refseq ID for UCSC canonical
        select = f"SELECT kgXref.refseq FROM {table}"
        join = f"LEFT JOIN kgXref ON kg.{fkid} = kgXref.kgID"
        query = f"{select} AS kg {join}"

        db_connection = self.set_up_db_conn()
        df = pd.read_sql(
            query,
            con=db_connection
        )

        # filter transcripts out if necessary
        if self.bedfile:
            df = self.filter_bedfile_regions(df)

        # only keep mRNA
        df = df[df.refseq.str.startswith('NM')]
        df = df.set_index('refseq').dropna()

        sys.stdout.write('OK!\n')
        return df

    def filter_bedfile_regions(self, df_all_transcripts):
        # fetch a list of transcripts in arg bedfile, filter df based on this
        bedf_tx = self.get_bedfile_transcripts()
        df = df_all_transcripts[df_all_transcripts.refseq.isin(bedf_tx)]

        # identify transcripts in bedfile but not UCSC to still fetch refseq
        missing_transcripts = list(set(bedf_tx) - set(df.refseq))
        df_missing = pd.DataFrame(missing_transcripts, columns=['refseq'])
        df = pd.concat(
            [df, df_missing], 
            sort=False
        ).reset_index(drop=True).fillna('missing')
        return df

    def get_bedfile_transcripts(self):
        # load the bedfile, read the name from each row
        bedfile = pybedtools.BedTool(self.bedfile)
        tx_names = list(set([
            re.findall(r'N[RM]_\d+', tx)[0] 
            for tx in pd.core.common.flatten(
                # split if there are mutliple transcripts per row
                [x.name.split(',') for x in bedfile]
            )
        ]))
        return tx_names

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

        sys.stdout.write('OK!\n')
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

        df_canonical_exons = df_canonical.join(df_exons, lsuffix='_ucsc')
        df_canonical_exons = df_canonical_exons.dropna()

        sys.stdout.write('OK!\n')
        return df_canonical_exons

    def create_df_bedfile(self, df_canonical_exons):
        sys.stdout.write(f"Fitting exon coordinates to transcripts...")
        sys.stdout.flush()

        #bytes to string columns
        decoded_s = df_canonical_exons.exonStarts.str.decode('utf-8')
        decoded_e = df_canonical_exons.exonEnds.str.decode('utf-8')
        df_canonical_exons.exonStarts = decoded_s
        df_canonical_exons.exonEnds = decoded_e

        df_bedfile = pd.concat(
            [
                pd.DataFrame(
                    data={
                        'chrom': row.chrom, 
                        'chromStart': map(
                            int, 
                            row.exonStarts.strip(',').split(',',)
                        ),
                        'chromEnd': map(
                            int,
                            row.exonEnds.strip(',').split(',')
                        ),
                        'name': map(
                            # two params to map: 
                            # 1) lambda func to append exon no. to tx name
                            lambda x: row.name + '_exon' + str(x + 1), 
                            # 2) range of exon numbers, in reverse if - strand
                            range(
                                len(row.exonStarts.strip(',').split(','))
                            ) if row.strand == '+'
                            else range(
                                len(row.exonStarts.strip(',').split(','))
                            )[::-1]
                        ), 
                        'strand': row.strand
                    },

                ) 
                for idx, row in df_canonical_exons.iterrows()
            ]
        ).drop_duplicates().reset_index(drop=True)

        # sorting by chrom
        chromosome_cat = pd.api.types.CategoricalDtype(
            categories=CHROMOSOMES,
            ordered=True
        )
        df_bedfile.chrom = df_bedfile.chrom.astype(chromosome_cat)
        df_bedfile = df_bedfile.sort_values(by='chrom')
        sys.stdout.write('OK!\n')
        sys.stdout.flush()


        # check if the users bedfile is missing any transcripts from result
        if self.bedfile:
            sys.stdout.write(f"Checking {self.bedfile} for transcripts...\n")
            found_transcripts = set(
                '_'.join(nm.split('_')[:2]) for nm in df_bedfile.name
            )
            desired_transcripts = set(self.get_bedfile_transcripts())
            missing_transcripts = desired_transcripts - found_transcripts
            if missing_transcripts:
                sys.stdout.write('Missing the following transcripts:\n')
                for missing_transcript in missing_transcripts:
                    sys.stdout.write(missing_transcript + '\n')
            else:
                sys.stdout.write('All transcripts from bedfile matched.\n')
        return df_bedfile

