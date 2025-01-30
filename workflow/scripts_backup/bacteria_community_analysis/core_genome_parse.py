#!/usr/bin/env python

import sys
import argparse
import re
import traceback
from itertools import groupby
from Bio import SeqIO
import pandas as pd
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


__author__ = "Malick Ndiaye"
__title__ = "Parse Core Genome"


def get_args():
    """Parse command line arguments"""
    desc = (
        """
        This script identifiees the core genome
        """
    )
    epi = """"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # input files
        parser.add_argument("-i",
            "--inog", action="store",
             help='orthogroups table'
        )
        parser.add_argument("-m",
            "--motupan", action="store",
             help='motupan table', nargs='+'
        )
        parser.add_argument("-g",
            "--gene", action="store",
             help='gene table'
        )

       # output files
        parser.add_argument("-o",
            "--outtab", action="store",
             help='output gene file'
        )
        parser.add_argument("-c",
            "--outcore", action="store",
             help='output gene file of core genome'
        )

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occurred with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return(parser.parse_args())

def concat_motupan(motu_list):
    """Concatenate motupan tables"""
    #use pd.concat to concatenate the dataframes
    # skip first 15 rows of each table
    df_list = []
    for motu in motu_list:
        df = pd.read_csv(motu, sep='\t', skiprows=16)

        # the column "genomes" is a list separated by ";", split it into multiple rows
        df= df.assign(genome=df['genomes'].str.split(';').tolist())
        df = df.explode('genome')

        # remove column "genomes"
        df = df.drop(columns=['genomes'])
        
        if "Lactobacillus" in motu:
            df.to_csv("Lactobacillus.csv", sep='\t', index=False)

        df_list.append(df)
    
    df_concat = pd.concat(df_list)

    # remove trailing and leading whitespaces all columns
    df_concat = df_concat.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    # retain only columns: "genomes", trait_name, mean_copy_per_genome,  type
    df_concat = df_concat[['genome', 'trait_name', 'mean_copy_per_genome', 'type']]
    
    return(df_concat)

def parse_og_file(og_file, df_motu):
    # read orthogroups file
    df_og = pd.read_csv(og_file, sep='\t')

    # remove trailing and leading whitespaces all columns
    df_og = df_og.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    # merge df_og and df_motu by genome and OG=trait_name
    df_merge = pd.merge(df_og, df_motu, left_on=['genome', 'OG'], right_on=['genome', 'trait_name'], how='left')

    # remove trailing and leading whitespaces all columns
    df_merge = df_merge.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    return(df_merge)


def parse_instrain_genes(gene_file, df_merge):
    # read gene file
    df_gene = pd.read_csv(gene_file, sep='\t')

    # merge df_gene and df_merge by gene
    df_merge2 = pd.merge(df_gene, df_merge, left_on=['gene'], right_on=['gene'], how='left')

    # remove colums: trait_name, mean_copy_per_genome
    df_merge2 = df_merge2.drop(columns=['trait_name', 'mean_copy_per_genome'])

    # move colmns genome, id, genus and type to the beginning
    # cols = list(df_merge2)
    # cols.insert(0, cols.pop(cols.index('genome')))
    # cols.insert(1, cols.pop(cols.index('id')))
    # cols.insert(2, cols.pop(cols.index('genus')))
    # cols.insert(3, cols.pop(cols.index('type')))
    # df_merge2 = df_merge2.loc[:, cols]

    
    return(df_merge2)


def main():
    args = get_args()

    # concatenate motupan tables
    df_motu = concat_motupan(args.motupan)

    # parse orthogroups file
    df_merge = parse_og_file(args.inog, df_motu)
    # parse instrain genes file
    df_merge2 = parse_instrain_genes(args.gene, df_merge)

    # get dataframe with only core genome
    df_core = df_merge2[df_merge2['type'] == 'core']

    # write output
    df_merge.to_csv("merge.csv", sep='\t', header=True, index=False)
    df_core.to_csv(args.outcore, sep='\t', header=True, index=False)
    df_merge2.to_csv(args.outtab, sep='\t', header=True, index=False)

if __name__ == "__main__":
    main()


