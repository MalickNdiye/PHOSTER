#!/usr/bin/env python

import sys
import argparse
import re
import traceback
from itertools import groupby
import pandas as pd
import os


__author__ = "Malick Ndiaye"
__title__ = "Parse Instrain"


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
            "--instrain", action="store",
             help='instrain output', nargs='+'
        )
        parser.add_argument("-g",
            "--genes", action="store",
             help='genes table', nargs='+'
        )
        parser.add_argument("-t",
            "--taxa", action="store",
             help='taxonomic table'
        )

       # output files
        parser.add_argument("-o",
            "--outdir", action="store",
             help='output directroy'
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

def get_filelist(dir, pattern):
    for file in os.listdir(dir):
        if re.search(pattern, file):
            path = os.path.join(dir, file)
    return(path)

def get_sample(file):
    sample=os.path.basename(file).split('_')[0]
    return(sample)

def concat_df(file_list, taxa_file):
    # concat dataframes and add sample name
    df_list = []
    for file in file_list:
        sample = get_sample(file)
        # if "mapping" in filename, skip fitst row
        if "mapping" in file:
            df = pd.read_csv(file, sep='\t', skiprows=1)
        else:
            df = pd.read_csv(file, sep='\t')

        df['sample'] = sample
        # move column "sample" to first column
        cols = df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df = df[cols]

        df_list.append(df)
    df_concat = pd.concat(df_list)

    # if genome in df colums, add column genus, id from taxa table. the genome will correspond to "bin Id" in taxa table
    if 'iRep' in df_concat.columns:
        taxa = pd.read_csv(taxa_file, sep='\t')
        df_concat = pd.merge(df_concat, taxa[['Bin Id', 'genus']], how='left', left_on='genome', right_on='Bin Id')
        df_concat = df_concat.drop('Bin Id', axis=1)
        # move genome, id, genus columns after sample column
        cols = df_concat.columns.tolist()
        cols = cols[:1] + cols[-3:] + cols[1:-3]
        df_concat = df_concat[cols]
        
    return(df_concat)

def filter_genome_info(df):
    # filter for row that have the column "breadth" >=0.5 and "coverge_median" >=5
    df = df[(df['breadth'] >= 0.5) & (df['coverage_median'] >= 5)]

    # for genome in each sample, caluclate "actual length" as length*breath
    df['actual_length'] = df['length']*df['breadth']
    # now calulate "norm_reads" as filtered_read_pair_count/actual_length
    df['norm_reads'] = df['filtered_read_pair_count']/df['actual_length']
    # now calulate "frequency" as norm_reads/sum(norm_reads) within each sample
    df['frequency'] = df.groupby('sample')['norm_reads'].apply(lambda x: x/x.sum())

def main():
    args = get_args()

    # get list of instrain genome info files
    instrain_output=[os.path.join(i, 'output') for i in args.instrain]
    instrain_genome_info = [get_filelist(i, 'genome_info.tsv') for i in instrain_output]
    

    # get list of instrain mapping info files
    instrain_mapping_info = [get_filelist(i, 'mapping_info.tsv') for i in instrain_output]

    all_genome_info = concat_df(instrain_genome_info, args.taxa)
    print(all_genome_info.head())
    print(all_genome_info.columns)
    all_genome_info_filtered = filter_genome_info(all_genome_info)

    all_mapping_info = concat_df(instrain_mapping_info, args.taxa)
    # retain only colums where scaffold is "all_scaffolds"
    all_mapping_info = all_mapping_info[all_mapping_info['scaffold'] == 'all_scaffolds']

    # concat gene table
    all_genes = concat_df(args.genes, args.taxa)

    # write output
    # if output directory does not exist, create it
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, exist_ok=True)
    all_genome_info.to_csv(os.path.join(args.outdir, 'all_genome_info.tsv'), sep='\t', index=False)
    all_mapping_info.to_csv(os.path.join(args.outdir, 'all_mapping_info.tsv'), sep='\t', index=False)
    all_genes.to_csv(os.path.join(args.outdir, 'all_genes_info.tsv'), sep='\t', index=False)

if __name__ == "__main__":
    main()