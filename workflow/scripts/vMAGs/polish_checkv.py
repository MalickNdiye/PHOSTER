#!/usr/bin/env python

import sys
import argparse
import re
import traceback
from itertools import groupby
from Bio import SeqIO
import pandas as pd

__author__ = "Malick Ndiaye"
__title__ = "Polish checkV"



def get_args():
    """Parse command line arguments"""
    desc = (
        """integrate data from checV pre and post binning, and rename contigs from prophages"""
    )
    epi = """Returns a fasta file of the representative sequences and a tabular file of the cluster membership"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)
        parser.add_argument("-f",
            "--infasta", action="store",
             help='input fasta file'
        )
        parser.add_argument("-b",
            "--pre_binning",
            action="store",
            help="pre binning checkV file",
        )
        parser.add_argument("-a",
            "--post_binning",
            action="store",
            help="post binning checkV file",
        )
        parser.add_argument("-m",
            "--metadata",
            action="store",
            help="contigs metadata file",
        )


        parser.add_argument("-t",
            "--checkv_out",
            action="store",
            help="checkV output file",
        )
        parser.add_argument("-o",
            "--outfasta",
            action="store",
            help="output fasta file",
        )

       

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occurred with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()

def main():
    """call functions"""

    args=get_args()

    # read fasta file using SeqIO as a list of SeqRecords
    fasta = list(SeqIO.parse(args.infasta, "fasta"))

    # read metadata file
    metadata = pd.read_csv(args.metadata, sep='\t', header=0)

    # read pre binning checkV file
    pre_binning = pd.read_csv(args.pre_binning, sep='\t', header=0)
    # add "contig" column of metadata to pre_binning dataframe in fuction of "old_name" column
    pre_binning = pd.merge(pre_binning, metadata[['old_name', 'contig', "bin"]], how='left', left_on='contig_id', right_on='old_name')    

    # read post binning checkV file
    post_binning = pd.read_csv(args.post_binning, sep='\t', header=0)

    



if __name__ == "__main__":
    main()
