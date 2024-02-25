#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import os
from Bio import SeqIO
import re


__author__ = "Malick Ndiaye"
__title__ = "filter assemblies CheckV"


def get_args():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser()

        # input files
        parser.add_argument("-i",
            "--infasta", action="store",
             help='input fasta file'
        )
        parser.add_argument("-d",
            "--checkv_data",
            action="store",
            help="binning data"
        )
        
       # output files
        parser.add_argument("-o",
            "--outfasta", action="store",
             help='output concatenated vMAGs fasta file'
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

def get_checkv_data(checkv_data):
    """
    This function reads the checkv data and returns a list of retained contig ids
    """
    checkv_df = pd.read_csv(checkv_data, sep="\t")

    # get all rows where column retained==1
    checkv_df_filt = checkv_df[checkv_df["retained"]==1]

    # get column contig_id to list
    contig_id = checkv_df_filt["contig_id"].tolist()

    return(contig_id)


def filter_assemblies(infasta, checkv_data, outfasta):
    """
    This function filters the input fasta file based on the checkv data
    """
    # get contig_id list
    contig_id = get_checkv_data(checkv_data)

    # filter fasta file
    SeqIO.write((seq for seq in SeqIO.parse(infasta, "fasta") if seq.id in contig_id), outfasta, "fasta")

def main():
    """Main function"""

    # get arguments
    args = get_args()

    # run function
    filter_assemblies(args.infasta, args.checkv_data, args.outfasta)

if __name__ == "__main__":
    main()