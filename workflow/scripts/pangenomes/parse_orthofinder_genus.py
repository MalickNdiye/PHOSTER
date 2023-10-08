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
__title__ = "Parse Orthofinder output"


def get_args():
    """Parse command line arguments"""
    desc = (
        """
        this script parses the orthofinder output to get various summary tables
        """
    )
    epi = """"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # input files
        parser.add_argument("-i",
            "--infile", action="store",
             help='parsed orthofinder file'
        )
        parser.add_argument("-g",
            "--genus", action="store",
             help='genus of interest'
        )

       # output files
        parser.add_argument("-o",
            "--outdir", action="store",
             help='output directory'
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



def split_group(file, outdir, genus):
    """
    this function splits the dataframe by genera and saves it in a directory named after the genus
    """

    # read the dataframe
    df = pd.read_csv(file, sep="\t", low_memory=False)

    # create a dictionary of dataframes with key the genus and value the dataframe
    df_dict = df[df["genus"] == genus]
    # create a directory for each genus
   
    print("Processing genus: ", genus, "\n")



        # create the directory
    genus_dir = os.path.join(outdir)
    if not os.path.exists(genus_dir):
        os.makedirs(genus_dir, exist_ok=True)
        # get the dataframe
    df_genus = df_dict
    # write the dataframe to file
    outfile = os.path.join(genus_dir, genus+"_Orthogroups_table.tsv")

    with open(outfile, "w") as f:
        df_genus.to_csv(f, sep="\t", index=False)

    return(df_dict)

def get_single_copy_ogs(df, outdir, genus):
    """"This function takes a  dictionary of dataframes and returns a list of dataframes with only single copy OGs"""


    # get the dataframe
    df = df[df["OG_size_in_genome"] == 1]
        

    genus_dir = os.path.join(outdir)

    # write the dataframe to file
    outfile = os.path.join(genus_dir, genus+"_single_copy_Orthogroups_table.tsv")

    with open(outfile, "w") as f:
        df.to_csv(f, sep="\t", index=False)

def create_og_file(df, outdir, genus):
    """this function creates a file with the following columns: Genome, OG. it does this for each genus"""



    print("create mOTUpan file for genus: ", genus, "\n")

    # get the dataframe
    df_genus = df[["genome", "OG"]]

    # transform the dataframe to a dictionary with the genome as key and a list of OGs as value
    df_dict = df_genus.groupby("genome")["OG"].apply(list).to_dict()

    genus_dir = os.path.join(outdir)

    # write the dataframe to file
    # ech row is a genome and the columns are the OGs separated by a tab
    outfile = os.path.join(genus_dir, genus+"_mOTUpan_table.tsv")

    with open(outfile, "w") as f:
        for k, v in df_dict.items():
            f.write(k + "\t" + "\t".join(v) + "\n")








def main():
    args=get_args()


    # split the dataframe by genera
    df_dict = split_group(args.infile, args.outdir, args.genus)

    # get single copy OGs
    get_single_copy_ogs(df_dict, args.outdir, args.genus)

    # create mOTUpan table
    create_og_file(df_dict, args.outdir, args.genus)



if __name__ == "__main__":
    main()



        

