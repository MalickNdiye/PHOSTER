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
            "--indir", action="store",
             help='input directory containing orthofinder output files'
        )
        parser.add_argument("-c",
            "--clust", action="store",
             help='MAGs clustering file'
        )

       # output files
        parser.add_argument("-o",
            "--outfile", action="store",
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

def create_og_table(file, clust_file):
    """
    this function parses the "Orthogroups.tsv" file to return a file with the following columns: Genome, id, genus, gene, OG, OG_size_in_genome
    """

    # open clust file
    clust = pd.read_csv(clust_file, sep="\t")
    # get dictionary of with key "Bin Id" and value "Genus"
    clust_dict_genus = dict(zip(clust["Bin Id"], clust["genus"]))
    # get dictionary of with key "Bin Id" and value "id"
    clust_dict_id = dict(zip(clust["Bin Id"], clust["id"]))

    # open orthogroups file
    df = pd.read_csv(file, sep="\t", low_memory=False)

    # the structure of the orthogroups file is as follows:
    # the first column is the orthogroup id
    # the other columns are the genomes
    # the cells contain the genes of the genomes that are in the orthogroup as a , separated list

    # pivot the dataframe to have the genomes and OGs as rows
    df = df.melt(id_vars=["Orthogroup"], var_name="Genome", value_name="Genes")
    # split the genes column into a list
    df["Genes"] = df["Genes"].str.split(",")
    # explode the list to have one gene per row
    df = df.explode("Genes")
    # remove the NaN values
    df = df.dropna()

    # get the size of each OG in each genome
    df["OG_size_in_genome"] = df.groupby(["Orthogroup", "Genome"])["Genes"].transform("size")

    # remove "_gene" from all genomes
    df["Genome"] = df["Genome"].str.replace("_genes", "")

    # get the genus of each genome
    df["genus"] = df["Genome"].map(clust_dict_genus)
    # get the id of each genome
    df["id"] = df["Genome"].map(clust_dict_id)

    # rename the columns
    df = df.rename(columns={"Orthogroup": "OG", "Genes": "gene", "Genome": "genome"})
    # reorder the columns
    df = df[["genome", "id", "genus", "gene", "OG", "OG_size_in_genome"]]

    
        
    return(df)




def main():
    args=get_args()

    orthofile=os.path.join(args.indir, "Orthogroups.tsv")
    clustfile=args.clust

    # create dataframe
    df = create_og_table(orthofile, clustfile)

    # write dataframe to file, it is very big so write it efficiently
    outfile=args.outfile
    # check if the output directory exists
    if not os.path.exists(os.path.dirname(outfile)):
        os.makedirs(os.path.dirname(outfile), exist_ok=True)
    # write the dataframe to file
    with open(outfile, "w") as f:
        df.to_csv(f, sep="\t", index=False)



if __name__ == "__main__":
    main()



        

