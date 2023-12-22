#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import os
from Bio import SeqIO


__author__ = "Malick Ndiaye"
__title__ = "Bin vRhyme MAGs"


def get_args():
    """Parse command line arguments"""
    desc = (
        """This script bins vRhyme contigs into vMAGs based on the cluster membership of the representative sequences.
        """
    )
    epi = """Returns a fasta file of the Binned Contigs sequences and a tabular file of the cluster membership"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # input files
        parser.add_argument("-i",
            "--infasta", action="store",
             help='input fasta file'
        )
        parser.add_argument("-d",
            "--bin_data",
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


def bin_contigs (infasta, bin_data, outfasta):
    """
    This function bin the contigs in the input fasta file based on the bin membership of the representative sequences
    """

    # open input fasta as a dictionary
    print("opening fasta file...")
    fasta_dict = SeqIO.to_dict(SeqIO.parse(infasta, "fasta"))

    # open binning data
    print("opening binning data...")
    binning_data = pd.read_csv(bin_data, sep='\t')

    outfasta_dict = fasta_dict.copy()


    print("binning contigs...")
    # remove all keys and values that are in binning data in the column "contig"
    for key in fasta_dict.keys():
        if key in binning_data["contig"].tolist():
            del outfasta_dict[key]
        else:
            outfasta_dict[key]=outfasta_dict[key].seq
    
    # add keys that are named as the column "bin_name" in binning data
    for bin in binning_data["bin_name"].tolist():
        outfasta_dict[bin] = ""

    # concatenate the sequences of the contigs in the same bin
    for key in outfasta_dict.keys():
        if key in binning_data["bin_name"].tolist():
            contigs = binning_data.query('bin_name==@key')["contig"].tolist()
            for contig in contigs:
                outfasta_dict[key] += fasta_dict[contig].seq

    # write the output fasta file
    print("writing output fasta file...")
    with open(outfasta, 'w') as out:
        for key in outfasta_dict.keys():
            out.write(">"+key+"\n"+str(outfasta_dict[key])+"\n")

def main():
    args = get_args()
    bin_contigs(args.infasta, args.bin_data, args.outfasta)

if __name__ == '__main__':
    main() 