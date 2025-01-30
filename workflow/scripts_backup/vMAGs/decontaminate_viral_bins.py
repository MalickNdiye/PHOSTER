#!/usr/bin/env python

import argparse
import sys
import pandas as pd
import os
from Bio import SeqIO

__author__ = "Malick Ndiaye"
__title__ = "Decontaminate viral bins"


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
            "--contam_data",
            action="store",
            help="contamination data"
        )
        
       # output files
        parser.add_argument("-o",
            "--outfasta", action="store",
             help='output decontaminated fasta file'
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

def decontam(fasta, data):
    # read fasta file as list of SeqIO objects
    fasta_list = list(SeqIO.parse(fasta, "fasta"))

    # read data file as dataframe
    data_df = pd.read_csv(data, sep="\t")

    for record in fasta_list:
        if record.id in data_df["contig"].tolist():
            print(record.id + " is contaminated. decomtaminating...")

            # trim the sequence according tot he start and end position in the data file
            start=list(data_df[data_df["contig"]==record.id]["start"])
            end=list(data_df[data_df["contig"]==record.id]["end"])
            new_name = data_df[data_df["contig"]==record.id]["new_name"].tolist()

            for i in range(len(start)):
                # add as many records as there are new names
                new_record = record[int(start[i])-1:int(end[i])]
                new_record.id = new_name[i] + "_prophage"
                new_record.description = ""
                fasta_list.append(new_record)
            
            # remove old record from the list, do not use remove() because it will remove only the first occurence
            fasta_list = [x for x in fasta_list if x.id != record.id]


            
        else:
            print(record.id + " is not contaminated. Keeping...")
            continue
    
    # remove "_prophage" from the record id
    for record in fasta_list:
        record.id = record.id.replace("_prophage","")

    return(fasta_list)

def write_fasta(fasta_list, outfasta):
    SeqIO.write(fasta_list, outfasta, "fasta")

def main():
    args = get_args()
    fasta_list = decontam(args.infasta, args.contam_data)
    write_fasta(fasta_list, args.outfasta)

if __name__ == "__main__":
    main()

        