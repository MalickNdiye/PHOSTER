#!/usr/bin/env python

import sys
import argparse
import re
import traceback
from itertools import groupby
from Bio import SeqIO
import pandas as pd

__author__ = "Malick Ndiaye"
__title__ = "ParseCDHIT"



def get_args():
    """Parse command line arguments"""
    desc = (
        """Retrieve cluster sequences from the result of CD-HIT using the .clstr file"""
    )
    epi = """Returns a fasta file of the representative sequences and a tabular file of the cluster membership"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)
        parser.add_argument("-i",
            "--clusterfile", action="store",
             help='CD-HIT output file (ends in ".clstr").'
        )

        parser.add_argument("-o",
            "--clstout",
            action="store",
            help="output ",
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

def parse_cdhit(clstrfile):
    """
    create open dataframe to store cluster membership
    dataframe will have following columns:
    contig_id, cluster_id, cluster_size, representative, percent_identity
    cluster_id is the cluster number from the cdhit output
    cluster_size is the number of contigs in the cluster
    representative is a boolean indicating whether the contig is the representative sequence for the cluster
    percent_identity is the percent identity of the contig to the representative sequence
    """
    df = pd.DataFrame(columns=['contig_id', 'cluster_id', 'representative', 'percent_identity'])

    # open cdhit output file
    with open(clstrfile, "r") as f:
        # iterate over lines in cdhit output file
        for line in f:
            if line.startswith(">"):
                line = line.split()
                cluster_id = line[1].strip(">") # get cluster id
                print("\nprocessing cluster " + cluster_id + "...")
            else:
                cluster_id = cluster_id # use previous cluster id
                # split line on spaces
                line = line.split()
                # get contig id
                contig_id = line[2].strip(">").strip("...")
                # get representative boolean
                if "*" in line:
                    representative = True
                else:
                    representative = False
                # get percent identity
                if representative:
                    percent_identity = 100
                else:
                    percent_identity = float(line[4].strip("%").strip("+/").strip("-/"))
                

                print("\tcontig " + contig_id + " is " + str(percent_identity) + "% identical to the representative sequence")

                # add row to dataframe
                df = df.append({'contig_id': contig_id, 'cluster_id': cluster_id, 'reference': representative, 'percent_identity': percent_identity}, ignore_index=True)
        
        # Add cluster size column to dataframe
        df['cluster_size'] = df.groupby('cluster_id')['cluster_id'].transform('count')

        # add sample column, it is the first part of the contig id and put it as the first column
        df['sample'] = df['contig_id'].str.split("_").str[0]
        cols = df.columns.tolist()
        cols = cols[-1:] + cols[:-1]
        df = df[cols]

    return(df)
    

def main():
    """Call functions"""

    args = get_args()
   
    print("\n\n################################################")
    print("parse cdhit output")
    print("################################################")
    
    df_clst = parse_cdhit(args.clusterfile)

    # write cluster membership dataframe to file
    df_clst.to_csv(args.clstout, sep="\t", index=False)

if __name__ == "__main__":
    main()