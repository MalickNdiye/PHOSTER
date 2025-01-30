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
__title__ = "Filter Bacteria from viral Contigs"


def get_args():
    """Parse command line arguments"""
    desc = (
        """This scripts remove from a fasta file all contigs that are not placd at the root by checkm"""
    )
    epi = """Returns a fasta file of the Binned Contigs sequences and a tabular file of the cluster membership"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # input files
        parser.add_argument("-f",
            "--infasta", action="store",
             help='input fasta file'
        )
        parser.add_argument("-c",
            "--checkm",
            action="store",
            help="checkm file"
        )
        parser.add_argument("-m",
            "--mtdata",
            action="store",
            help="metadatafile"
        )
        
       # output files
        parser.add_argument("-o",
            "--outfasta", action="store",
             help='output fasta file, filtered for viral contigs'
        )
        parser.add_argument("-t",
            "--outtab", action="store",
             help='updated metadata file'
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

def get_bacterial_contigs(checkm_file):
    """Parse checkm file to get bacterial contigs"""

    # get sample name
    sample_name = os.path.basename(checkm_file).split('_')[0]
    
    print('\nGetting bacterial contigs for sample {}'.format(sample_name))

    df = pd.read_csv(checkm_file, sep='\t')
    
    # get elements of column "Bin Id" that have "Marker lineage" != "root (UID1)" and Genome size (bp)/1000000 > Completeness
    bacterial_contigs = df.loc[df['Genome size (bp)']/1000000 < df['Completeness'], 'Bin Id'].tolist()

    print('-->Found {} bacterial contigs'.format(len(bacterial_contigs)))

    return(bacterial_contigs)

def update_mtadata_file(mtdata, bacterial_contigs):
    """Update metadata file with bacterial contigs"""
    print('\nUpdating metadata file')

    # read metadata file
    df = pd.read_csv(mtdata, sep='\t')

    # if and element of the column "genome_id" is in bacterial_contigs, then update the column "retained" to "no" and the column "reason" to "bacterial"
    df.loc[df['genome_id'].isin(bacterial_contigs), 'retained'] = 'no'
    df.loc[df['genome_id'].isin(bacterial_contigs), 'reason'] = 'bacterial'

    return(df)

def filter_fasta_file(infasta, bacterial_contigs):
    """Filter fasta file for bacterial contigs"""

    print('\nFiltering fasta file for bacterial contigs')

    # read fasta file
    records = list(SeqIO.parse(infasta, 'fasta'))

    # filter records
    filtered_records = [record for record in records if record.id not in bacterial_contigs]

    return(filtered_records)

def main():
    """Main function"""

    # get arguments
    args = get_args()

    print('\n#######################################################################')
    print('Filter Bacteria from viral Contigs')
    print("input fasta file: {}".format(args.infasta))
    print("checkm_file: {}".format(args.checkm))
    print("detadata file {} \n will be updated as".format(args.mtdata, args.outtab))
    print("output fasta file: {}".format(args.outfasta))
    print('#######################################################################\n')

    # get bacterial contigs
    bacterial_contigs = get_bacterial_contigs(args.checkm)

    # update metadata file
    df = update_mtadata_file(args.mtdata, bacterial_contigs)

    # write updated metadata file
    df.to_csv(args.outtab, sep='\t', index=False)

    # filter fasta file, if path to output fasta file does not exist, create it
    print('\nWriting filtered fasta file to {}'.format(args.outfasta))
    if not os.path.exists(os.path.dirname(args.outfasta)):
        os.makedirs(os.path.dirname(args.outfasta))
    filtered_records = filter_fasta_file(args.infasta, bacterial_contigs)

    # write filtered fasta file, if path to output fasta file does not exist, create it
    print('\nWriting filtered fasta file to {}'.format(args.outfasta))
    if not os.path.exists(os.path.dirname(args.outfasta)):
        os.makedirs(os.path.dirname(args.outfasta))
    SeqIO.write(filtered_records, args.outfasta, 'fasta')

    print('\nDONE!')

if __name__ == "__main__":
    main()
