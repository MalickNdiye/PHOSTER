import argparse
import os
from Bio import SeqIO
import pandas as pd

def get_args():
    parser = argparse.ArgumentParser(description='Parse Crass output')
    parser.add_argument('-i', '--input', help='Input fasta', required=True)
    parser.add_argument('-o', '--output', help='output table', required=True)
    args = parser.parse_args()
    return args


def main():
    args = get_args()

    df = pd.DataFrame(columns=['contig', 'length'])
    for record in SeqIO.parse(args.input, 'fasta'):
        df = df.append({'contig': record.description, 'length': len(record.seq)}, ignore_index=True)
    

if __name__ == '__main__':
    main()