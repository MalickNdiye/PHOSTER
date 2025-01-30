import pandas as pd
import argparse
import os
from Bio import SeqIO


def get_args():
    parser = argparse.ArgumentParser(description='Parse Crass output')
    parser.add_argument('-d', '--indr', help='Input fasta DR', required=True)
    parser.add_argument('-s', '--insp', help='Input fasta spacers', required=True)
    parser.add_argument('-t', '--intab', help='input crass stat table', required=True)
    parser.add_argument('-o', '--out', help='output table', required=True)
    parser.add_argument('-f', '--outspacer', help='output spacer fasta', required=True)
    parser.add_argument('-r', '--outdr', help='output dr fasta', required=True)
    args = parser.parse_args()
    return args

def reformat_fasta(fasta, outfile):
    """
    This fuction takes a fasta files and add the sample name at the beginning of the header
    (the sample name in the charachter before the first underscore in the filename)
    """
    # check if the output directory exists, if not create it
    if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))
                
    print("reformatting fasta...")
    sample = os.path.basename(fasta).split('_')[0]
    with open(outfile, 'w') as out:
        for record in SeqIO.parse(fasta, 'fasta'):
            record.id = sample + '_' + record.id
            record.description = sample + '_' + record.description

            SeqIO.write(record, out, 'fasta')

def reformat_table(tab, outfile):
    """
    this function adds a column with the sample name to the crass stat table. also, it remove the last row of the table (contains sum of all samples)
    """

    if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))

    print("reformatting crass table...")
    sample = os.path.basename(tab).split('_')[0]
    df = pd.read_csv(tab, sep='\t')
    df['sample'] = sample
    cols = df.columns.tolist()
    cols = cols[-1:] + cols[:-1]
    df = df[cols]
    df = df[:-1]

    if not os.path.exists(os.path.dirname(outfile)):
                os.makedirs(os.path.dirname(outfile))
    df.to_csv(outfile, sep='\t', index=False)

def main():
    args = get_args()
    reformat_fasta(args.indr, args.outdr)
    reformat_fasta(args.insp, args.outspacer)
    reformat_table(args.intab, args.out)

if __name__ == '__main__':
    main()