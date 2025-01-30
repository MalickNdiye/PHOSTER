import sys
import argparse
import re
import traceback
from itertools import groupby
import pandas as pd
from Bio import SeqIO
import glob

__author__ = "Malick Ndiaye"
__version__ = "1.0"
__title__ = "concat viral contigs that are shorter than 500kb"


def get_args():
    """Parse command line arguments"""
    desc = (
        """concat viral contigs that are shorter than 500kb"""
    )
    epi = """concat viral contigs that are shorter than 500kb
          """

    parser = argparse.ArgumentParser(description=desc, epilog=epi)

    # Required Inputs
    parser.add_argument("-i", "--infasta",  help='input fasta files', nargs='+')

    # Required Outputs
    parser.add_argument("-o",
        "--fastout",
        help="output fasta file",
    )

    if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    return(parser.parse_args())

def concatenate_fasta(infastas, fastout):
    """concatenate fasta files"""

    with open(fastout, "w") as out:
            for i in infastas:
                with open(i, "r") as f:
                    for record in SeqIO.parse(f, "fasta"):
                        if len(record.seq) < 500000:
                            SeqIO.write(record, out, "fasta")
                        else:
                            print("contig " + record.id + " longer than 500kb. skipping...")

        


def main():
    """Call functions"""

    print("################################################")
    print("concatenate viral contigs shorter than 500kb")
    print("################################################")

    args = get_args()

    print("\ninput fasta files: " + "\n" "\n\t".join(args.infasta))
    print("\noutput fasta file: " + str(args.fastout))

    concatenate_fasta(args.infasta, args.fastout)
                            

if __name__ == "__main__":
    main()