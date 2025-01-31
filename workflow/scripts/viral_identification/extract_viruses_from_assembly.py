#!/usr/bin/env python

import sys
import os
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="""Given a table of identified viral contigs and a FATSA file, this
    script extract the viral contigs from the FASTA and creates a new fatsa file of only viral contigs.
    if a prophage is identified, the script also cut the sequence to correspond only to the prophage""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("-i", "--infasta", help="Paths to input FASTA file")
parser.add_argument("-t", "--intab", help="Paths to input viral identification table")
parser.add_argument("-o", "--outfasta", help="Path to output FASTA file")

if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
    parser.print_help()
    sys.exit()

args = parser.parse_args()

print("**** EXTRACTING VIRAL CONTIGS FROM "+   args.infasta + " ****\n")

# get sample name
sam=args.infasta.split("/")[-1].split("_")[0]

# open inputs
tab=pd.read_csv(args.intab, delimiter="\t", low_memory=False).query('sample==@sam')
contigs=tab["contig"].tolist() # list of contigs that are viruses

outfatsa=open(args.outfasta, "w")
for fasta in SeqIO.parse(args.infasta, "fasta"):
    contig=fasta.description
    if contig in contigs:
        print("\n---contig " + contig + " is a virus")
        sequence=fasta.seq

        if "yes" in tab.query('contig==@contig')["prophage"].tolist():
            print("\tit is in fact a prophage! Trimming the sequence.")

            start=tab.query('contig==@contig')["start"].dropna().tolist()
            end=tab.query('contig==@contig')["stop"].dropna().tolist()
            new_name=tab.query('contig==@contig')["new_contig"].dropna().tolist()

            # remove white spaces from the list
            new_name = [x.strip(' ') for x in new_name]

            for i in range(len(start)):
                print(int(start[i])-1)
                print(int(end[i]))
                outfatsa.write(">"+str(new_name[i])+"\n"+str(sequence[int(start[i])-1:int(end[i])])+"\n")

        else:
            new_name=tab.query('contig==@contig')["new_contig"].dropna().tolist()
            new_name = [x.strip(' ') for x in new_name]
            outfatsa.write(">"+str(new_name[0])+"\n"+str(sequence)+"\n")

    else:
        print("\n---contig " + contig + " is NOT a virus")
        continue

print("\nDONE!")
