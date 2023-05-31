import sys
import os
import argparse
from Bio import SeqIO
import pathlib
import pandas as pd
from os import system

parser = argparse.ArgumentParser(
    description="""This script takes the dereplicated Bins generated after vRhyme dereplication and filter away those than have less than 3kb length.
    plus it produces a table smmarizing the dereplication step.""", #TODO in future will also filter based on viral prediction tools
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("infile", help="Paths to input FASTA file")
parser.add_argument("intab", help="Path table composite sequence table")
parser.add_argument("outfile", help="Path to output FASTA file")
parser.add_argument("outtab", help="Path to output table filter")
parser.add_argument("outderep", help="Path to output derep table")
parser.add_argument("outsing", help="Path to output single genomes")
args = parser.parse_args()

# Parse dereplication tables
# open composite contigs as dictionary
#args.intab
composite_dict={}
with open(args.intab, 'r') as file:
    first_line = file.readline()
    for line in file:
        line_list=line.rstrip("\n").split("\t")
        composite_dict[line_list[0]]=line_list[1:]
        print("\n---sequence " + line_list[0] + " is composed of:\n" + "\t"+  "\n\t".join(composite_dict[line_list[0]]))

df=pd.DataFrame(dict([ (k,pd.Series(v)) for k,v in composite_dict.items() ])). melt(var_name="composite_contig", value_name="contig").dropna()

df.to_csv(args.outderep, index=False, sep="\t")

#loop through contigs and filter them
outtab=open(args.outtab, "w")
outtab.write("contig\tlength\tretained\n")

outfasta=open(args.outfile, "w")
pathlib.Path(args.outsing).mkdir(parents=True, exist_ok=True)

for record in SeqIO.parse(args.infile, "fasta"):
    length=len(record.seq)
    contig=(record.id)
    entry=[str(contig), str(length)]
    if length > 3000:
        print("\n--- Seqeunce: " + str(record.id) + " passes filter with length=" + str(length))

        # write filtered genome in single fasta file
        sing_fasta_p=os.path.join(args.outsing, contig+".fasta")
        sing_fasta=open(sing_fasta_p, "w")
        sing_fasta.write(">"+str(record.id)+"\n"+str(record.seq)+"\n")

        # write filtered genome filtering table
        entry.append("yes")
        outtab.write("\t".join(entry)+"\n")

        # write filtered genome in conct file
        outfasta.write(">"+str(record.id)+"\n"+str(record.seq)+"\n")
        
    else:
        print("\n--- Seqeunce: " + str(record.id) + " DOES NOT pass filter with length=" + str(length))
        entry.append("no")
        outtab.write("\t".join(entry)+"\n")
