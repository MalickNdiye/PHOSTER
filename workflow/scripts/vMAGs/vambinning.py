import sys
import os
import argparse
from Bio import SeqIO
import pathlib
import pandas as pd
from os import system

parser = argparse.ArgumentParser(
    description="""This script takes the VAMB output as parsed and creates individual files for each bin as well as a concat file with all bins. 
    The bins take the name of the bin or of the ontig i they were not binned by VAMB""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("intab", help="Path table containing clustering info")
parser.add_argument("infile", help="Paths to input FASTA file")
parser.add_argument("outpath", help="Path to output directory")
args = parser.parse_args()


vamb_clust=pd.read_csv(args.intab, delimiter="\t")
vamb_c=list(vamb_clust["vamb_name"])


#loop through contigs and assign them to correct bin
seqs={}
for record in SeqIO.parse(args.infile, "fasta"):
    contig=str(record.id)
    sequence=str(record.seq)
    if contig in vamb_c:
        vae=vamb_clust.query('vamb_name == @contig')["binname"].tolist()[0]
        print(vae)

        if vae in seqs.keys():
            seqs[vae]+=sequence
        else:
            seqs[vae]=sequence
    else:
        seqs[contig]=sequence


# write each bin concat sequence in a single fasta file
outsingle=os.path.join(args.outpath, "single_bins/")
pathlib.Path(outsingle).mkdir(parents=True, exist_ok=True)
for id in seqs.keys():
    filename=os.path.join(outsingle,id)+".fasta"
    file=open(filename, "w")
    print("creating bin: " + filename)
    file.write(">"+id+"\n"+seqs[id]+"\n")
    file.close()

out_concat="viral_bins.fasta"
