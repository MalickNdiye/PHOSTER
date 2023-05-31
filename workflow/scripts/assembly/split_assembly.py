import os
from Bio import SeqIO


def split(fastafile, outfastadir   =   "splitoutput"):

    #"""Extract multiple sequence fasta file and write each sequence in separate file"""
    os.makedirs(outfastadir, exist_ok=True)

    sample_name=fastafile.split("/")[-1].split("_")[0]

    with open (fastafile) as FH:
        record = SeqIO.parse(FH, "fasta")
        file_count =  0

        for seq_rec in record:
            name= str(seq_rec.description)
            file_count  =   file_count  +   1
            with open("%s/%s.fasta" % (outfastadir,str(name)), "w") as FHO:
                SeqIO.write(seq_rec, FHO, "fasta")

    if file_count== 0:
        raise Exception("No valid sequence in fasta file")

    return("DONE")

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description="Extract multiple sequence fasta file and write each sequence in separate file")

    parser.add_argument('-f','--fastafile',
                        action  ="store",
                        help="Fasta File for parsing")
    parser.add_argument('-d','--outfastadir',
                        action  ="store",
                        help    ="Fasta File output directory")

    args = parser.parse_args()
    split(fastafile = args.fastafile,
          outfastadir= args.outfastadir)