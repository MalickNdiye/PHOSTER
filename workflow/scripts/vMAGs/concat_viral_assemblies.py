import sys
import os
import argparse
from Bio import SeqIO

# This script parses the viral asemblies for VAMB
# It creaes a concat file with contigs renamed in the VAMB fashion and a metadata table of all contigs and wether they are retained or not

parser = argparse.ArgumentParser(
    description="""Creates the input FASTA file for Vamb.
Input should be one or more FASTA files, each from a sample-specific assembly.
If keepnames is False, resulting FASTA can be binsplit with separator 'C'.""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("outpath", help="Path to output FASTA file")
parser.add_argument("outfile", help="Path to output contig metadata table")
parser.add_argument("inpaths", help="Paths to input FASTA file(s)", nargs="+")
args = parser.parse_args()


concat= open(args.outpath, "w+")

tab= open(args.outfile, "w+")
header="\t".join(["old_name", "new_name", "length", "retained\n"])
tab.write(header)

for file in args.inpaths:
    print("\nParsing " + file)

    for entry in SeqIO.parse(file, "fasta"):
        length=len(entry.seq)
        print("-------------sequence: " + str(entry.id) + " Length="+str(length))
        
        if length>=2000:
            old=entry.id
            if "SAMPLE" in old:
                sam=old.split("_")[3]
            else:
                sam=old.split("_")[0]
            new=sam + "C" + old

            concat.write(">" + new + "\n" + str(entry.seq) +"\n")

            tab_new="\t".join([old, new, str(length), "yes\n"])
            tab.write(tab_new)
        else:
            old=entry.id
            tab_new="\t".join([old, "NaN", str(length), "no\n"])
            tab.write(tab_new)

tab.close()
concat.close()