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
__title__ = "PolishVMAGs"


def get_args():
    """Parse command line arguments"""
    desc = (
        """filter dereplicated viral contigs to remove bacteria sequences and create summary table of filtering"""
    )
    epi = """Bacterial contigs are defied as:
            1. Contigs with a bacterial marker genes as dfined by CheckM
            2. contigs longer than 500kb
          """

    parser = argparse.ArgumentParser(description=desc, epilog=epi)

    # Required Inputs
    parser.add_argument("-d","--rep", help="dRep output directory")
    parser.add_argument("-i", "--infasta",  help='input fasta file')
    parser.add_argument("-c",
        "--oldfasta", nargs='+',
        help="old fasta file"
    )

    # Required Outputs
    parser.add_argument("-o",
        "--fastout",
        help="output fasta file",
    )
    parser.add_argument("-t",
        "--outtab",
        help="output fasta file",
    )



    if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    return(parser.parse_args())
        


def main():
    """Call functions"""

    print("################################################")
    print("parsing dRep output")
    print("################################################")

    args = get_args()

    # read file data_tables/Cdb.csv in drp output directory using pandas
    # this file contains the dereplicated genomes and their associated clusters
    print("\nreading dRep output file " + args.rep + "/data_tables/Cdb.csv")
    clst_info=pd.read_csv(args.rep + "/data_tables/Cdb.csv", sep=",")
    clst_list=clst_info["genome"].tolist()
    clst_list=[clst.replace(".fasta", "") for clst in clst_list]

  
    # read the input fasta file using SeqIO as a list of SeqRecords
    print("\nreading input fasta file " + args.infasta)
    fasta_list=list(SeqIO.parse(args.infasta, "fasta"))

    # read the old fasta file list and concatenate it as a list of SeqRecords
    print("\nreading old fasta files " + str(args.oldfasta))
    old_fasta_list=[]
    for fasta in args.oldfasta:
        old_fasta_list+=list(SeqIO.parse(fasta, "fasta"))
    
    old_fasta_headers=[fasta.id for fasta in old_fasta_list]
   

    # find element in clst_list that are not in old_fasta_list headers
    # these are the contigs that were removed by erroneously by dRep
    # and need to be added back to the fasta file
    missing_contigs=[clst for clst in old_fasta_headers  if clst not in clst_list ]

    # add the missing contigs to the fasta list
    print("\nadding " + str(len(missing_contigs)) + " missing contigs to the fasta file")
    for contig in missing_contigs:
        fasta_list.append(old_fasta_list[old_fasta_headers.index(contig)])

    # if an header in fasta_list contains the strin "||", replace it and everything that comes after with "_" plus the sum of all numbers in the header
    # this is to avoid problems with the headers in future
    print("\nreplacing || in fasta headers")
    # create empty dataframe to store the name changes
    name_changes=pd.DataFrame(columns=["name_pre", "name_post"])

    for fasta in fasta_list:
        if "||" in fasta.description:
            print("\tfound || in header " + fasta.description)

            new_header=re.sub(r"\|\|.*", "", fasta.description) + "_" + str(sum([int(s) for s in re.findall(r'\d+', fasta.description)]))
            name_changes=name_changes.append({"name_pre":fasta.description, "name_post":new_header}, ignore_index=True)

            print("\t--->replacing with " + new_header)
            fasta.description=""
            fasta.id=""
            fasta.id=new_header
            fasta.description=new_header
        
        # if the header is not in the old fasta file, add it to the name changes dataframe
        else:
            name_changes=name_changes.append({"name_pre":fasta.description, "name_post":fasta.description}, ignore_index=True)

            
    # write the fasta list to the output fasta file
    print("writing output fasta file " + args.fastout)
    SeqIO.write(fasta_list, args.fastout, "fasta")

    # write the name changes to the output table
    print("writing output table " + args.outtab)
    name_changes.to_csv(args.outtab, sep="\t", index=False)

if __name__ == "__main__":
    main()



    
    



    




