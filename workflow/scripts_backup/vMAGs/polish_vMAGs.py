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

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # Required Inputs
        parser.add_argument("-i",
            "--infasta", action="store", 
            help='input fasta file'
        )
        parser.add_argument("-f",
            "--intab",
            action="store",
            help="input post-vRhyme filtering table"
        )
        parser.add_argument("-m",
            "--checkm",
            action="store",
            help="checkM output table"
        )
        parser.add_argument("-v",
            "--checkv",
            action="store",
            help="checkV output table"
        )
        parser.add_argument("-d",
            "--checkv_old",
            action="store",
            help="old checkV directory"
        )
        parser.add_argument("-n",
            "--names_drep",
            action="store",
            help="tsv file with old and new names of vMAGs after dRep dereplication"
        )


        # Required Outputs
        parser.add_argument("-o",
            "--fastout",
            action="store",
            help="output fasta file",
        )
        parser.add_argument("-t",
            "--tabout",
            action="store",
            help="output summary table",
        )

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occurred with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()


def main():
    """Call functions"""

    args = get_args()

    print("Starting to filter viral contigs...")

    # read names_drep table
    print("\tReading names_drep table: " + args.names_drep + "...")
    names_drep = pd.read_csv(args.names_drep, sep='\t', index_col=0)

    # read filtering table
    print("\tReading filtering table: " + args.intab + "...")
    filter_table = pd.read_csv(args.intab, sep='\t', index_col=0)
    filter_table = filter_table.drop(columns=['filter', 'length'])
    filter_table = filter_table.merge(names_drep, how='left', left_on='bin', right_on='name_pre')
    filter_table = filter_table.rename(columns={'name_post': 'names_drep'})

    # read checkM table and get elements in column "Bin Id" to list for rows where "Marker lineage" is not "root (UID1)" (i.e. most likely bacteria)
    print("\tReading checkM table: " + args.checkm + "...")
    checkm_table = pd.read_csv(args.checkm, sep='\t', index_col=0)
    bacterial_bins = checkm_table.loc[checkm_table['Marker lineage'] != 'root (UID1)', 'Bin Id'].tolist()
   

    # read checkV table
    print("\tReading checkV table: " + args.checkv + "...")
    checkv_table = pd.read_csv(args.checkv, sep='\t', index_col=0)

    # open and concatenate all the tables in the old checkV directory
    print("\tReading old checkV directory: " + args.checkv_old + "...")
    checkv_old_table = pd.concat([pd.read_csv(f, sep='\t', index_col=0) for f in glob.glob(args.checkv_old + "/*.tsv")], ignore_index=True)
    checkv_old_table['bin'] = checkv_old_table.apply(lambda row: filter_table.loc[filter_table['old_name'] == row['contig_id'], 'bin'].values[0], axis=1)

    # remove all rows where element "bin" is presnt multiple times
    checkv_old_table = checkv_old_table.drop_duplicates(subset=['bin'], keep=False)

    # remove "contig_id" column and rename "bin" column to "contig_id"
    checkv_old_table = checkv_old_table.drop(columns=['contig_id'])
    checkv_old_table = checkv_old_table.rename(columns={'bin': 'contig_id'})

    # add to checkv_old_table all rows in checkv_table where element "contig_id" is not present in checkv_old_table
    checkv_final_table = checkv_old_table.append(checkv_table[~checkv_table['contig_id'].isin(checkv_old_table['contig_id'])], ignore_index=True)


    # Read in the fasta file as list of SeqRecords
    print("\tReading fasta file: " + args.infasta + "...")
    viral_contigs = list(SeqIO.parse(args.infasta, "fasta"))

    print("\n############################################################")
    print("filtering viral contigs")
    print("############################################################")

    print("\nRemoving contig of bacterial origin")
    # remove contigs longer than 500kb and contigs in bacterial bins
    print("\tRemoving contigs longer than 500kb...")
    viral_contigs_nolong = [contig for contig in viral_contigs if len(contig.seq) < 500000]
    
    print("\tRemoving contigs in bacterial bins...")
    viral_contigs_nolong_nobact = [contig for contig in viral_contigs_nolong if contig.description not in bacterial_bins]

    # store names of removed contigs in list
    removed_contigs_bact = [contig.description for contig in viral_contigs if contig not in viral_contigs_nolong]

    # remvoe rows in intab if element of "bin" column is not in viral_contigs_nolong_nobact
    print("\nprocessing filtering table...")
    filter_table = filter_table[filter_table['names_drep'].isin([contig.description for contig in viral_contigs_nolong_nobact])]
    filter_table = pd.merge(filter_table, checkv_final_table, on=['names_drep', 'contig_id'], how='left')

    # rename names_drep column to bin
    filter_table = filter_table.rename(columns={'names_drep': 'drep_bin'})

    # if header in viral_contigs_nolong_nobact contains a whitespace, remove it and everything after it and eplace it with "_" + sum of all numbers in the old header
    # then change the "bin" column in filter_table to the new header
    for contig in viral_contigs_nolong_nobact:
        if " " in contig.description:
            new_header = contig.description.split(" ")[0] + "_" + str(sum([int(s) for s in re.findall(r'\d+', contig.description)]))
            filter_table.loc[filter_table['new_bin'] == contig.description, 'new_bin'] = new_header

            contig.description=""
            contig.id = ""
            contig.description = new_header
            contig.id = new_header

    n_low=filter_table.loc[filter_table['checkv_quality'] == 'Low-quality', 'checkv_quality'].count()
    n_high=filter_table.loc[filter_table['checkv_quality'] == 'High-quality', 'checkv_quality'].count()
    n_complete=filter_table.loc[filter_table['checkv_quality'] == 'Complete', 'checkv_quality'].count()
    n_medium=filter_table.loc[filter_table['checkv_quality'] == 'Medium-quality', 'checkv_quality'].count()
    n_na=filter_table.loc[filter_table['checkv_quality'] == 'Not-determined', 'checkv_quality'].count()

    print("\n################# SUMMARY #################")
    print("\nnumber of retained viral contigs: " + str(len(viral_contigs_nolong_nobact)))
    print("\tnot-determined checkV score: " + str(n_na) + " Genomes")
    print("\tlow-quality checkV score: " + str(n_low) + " Genomes")
    print("\tmedium-quality checkV score: " + str(n_medium) + " Genomes")
    print("\thigh-quality checkV score: " + str(n_high) + " Genomes")
    print("\tcomplete checkV score: " + str(n_complete) + " Genomes")
    print("\n###########################################")

    # write fasta file
    print("\nWriting filtered fasta file to: " + args.fastout)
    SeqIO.write(viral_contigs_nolong_nobact, args.fastout, "fasta")

    # write summary table
    print("\nWriting summary table to: " + args.tabout)
    filter_table.to_csv(args.tabout, sep='\t', index=False)


if __name__ == "__main__":
    main()