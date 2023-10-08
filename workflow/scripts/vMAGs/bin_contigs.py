#!/usr/bin/env python

import sys
import argparse
import re
import traceback
from itertools import groupby
from Bio import SeqIO
import pandas as pd
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


__author__ = "Malick Ndiaye"
__title__ = "Bin vRhyme MAGs"


def get_args():
    """Parse command line arguments"""
    desc = (
        """This script bins vRhyme contigs into vMAGs based on the cluster membership of the representative sequences.
        """
    )
    epi = """Returns a fasta file of the Binned Contigs sequences and a tabular file of the cluster membership"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # input files
        parser.add_argument("-i",
            "--infasta", action="store",
             help='input fasta file'
        )
        parser.add_argument("-d",
            "--vrhyme_dir",
            action="store",
            help="output vRhyme directory"
        )
        
       # output files
        parser.add_argument("-f",
            "--outfasta", action="store",
             help='output concatenated vMAGs fasta file'
        )
        parser.add_argument("-t",
            "--outtab", action="store",
             help='output membership table file'
        )


        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occurred with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return(parser.parse_args())


def get_dereplicated_contigs(dir):
    """Read fasta file and return a list of dereplicated contigs"""
    # get semple name from dir name
    sample_name=dir.split("/")[-1]
    if sample_name=="":
        sample_name=dir.split("/")[-2]

    drep_dir=os.path.join(dir,"vRhyme_dereplication")
    # find file ending with "fa" or "fasta" in dir
    drep_list=[x for x in os.listdir(drep_dir) if x.endswith(".fa") or x.endswith(".fasta")]
    if len(drep_list)==0:
        print("No dereplicated fasta file for sample {}".format(sample_name))
        return([])
    drep_file=drep_list[0]

    # read drepliacted fasta file
    dereplicated_fasta=os.path.join(drep_dir,drep_file)

    # read fasta file
    records = list(SeqIO.parse(dereplicated_fasta, "fasta"))

    # get contig ids in list
    contigs = [x.description for x in records]

    return(contigs)


def get_binning_data(dir):
    """Read binning data from vRhyme output directory"""

    # get semple name from dir name
    sample_name=dir.split("/")[-1]
    if sample_name=="":
        sample_name=dir.split("/")[-2]

    # read file ending with membership.tsv in dir
    print("\nReading binning data for sample {}".format(sample_name))
    print(os.listdir(dir))

    mem_list=[x for x in os.listdir(dir) if "membership.tsv" in x]
    
    if len(mem_list)==0:
        print("No binning data for sample {}".format(sample_name))
        #return empt dataframe with columns sample, contig, bin, members, proteins, redundancy 
        membership=pd.DataFrame(columns=["sample","contig", "bin"])
    else:
        mem_file=mem_list[0]
        mem_path=os.path.join(dir,mem_file)
        membership=pd.read_csv(mem_path,sep="\t", low_memory=False)
        membership=membership.assign(sample=sample_name)

        # move sample to first column and rename scaffodl to contig
        membership=membership[["sample","scaffold", "bin"]]
        membership=membership.rename(columns={"scaffold":"contig"})

    # read file ending with summary.tsv in dir
    sum_list=[x for x in os.listdir(dir) if "summary.tsv" in x]
    if len(sum_list)==0:
        print("No binning data for sample {}".format(sample_name))
        #return empt dataframe with columns sample, contig, bin, members, proteins, redundancy 
        summary=pd.DataFrame(columns=["bin","members","proteins","redundancy"])
    else:
        sum_file=sum_list[0]
        summary=pd.read_csv(os.path.join(dir,sum_file),sep="\t", low_memory=False)

    # merge membership and summary
    joined=pd.merge(membership,summary,on="bin")
    # rename column "contig" to "old_contig_id"
    joined=joined.rename(columns={"contig":"old_contig_id"})
    # bin is a float, convert to str with no decimal
    joined["bin"]=joined["bin"].astype(int).astype(str)
    # rename bin as sample_vBin_bin
    joined["bin"]=joined["sample"]+"_vBin_"+joined["bin"].astype(str)
    
    # get dereplicated contigs
    all_contigs=get_dereplicated_contigs(dir)
    all_contigs=pd.DataFrame(all_contigs,columns=["old_contig_id"])

    # merge all_contigs and joined
    joined_final=pd.merge(all_contigs,joined,on="old_contig_id",how="left")
    # fill NA values in bin column with old_contig_id
    joined_final["bin"]=joined_final["bin"].fillna(joined_final["old_contig_id"].astype(str))
    # fill NA values in sample column with sample_name
    joined_final["sample"]=joined_final["sample"].fillna(sample_name)
    # fill NA values in members column with 1
    joined_final["members"]=joined_final["members"].fillna(1)

    joined_final=joined_final[["sample","old_contig_id","bin","members","proteins","redundancy"]]
    

    return(joined_final)

def bin_contigs(fasta, binning_data, derplicated_contigs):
    """create a SeqIo object of binned contigs"""
    print("\nConcatenating fasta files based on bin")

    # get sample name from fasta file name
    sample_name=fasta.split("/")[-1].split(".")[0]

    # read fasta file
    records = list(SeqIO.parse(fasta, "fasta"))

    # remove contigs that are not in dereplicated contigs
    records_drep=[x for x in records if x.description in derplicated_contigs]

    # get contig ids in list
    contigs = [x.description for x in records_drep]

    # empty dictionary to store binned contigs
    binned_contigs={}
  
    # loop through records_drep, if contig in binning_data, add to binned_contigs, headr as key and sequence as value
    # if header is already in binned_contigs, concatenate sequence
    for record in records_drep:
        if record.description in binning_data["old_contig_id"].tolist():
            # get bin name
            bin_name=binning_data[binning_data["old_contig_id"]==record.description]["bin"].tolist()[0]
            # add to binned_contigs
            if bin_name in binned_contigs.keys():
                binned_contigs[bin_name]=binned_contigs[bin_name]+record.seq
            else:
                binned_contigs[bin_name]=record.seq

    # convert binned_contigs to a list of SeqRecord objects
    binned_contigs_record=[SeqRecord(seq, id=header, description="") for header,seq in binned_contigs.items()]

    return(binned_contigs_record)

def main():
    """call all functions and write output files"""

    args=get_args()

    print("\n##############################################################")
    print("vRhyme binning data to vMAGs")
    print("input fasta file: {}".format(args.infasta))
    print("input vRhyme output directory: {}".format(args.vrhyme_dir))
    print("output binned fasta: {}".format(args.outfasta))
    print("output bin memebrship table: {}".format(args.outtab))
    print("##############################################################")

    binnig_data=get_binning_data(args.vrhyme_dir)
    binned_contigs=bin_contigs(args.infasta,binnig_data,get_dereplicated_contigs(args.vrhyme_dir))

    # write binned fasta, if path does not exist, create it
    if not os.path.exists(os.path.dirname(args.outfasta)):
        os.makedirs(os.path.dirname(args.outfasta))
    print("\nWriting binned fasta to {}".format(args.outfasta))
    SeqIO.write(binned_contigs,args.outfasta,"fasta")

    # write bin membership table, if path does not exist, create it
    if not os.path.exists(os.path.dirname(args.outtab)):
        os.makedirs(os.path.dirname(args.outtab))
    print("\nWriting bin membership table to {}".format(args.outtab))
    binnig_data.to_csv(args.outtab,sep="\t",index=False)

    print("\n##############################################################")
    print("Done!")
    print("##############################################################")

if __name__ == '__main__':
    main()
