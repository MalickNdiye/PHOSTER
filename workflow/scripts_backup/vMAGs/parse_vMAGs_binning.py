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
import difflib


__author__ = "Malick Ndiaye"
__title__ = "Parse_VMAGs_binning"


def get_args():
    """Parse command line arguments"""
    desc = (
        """This script integrate the binning files from vMAGs pipeline (vRhyme is the binning tool).
        """
    )
    epi = """Returns a fasta file of the representative sequences and a tabular file of the cluster membership"""
         

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # input files
        parser.add_argument("-i",
            "--infasta", action="store",
             help='input fasta file'
        )
        parser.add_argument("-b",
            "--binning_data",
            action="store",
            help="binning data as a tabular file"
        )
        parser.add_argument("-v",
            "--vi",
            action="store",
            help="table of viral identification scores"
        )
        parser.add_argument("-c",
            "--checkv",
            action="store",
            help="checkV output file"
        )

        # general options
        parser.add_argument("-m",
            "--minlength",
            action="store",
            help="minimum length of contigs to keep",
            default=3000
        )
        parser.add_argument("-a",
            "--maxlength",
            action="store",
            help="maximum length of contigs to keep",
            default=500000
        )
       
       # output files
        parser.add_argument("-f",
            "--outfasta", action="store",
             help='output concatenated fasta file'
        )
        parser.add_argument("-t",
            "--outtab", action="store",
             help='output filtering tabular file'
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


def rename_contigs(contig_list, sample): # TODO : rename contigs with unique names in a mor rubust way
    """Rename contigs with unique names"""

    # the new name of the contig will be sample_viral_contig_#
    # where # is the index of the contig in the list
    new_contig_list = []
    for i, contig in enumerate(contig_list):
        if "Bin" not in str(contig):
            new_contig_list.append("{}_viral_contig_{}".format(sample, i))
        else:
            new_contig_list.append(contig)
    
    # if some elements of new_contig_list are duplicated (it can happen if two prophages are dtected in the same contig), add "_" plus a number at the end of the name of all duplicates
    for i, contig in enumerate(new_contig_list):
        if new_contig_list.count(contig) > 1 and "Bin" not in str(contig):
            new_contig_list[i] = contig + "_" + str(new_contig_list[:i].count(contig))
    
    return(new_contig_list)


def parse_ckv(checkv, infasta):
    """parse checkV file to obtain new headers for the contigs"""
    print("\nParsing checkV file")

    # get fasta headers single column dataframe
    fasta_headers = pd.DataFrame([seq_record.id for seq_record in SeqIO.parse(infasta, "fasta")], columns=["trimmed_id"])
    # add empty column "contig_id"
    fasta_headers["contig_id"] = ""

    for i, element in fasta_headers["trimmed_id"].iteritems():
        if "vBin" in element:
            if len(element.split("_"))>3:
                contig_id = "_".join(element.split("_")[:-1])
            else:
                contig_id=element
        else:
            if "cutoff" in element and "parital" in element:
                if len(element.split("_"))>11:
                    contig_id = "_".join(element.split("_")[:-1])
                else:
                    contig_id=element
            elif "partial" in element:
                if len(element.split("_"))>7:
                    contig_id = "_".join(element.split("_")[:-1])
                else:
                    contig_id=element
            elif "cutoff" in element:
                if len(element.split("_"))>10:
                    contig_id = "_".join(element.split("_")[:-1])
                else:
                    contig_id=element
            elif len(element.split("_"))>6:
                contig_id = "_".join(element.split("_")[:-1])   
            else: 
                contig_id=element 

        # add contig_id to the dataframe
        fasta_headers.loc[i, "contig_id"] = contig_id

    # merge checkv and fasta_headers dataframes
    checkv = pd.merge(fasta_headers, checkv, how="left", on="contig_id")

    return(checkv)

    
def read_viral_mtdata(checkv_file, vi_file, binning_data, infasta):
    """Read viral metadata files and merge them"""
    print("\nReading viral metadata files")

    # create a dictionary with keys the bin and values the old contig names
    binning_data = binning_data[["bin", "old_contig_id"]]
    # rename column "bin" to "contig_id"
    binning_data = binning_data.rename(columns={"bin": "contig_id"})


    # get sample name
    sample = checkv_file.split("/")[-1].split("_")[0]
    print("\nprocessing sample: {}".format(sample))

    # read checkV file
    checkv = pd.read_csv(checkv_file, sep="\t", low_memory=False)
    checkv= parse_ckv(checkv, infasta)
    checkv["sample"] = sample

    # add old_contig_id column using the dictionary as the column "contig_id" is the same as "bin"
    checkv = pd.merge(checkv, binning_data, how="left", on="contig_id")
    # rename column "contig_id" to "bin"
    checkv = checkv.rename(columns={"contig_id": "bin"})
    # rename column "old_contig_id" to "contig_id"
    checkv = checkv.rename(columns={"old_contig_id": "contig_id"})

    # read vi file
    vi = pd.read_csv(vi_file, sep="\t", low_memory=False)
    vi = vi[vi["sample"] == sample]
   
    # merge columns "contig" anf "fragment" to get column "genome_id"
    # if fragment is "NaN", then contig_id = contig, else contig_id = fragment
    vi["contig_id"] = vi.apply(lambda x: x["contig"] if pd.isnull(x["fragment"]) else x["fragment"], axis=1)
    vi = vi.drop(["contig", "fragment"], axis=1)


    # merge checkv and vi by sample and "genome_id"
    viral_mtdata = pd.merge(checkv, vi, how="left", on=["sample", "contig_id"])
    viral_mtdata = viral_mtdata.rename(columns={"contig_id": "old_contig_id"})

    # merge columns "provirus" and "prophage" to get column "is_prophage"
    # if one of the two columns is "yes" or "Yes", then is_prophage = "yes", else is_prophage = "No"
    viral_mtdata["is_prophage"] = viral_mtdata.apply(lambda x: "yes" if x["provirus"] == "Yes" or x["prophage"] == "yes" else "no", axis=1)
    viral_mtdata = viral_mtdata.drop(["provirus", "prophage"], axis=1)
    # if is_prophage = "yes", then "type"=lysogenic, else type remains the same
    viral_mtdata["type"] = viral_mtdata.apply(lambda x: "lysogenic" if x["is_prophage"] == "yes" else x["type"], axis=1)
    # merge columns "proviral_length" and "frag_len" to get column "prophage_length"
    # if one both column are "NA", then prophage_length = "NA", if only one is "NA", then prophage_length = the other, else prophage_length = min(proviral_length, frag_len)
    viral_mtdata["prophage_length"] = viral_mtdata.apply(lambda x: x["proviral_length"] if pd.isnull(x["frag_len"]) else x["frag_len"] if pd.isnull(x["proviral_length"]) else min(x["proviral_length"], x["frag_len"]), axis=1)
    viral_mtdata = viral_mtdata.drop(["proviral_length", "frag_len"], axis=1)
    # merge columns "contig_length" and "len" to get column "contig_length"
    viral_mtdata["original_contig_length"] = viral_mtdata.apply(lambda x: x["len"] if x["contig_length"] == "NA" else x["contig_length"], axis=1)
    viral_mtdata = viral_mtdata.drop(["contig_length", "len"], axis=1)
    # merge columns "original_contig_length" and "prophage_length" to get column "final_contig_length"
    # if "prophage_length" is "NA", then final_contig_length = original_contig_length, else final_contig_length = prophage_length
    viral_mtdata["final_contig_length"] = viral_mtdata.apply(lambda x: x["original_contig_length"] if pd.isnull(x["prophage_length"]) else x["prophage_length"], axis=1)

    # reorder rows of viral_mtadata based on final_contig_length (descending order)
    viral_mtdata = viral_mtdata.sort_values(by=["final_contig_length"], ascending=False)

    # reorder columns of viral_mtadata 
    viral_mtdata=viral_mtdata[['sam_type', 'sam_name', 'sample', 'old_contig_id', "bin", "trimmed_id", 'original_contig_length', 'prophage_length',
     'final_contig_length', 'is_prophage', 'type', 'start', 'stop', 'quality', 'vv_quality', 'vs_quality', 'vibrant_quality',
      'vibrant_score', 'checkv_quality', 'miuvig_quality', 'completeness', 'completeness_method', 'Circular', 'max_score_group',
       'gene_count', 'viral_genes', 'host_genes', 'contamination', 'kmer_freq', "warnings"]]

    # rename column "qualuty" to "identification_tools_quality" and round values to 0 decimals
    viral_mtdata = viral_mtdata.rename(columns={"quality": "identification_tools_quality"})
    # if checkv_quality is Low-quality, then identification_tools_quality = (vv_quality + vs_quality + vibrant_quality + 1)/4
    # if checkv_quality is Medium-quality, then identification_tools_quality = (vv_quality + vs_quality + vibrant_quality + 2)/4
    # if checkv_quality is High-quality or Complete , then identification_tools_quality =  (vv_quality + vs_quality + vibrant_quality + 3)4
    # else identification_tools_quality = identification_tools_quality/4
    viral_mtdata["identification_tools_quality"] = viral_mtdata.apply(lambda x: (x["vv_quality"] + x["vs_quality"] + x["vibrant_quality"] + 1)/4 if x["checkv_quality"] == "Low-quality" else (x["vv_quality"] + x["vs_quality"] + x["vibrant_quality"] + 2)/4 if x["checkv_quality"] == "Medium-quality" else (x["vv_quality"] + x["vs_quality"] + x["vibrant_quality"] + 3)/4 if x["checkv_quality"] == "High-quality" or x["checkv_quality"] == "Complete" else x["identification_tools_quality"]/4, axis=1)

    viral_mtdata["identification_tools_quality"] = viral_mtdata["identification_tools_quality"].round(0)

    return(viral_mtdata)


def merge_binning_to_mtdata(binning, viral_mtdata):
    """Merge binning and viral_mtdata"""
    print("\nMerging binning and viral_mtdata")

    # merge binning and viral_mtdata by sample and "genome_id"
    merged = pd.merge(viral_mtdata, binning, how="left", on=["sample", "old_contig_id", "bin"])

    # if bin is empty, then bin="dereplicated"
    merged["bin"] = merged.apply(lambda x: "dereplicated" if pd.isnull(x["bin"]) else x["bin"], axis=1)

    # add column genome_id by renameing column "bin" using function rename_contigs
    merged["genome_id"] = rename_contigs(merged["bin"].to_list(), merged["sample"].tolist()[0])

    # add column bin_length wich is the sum of the final_contig_length of all contigs in the bin
    merged["bin_length"] = merged.groupby("bin")["final_contig_length"].transform("max")

    # rename column "members" to "bin_members"
    merged = merged.rename(columns={"members": "bin_members"})
    # rename column "proteins" to "bin_proteins"
    merged = merged.rename(columns={"proteins": "bin_proteins"})
    # rename column "redundancy" to "bin_redundancy"
    merged = merged.rename(columns={"redundancy": "bin_redundancy"})

    # change identification_tools_quality elements to max of the values in the bin
    merged["identification_tools_quality"] = merged.groupby("bin")["identification_tools_quality"].transform("max")

    # reorder columns of merged
    merged=merged[['sam_type', 'sam_name', 'sample', 'old_contig_id',  "bin", "trimmed_id", "genome_id", 'original_contig_length', 'prophage_length',
     'final_contig_length', "bin_length", 'is_prophage', 'type', 'start', 'stop', 'identification_tools_quality', 'vv_quality', 'vs_quality', 'vibrant_quality',
      'vibrant_score', 'checkv_quality', 'miuvig_quality', 'completeness', 'completeness_method', 'Circular', 'max_score_group',
       'gene_count', 'viral_genes', 'host_genes', 'contamination', "bin_members", "bin_proteins", "bin_redundancy", 'kmer_freq', "warnings"]]

    return(merged)

    
def filter_contigs_tab(viral_mtdata, maxlength, minlength):
    """filter contigs based on quality and length"""
    print("\nFiltering contigs metadata based on quality and length in table")

    viral_mtadata_filtered = viral_mtdata.copy()

    # add empty columns "retained" and "reason" to viral_mtdata
    viral_mtadata_filtered["retained"] = ""
    viral_mtadata_filtered["reason"] = ""

    # filter contigs based on quality
    viral_mtadata_filtered["retained"] = viral_mtadata_filtered.apply(lambda x: "no" if x["identification_tools_quality"] < 1 else "yes", axis=1)
    viral_mtadata_filtered["reason"] = viral_mtadata_filtered.apply(lambda x: "quality" if x["identification_tools_quality"] < 1 else x["reason"], axis=1)

    # filter cotigs based on dereplicated
    viral_mtadata_filtered["retained"] = viral_mtadata_filtered.apply(lambda x: "no" if x["bin"] == "dereplicated" else x["retained"], axis=1)
    viral_mtadata_filtered["reason"] = viral_mtadata_filtered.apply(lambda x: "dereplicated" if x["bin"] == "dereplicated" else x["reason"], axis=1)

    # filter contigs based on length
    viral_mtadata_filtered["retained"] = viral_mtadata_filtered.apply(lambda x: "no" if x["bin_length"] >= maxlength or x["bin_length"] < minlength else x["retained"], axis=1)
    viral_mtadata_filtered["reason"] = viral_mtadata_filtered.apply(lambda x: "length" if x["bin_length"] >= maxlength or x["bin_length"] < minlength else x["reason"], axis=1)

    return(viral_mtadata_filtered)

def filter_contigs_fasta(fasta, viral_mtdata_filtered):
    """filter fasta file based on filtered metadata"""
    print("\nFiltering fasta file based on filtered metadata")

    # filter vital_mtdata_filtered based on "retained" column
    viral_mtdata_filtered = viral_mtdata_filtered[viral_mtdata_filtered["retained"] == "yes"]

    # open fasta file as list using biopython
    fasta = list(SeqIO.parse(fasta, "fasta"))

    # reorder fasta file based on length of contigs (descending)
    fasta = sorted(fasta, key=lambda x: len(x.seq), reverse=True)

    # filter fasta file based on "old_contig_id" column
    fasta_filtered = [x for x in fasta if x.id in viral_mtdata_filtered["trimmed_id"].tolist()]

    # rename header of fasta file based on "genome_id" column
    for contig in fasta_filtered:
        contig.id = viral_mtdata_filtered[viral_mtdata_filtered["trimmed_id"] == contig.id]["genome_id"].tolist()[0]
        contig.description = contig.id

    # some sequences in fasta_filtered have the same id, we need to concatenate the sequences

    # create empty list to store fasta sequences
    fasta_filtered_concat = []

    # create empty list to store fasta sequences ids
    fasta_filtered_concat_ids = []

    # loop over fasta_filtered
    for contig in fasta_filtered:
            # if contig id is not in fasta_filtered_concat_ids, append contig to fasta_filtered_concat and append contig id to fasta_filtered_concat_ids
            if contig.id not in fasta_filtered_concat_ids:
                fasta_filtered_concat.append(contig)
                fasta_filtered_concat_ids.append(contig.id)
    
            # if contig id is in fasta_filtered_concat_ids, concatenate contig sequence to fasta_filtered_concat
            else:
                fasta_filtered_concat[fasta_filtered_concat_ids.index(contig.id)].seq += contig.seq

    fasta_filtered=fasta_filtered_concat



    return(fasta_filtered)

def main():
    """Call functions"""

    args = get_args()
   
    print("\n\n################################################")
    print("running script: " + os.path.basename(__file__))
    print("input fasta file: " + args.infasta)
    print("input binning membership directory: " + args.binning_data)
    print("input contigs metadata table: " + args.vi)
    print("input checkV output file: " + args.checkv)
    print("output filtered fasta file: " + args.outfasta)
    print("output filtered contigs metadata table: " + args.outtab)
    print("################################################")

    binning_data=pd.read_csv(args.binning_data, sep="\t")

    print(binning_data.head())

    viral_mtdata=read_viral_mtdata(args.checkv, args.vi, binning_data, args.infasta)
    merged=merge_binning_to_mtdata(binning_data, viral_mtdata)

    viral_mtdata_filtered=filter_contigs_tab(merged, args.maxlength, args.minlength)

    fasta_filtered=filter_contigs_fasta(args.infasta, viral_mtdata_filtered)

    # write filtered fasta file
    print("\nWriting filtered fasta file...")
    if not os.path.exists(os.path.dirname(args.outfasta)):
        os.makedirs(os.path.dirname(args.outfasta))
    SeqIO.write(fasta_filtered, args.outfasta, "fasta")

    # write filtered metadata table
    print("\nWriting filtered metadata table...")
    if not os.path.exists(os.path.dirname(args.outtab)):
        os.makedirs(os.path.dirname(args.outtab))
    viral_mtdata_filtered.to_csv(args.outtab, sep="\t", index=False)


    # fiöter vrial_mtdata based on "retained" column
    viral_mtdata_check = viral_mtdata_filtered[viral_mtdata_filtered["retained"] == "yes"]

    #check that all genome_id of viral_mtdata_check are in fasta_filtered
    # if not print the missing genome_id
    if len(set(viral_mtdata_check["genome_id"].tolist())) != len(fasta_filtered):
        print("\n\n################################################")
        print("ERROR: not all genome_id of viral_mtdata_check are in fasta_filtered")
        print("missing genome_id:")
        print(set(viral_mtdata_check["genome_id"].tolist()) - set([x.id for x in fasta_filtered]))
        print("################################################\n\n")
        sys.exit(1)


    
    

    print("\n\n################################################")
    print("DONE!")
    print("################################################\n\n")


if __name__ == "__main__":
    main()
        