#!/usr/bin/env python

import sys
import argparse
import re
import traceback
from itertools import groupby
from Bio import SeqIO
import pandas as pd
import os
import subprocess
import matplotlib.pyplot as plt
import time
import seaborn as sns
import shutil

__author__ = "Malick Ndiaye"
__title__ = "Align Viral Contigs"

# functions
# get command line arguments
def get_args():
    """Parse command line arguments"""
    desc = (
        """This script align viral contigs all vs all using fastANI. and generate a heatmap of ANI vs aligmnet fraction."""
    )
    epi = """returns directory with alignment files and a heatmap of ANI vs alignment fraction"""

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)

        # input
        parser.add_argument("-i",
            "--infasta", action="store",
             help='input fasta file'
        )

        # output
        parser.add_argument("-o",
            "--outdir",
            action="store",
            help="output directory"
        )

        # optional arguments
        parser.add_argument("-t",
            "--threads",
            action="store",
            help="threads",
            default=1,
            type=int,
            required=False
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


def split_contigs(infasta, dir):
    """split multifasta file into individual fasta files"""
    print("\nsplitting fasta file into individual fasta files")

    outdir = os.path.join(dir, "single_fasta_files")

    # create directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # split fasta file
    for seq_record in SeqIO.parse(infasta, "fasta"):
        filename = seq_record.id
        if len(seq_record.seq) > 10000:
            with open(os.path.join(outdir, filename + ".fasta"), "w") as output_handle:
                SeqIO.write(seq_record, output_handle, "fasta")


def create_fastANI_path_file(dir, outfile):
    """create a file with the path to all fasta files"""
    print("\ncreating a file with the path to all fasta files")

    # create a list of all absolute paths to fasta files
    fasta_files = []
    for root, dirs, files in os.walk(dir):
        for file in files:
            if file.endswith(".fasta"):
                fasta_files.append(os.path.join(root, file))

    # write list to file
    with open(outfile, "w") as f:
        for file in fasta_files:
            f.write(file + "\n")

    
def run_fastANI(fastANI_path_file, outdir, threads):
    """run fastANI"""

    # create directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # run fastANI
    cmd=["fastANI", "--ql", fastANI_path_file, "--rl", fastANI_path_file,  "-o", os.path.join(outdir, "fastANI_output.txt"), "-t", str(threads)]
    print("\nrunning fastANI")
    print("command:")
    print("\t" + " ".join(cmd))
    subprocess.run(" ".join(cmd), check=True, shell=True)

def format_fastani_file(fastani_file, outdir):
    """format fastANI output file for heatmap"""
    print("\nformatting fastANI output file for heatmap")

    # open file 
    df = pd.read_csv(fastani_file, sep="\t")

    # rename columns as suggested by fastANI
    df.columns = ["query", "reference", "ANI", "aligned_fragments", "total_fragments"]

    # query and reference are paths to fasta files, we only need the file names, without the path and extension
    df["query"] = df["query"].apply(lambda x: os.path.basename(x).split(".")[0])
    df["reference"] = df["reference"].apply(lambda x: os.path.basename(x).split(".")[0])

    # remove rows where query == reference
    df = df[df["query"] != df["reference"]]

    # crate a new column with aligmnet fraction (AF)
    df["AF"] = (df["aligned_fragments"] / df["total_fragments"])*100

    # remove rows with AF < 20
    df = df[df["AF"] >= 20]

    # remove rows with ANI<60
    df = df[df["ANI"] >= 60]

    # add wgANI column, which is (ANI*AF)/100
    df["wgANI"] = (df["ANI"] * df["AF"]) / 100

    # add same_cluster column for all rows where wgANI is >= 80.75
    df["same_cluster"] = df["wgANI"].apply(lambda x: True if x >= 80.75 else False)

    # save df to file named fastANI_output_tab.txt
    print("-->saving formatted fastANI output file in {}".format(os.path.join(outdir, "fastANI_output_formatted.txt")))
    df.to_csv(os.path.join(outdir, "fastANI_output_formatted.txt"), sep="\t")



def plot_figure(tab_file, outdir):
    """plot heatmap"""
    print("\nplotting Figure")

    #open tab file
    df_tab = pd.read_csv(tab_file, sep="\t", index_col=0)

    # take 100k random samples from df_tab
    if len(df_tab) > 100000:
        df_tab = df_tab.sample(n=100000, random_state=1)

    # plot scatter plot from tab file
    plt.figure(figsize=(40, 20))

    # add density plot with color gradient in function of density of points
    sns.kdeplot(data=df_tab, x="AF", y="ANI", thresh=0, fill=True, levels=100, cmap="magma", cbar=False)

    # add vertical line at x= 85 AF ymin=95 ymax=100
    plt.vlines(x=85, ymin=95, ymax=100, color="red", linewidth=5)

    # add diagonal line connecting the two vertical lines
    plt.plot([85, 100], [95, 80], color="red", linewidth=5)

    # add title
    plt.title("ANI vs AF Viral vMAGs", fontsize=20)

    # add x and y labels
    plt.xlabel("AF%", fontsize=20)
    plt.ylabel("ANI%", fontsize=20)

    # add x and y ticks every 5 points
    plt.xticks(range(20, 101, 5))
    plt.yticks(range(70, 101, 5))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)

    # save figure
    print("-->saving figure in {}".format(os.path.join(outdir, "fastANI_output_density.png")))
    plt.savefig(os.path.join(outdir, "fastANI_output_density.png"), dpi=100)


def main():
    """main function"""

    args=get_args()

    print("\n################################################################")
    print("arguments:")
    print("\tinput fasta file: " + args.infasta)
    print("\toutput directory: " + args.outdir)
    print("\tthreads: " + str(args.threads))
    print("################################################################\n")

    # split fasta file into individual fasta files
    split_contigs(args.infasta, args.outdir)

    # create a file with the path to all fasta files
    create_fastANI_path_file(os.path.join(args.outdir, "single_fasta_files"), os.path.join(args.outdir, "fastANI_path_file.txt"))

    # run fastANI
    run_fastANI(os.path.join(args.outdir, "fastANI_path_file.txt"), args.outdir, args.threads)

    #  wait for fastANI to finish
    print("waiting for fastANI to finish")
    while not os.path.exists(os.path.join(args.outdir, "fastANI_output.txt")):
        time.sleep(1)

    # remove single fasta files directory
    shutil.rmtree(os.path.join(args.outdir, "single_fasta_files"))

    # remove fastANI path file
    os.remove(os.path.join(args.outdir, "fastANI_path_file.txt"))
        
    # format fastANI output file for heatmap
    format_fastani_file(os.path.join(args.outdir, "fastANI_output.txt"), args.outdir)

    # plot heatmap
    plot_figure(os.path.join(args.outdir, "fastANI_output_formatted.txt"), args.outdir)

if __name__ == "__main__":
    main()

