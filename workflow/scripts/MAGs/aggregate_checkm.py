import pandas as pd
import os
import shutil
import sys
import os
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="""this script takes as input a list of check stats tables (--tab-file) and the corresponding bins,
     cleans the table and filters the bins based on completeness and contamination. It returns a combined stat table for all MAGs and one for the good MAGs. 
     Also copies all "good" MAGs in a single directory. The script is called by the snakemake workflow. The input and output files are defined in the snakemake file. 
     The script can also be run independently by providing the input and output files as arguments. The input files are:""",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)

parser.add_argument("-i", "--input_stats", help="Paths to input checkm stat files", nargs="+")
parser.add_argument("-b", "--bins", help="paths to bin files", nargs="+")
parser.add_argument("-o", "--full_stats", help="paths to full stat output")
parser.add_argument("-f", "--filt_stats", help="paths to filtered stat output")
parser.add_argument("-m", "--mags", help="paths to filtered MAGs directory")
parser.add_argument("-c", "--completeness", help="completeness threshold")
parser.add_argument("-e", "--contamination", help="contamination threshold")

if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
    parser.print_help()
    sys.exit()

args = parser.parse_args()


################################################################################
################################################################################
# This scripts take as input a list of check stats tables (--tab-file) and a
# directory where several subdirectories containing MAGs from multiple samples
# are stored. It parses the stats tables to chose "good" MAGs (completeness > 50%
# and contamination < 10%). Then, it return a combined stat table for all MAGs
# and one for the good MAGs. Also copies all "good" MAGs in a single directory.

#inputs
stats_list=args.input_stats

#outputs
out_full_stats=args.full_stats # checkm stats of all MAGs
out_filtered_stats=args.filt_stats # checkm stats of "good" MAGs
chosen_mags_dir=args.mags # directory where all good MAGs will be stored


def create_df_list(list):
# opens dataframes in file list and creates list of df
    df_list=[]
    for file in list:
        sam=file.split("/")[-1]
        sam=sam.split("_")[0]
        df=pd.read_csv(file, delimiter="\t")

        df.insert(0, "sample", sam, True)


        df_list.append(df)

    return(df_list)

def get_good_mags(dir, good_mags):
# given a directory with subfirectories containings MAGs, look for fasta files
# that correspond to an input list of identifiers of good MAGs.
# Returns the paths to the good MAGs
    filelist=os.listdir(dir)
    good_mags_paths=[]

    for file in filelist:
        if os.path.basename(file).strip(".fa") in good_mags:
            entry= os.path.abspath(os.path.join(dir,file))
            path_mag=entry
            print(path_mag)
            good_mags_paths.append(path_mag)

    return(good_mags_paths)


dfs=create_df_list(stats_list)
full_df=pd.concat(dfs)
filter="Completeness>"+str(args.completeness)+ " and Contamination<" +str(args.contamination)
filtred_df=full_df.query(filter) # MAGs with completeness < 50% and contamination > 10% are not even considered for dereplication by dRep

full_df.to_csv(out_full_stats, sep="\t", index=False)
filtred_df.to_csv(out_filtered_stats, sep="\t", index=False)

bins_sam=args.bins
all_filtered_mags=[]
gd_mags=list(filtred_df["Bin Id"])

for dir in bins_sam:
    print("looking for good MAGs in " + dir)
    inpath=dir
    new_mags=get_good_mags(inpath, gd_mags)
    for entry in new_mags:
        source=entry
        basename=os.path.basename(entry)
        destination=os.path.join(chosen_mags_dir,basename)

        if not os.path.exists(chosen_mags_dir):
            os.makedirs(chosen_mags_dir)

        shutil.copy(source, destination)
