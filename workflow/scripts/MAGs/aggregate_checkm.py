import pandas as pd
import os
import shutil


################################################################################
################################################################################
# This scripts take as input a list of check stats tables (--tab-file) and a
# directory where several subdirectories containing MAGs from multiple samples
# are stored. It parses the stats tables to chose "good" MAGs (completeness > 50%
# and contamination < 10%). Then, it return a combined stat table for all MAGs
# and one for the good MAGs. Also copies all "good" MAGs in a single directory.

#inputs
stats_list=list(snakemake.input["stats"])
bins_dir=snakemake.input["bins"]

#outputs
out_full_stats=snakemake.output["full_stats"] # checkm stats of all MAGs
out_filtered_stats=snakemake.output["filtered_stats"] # checkm stats of "good" MAGs
chosen_mags_dir=snakemake.output["filered_mags"] # directory where all good MAGs will be stored


def create_df_list(list):
# opens dataframes in file list and creates list of df
    df_list=[]
    for file in list:
        sam=file.split("_")[0]
        sam=sam.split("/")[-1]
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
filtred_df=full_df.query("Completeness>50 and Contamination<10") # MAGs with completeness < 50% and contamination > 10% are not even considered for dereplication by dRep

full_df.to_csv(out_full_stats, sep="\t", index=False)
filtred_df.to_csv(out_filtered_stats, sep="\t", index=False)

bins_sam=os.listdir(bins_dir)
all_filtered_mags=[]
gd_mags=list(filtred_df["Bin Id"])

for dir in bins_sam:
    inpath=os.path.join(bins_dir, dir)
    new_mags=get_good_mags(inpath, gd_mags)
    for entry in new_mags:
        source=entry
        basename=os.path.basename(entry)
        destination=os.path.join(chosen_mags_dir,basename)

        if not os.path.exists(chosen_mags_dir):
            os.makedirs(chosen_mags_dir)

        shutil.copy(source, destination)
