from Bio import SeqIO
import os
import pandas as pd

#open Orthofinder output
cluster=snakemake.params["clust"]
print("#####################################################################")
print("\t\tWORKING ON CLUSTER " + cluster + "\n")
print("#####################################################################")

# Open all Instrain Profiles
IS_list=list(snakemake.input["IS_profile"])


IS_genes_files=[]
for prof in IS_list:
    IS_profiles_dir=os.path.normpath(prof+ "/output")
    IS_genes_f=os.listdir(IS_profiles_dir)
    IS_genes_f=[f for f in IS_genes_f if "gene_info.tsv" in f]
    IS_genes_f=[f for f in IS_genes_f if not f.startswith(".")]
    IS_genes_f=[IS_profiles_dir + "/" + f for f in IS_genes_f]
    IS_genes_files.extend(IS_genes_f)

print("profiles: ")
print(IS_genes_files)

# open scaffold to genome file
stb=snakemake.input["stb"]

stb_df=pd.read_csv(stb, delimiter="\t",  header=None)
stb_df.columns=["scaffold", "genome"]


single_copy_dir=os.path.normpath(snakemake.input["og"] + "/Results_results/Single_Copy_Orthologue_Sequences")
single_copy_f=os.listdir(single_copy_dir)
OG_path=[os.path.join(single_copy_dir,f) for f in single_copy_f]

#set outfile
out_tab=snakemake.output[0]


og_df=pd.DataFrame(columns=["secondary_cluster", "OG", "gene"])
for og in OG_path:
    og_id=og.split("/")[-1]
    og_id=og_id.split(".")[0]
    print("---parsing OG: " + og_id)
    for record in SeqIO.parse(og, 'fasta'):
        header=record.id
        new_row=[cluster, og_id, header]
        og_df.loc[len(og_df)] = new_row



single_copy_df=pd.DataFrame()
for f in IS_genes_files:
    print(f)
    sam=f.split("/")[-1]
    sam=sam.split("_")[0]
    print("Recovering Signle copy orthologues in sample: " + sam)

    path=f
    profile=pd.read_csv(path, skiprows=0, delimiter="\t")
    profile["sample"]=[sam]*profile.shape[0]

    a=profile.merge(og_df, how="inner",on="gene")
    a=a.merge(stb_df,how="inner", on="scaffold")
    single_copy_df=pd.concat([single_copy_df, a])

print("saving summmary table in: " + out_tab)
single_copy_df.head()
single_copy_df.to_csv(out_tab, sep="\t", index=False)
