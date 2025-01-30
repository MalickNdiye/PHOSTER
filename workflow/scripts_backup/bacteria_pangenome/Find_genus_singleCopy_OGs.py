import os
import pandas as pd

os.getcwd()
target_genus=snakemake.params["genus"]

print("########################################################################")
print("\tFinding single copy orthiologues for all the genomes in the genus" + target_genus)
print("########################################################################")


# Set paths
og_tab=os.path.join(snakemake.input["og"], "Results_results/Orthogroups/Orthogroups.tsv")
clust_info_p=snakemake.input["clust_info"]


# Open clusering info
clust_info=pd.read_csv(clust_info_p, delimiter="\t")
clust_info=clust_info.assign(genus=[f.split(" ")[0] for f in clust_info.species]) #TODO Remove after code update

## Select only isolates
clust_info_isolate=clust_info.query("isolate=='yes'").query("genus==@target_genus")

# Format Orthofider Single copy OGs data
og_df=pd.read_csv(og_tab, delimiter="\t")


og_df_melt=og_df.melt(id_vars="Orthogroup",var_name="genome", value_name="gene").dropna()
og_df_melt=og_df_melt.assign(length_og=[len(f) for f in og_df_melt.gene.str.split(",")]).query("length_og==1")


og_df_melt["genome"]=[f.strip("_genes") for f in og_df_melt["genome"]]
og_df_melt=og_df_melt.drop(columns=["length_og"])

clust_add=clust_info[["Bin Id", "secondary_cluster", "species", "genus", "secondary_cluster_n"]]

og_df_melt=og_df_melt.merge(clust_add, how="left", left_on="genome", right_on="Bin Id")

## select OG only from isolates
og_df_isolates=og_df_melt.query("genome in @clust_info_isolate['Bin Id']")

og_df_melt.to_csv(snakemake.output["genus_OGs"], sep="\t", index=False)
og_df_isolates.to_csv(snakemake.output["isolates_OGs"], sep="\t", index=False)
