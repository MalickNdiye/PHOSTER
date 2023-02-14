import os
import pandas as pd
import numpy as np
import subprocess

# Set global variables
genus=snakemake.params["genus"]
outdir=snakemake.output["pan"]
initial_p=snakemake.input["faas"]
threads=snakemake.threads
bootstrap=snakemake.params["boots"]
clust_p=snakemake.input["clust_info"]
checkm_df=snakemake.input["checkm_out"]

# open clustering info and look for genomes that correspond to the wanted genus
clust_df=pd.read_csv(clust_p, delimiter="\t")
filt_clust=clust_df.query("genus==@genus")

# Get all clusters for a given genome
cluster_list=np.unique(filt_clust["secondary_cluster"])

# Loop through the clusters and launch motupan for the genomes in the cluster
for cluster in cluster_list:
    print("\n----Finding genomes for cluster " + cluster)
    df_sec_clust=filt_clust.query("secondary_cluster==@cluster")
    genomes_list=df_sec_clust["Bin Id"]

    path_list=[os.path.join(initial_p, g + "_genes.faa") for g in genomes_list]
    outfile=os.path.join(outdir, cluster + "_mOTUpan.tsv")

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    cmd=["mOTUpan.py", "--faas", " ".join(path_list), "-o", outfile, "--checkm", checkm_df, "--boots", str(bootstrap), "--threads", str(threads)]
    print("----Launching mOTUpan with the following command: " + " ".join(cmd))
    subprocess.run(" ".join(cmd), shell=True)
