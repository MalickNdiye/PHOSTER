import pickle
import pandas as pd
import os
from drep import d_analyze
import matplotlib.pyplot as plt
import shutil


os.getcwd()
os.listdir("../results/reference_db/data_tables/")

# Open dRep tables

clust_info_p="../results/reference_db/data_tables/clust_info.tsv"
clust_assign_p="../results/reference_db/data_tables/clust_assign.tsv"
mtdata_p= "../data/metadata/RefGenomes_isolates_mtdata.csv"
to_delete_p="../results/reference_db/data_tables/to_del.tsv"
pickle_p="../results/reference_db/data/Clustering_files/"

clust_info=pd.read_csv(clust_info_p, delimiter="\t")
clust_assign=pd.read_csv(clust_assign_p, delimiter="\t")
mtdata=pd.read_csv(mtdata_p, delimiter="\t")
to_delete=pd.read_csv(to_delete_p, delimiter="\t")

# prepare dendogram
def plot_clust(file):
    prim_clust=file.split("/")[-1]
    prim_clust=prim_clust.split("_")[-1]
    prim_clust=int(prim_clust.split(".")[0])
    
    
        
    
    print("making dendogram for cluster " + str(prim_clust))
    
    clust_filt=clust_assign.query('primary_cluster == @prim_clust')
    
    gen_to_sp=clust_filt.set_index("genome").to_dict()["species"]
    gen_to_sdp=clust_filt.set_index("genome").to_dict()["SDP"]
    gen_to_ref=clust_filt.set_index("genome").to_dict()["ref"]
    
    
    f = open(file, 'rb')
    linkage = pickle.load(f)
    db = pickle.load(f)
    
    db_df=pd.DataFrame(db)
    
    
    genomes=list(db_df.index)
    SDP_genome=[gen_to_sdp[x] for x in genomes]
    species_genome=[gen_to_sp[x] for x in genomes]
    ref_genome=[gen_to_ref[x] for x in genomes]
    ref_genome = ["" if pd.isna(x) else x for x in ref_genome]
    
    
    names=[g + " (" + sdp + ") " for g,sdp in zip(genomes, SDP_genome)]
    names=[r + " " + n for r,n in zip(ref_genome, names)]
    names=[n + " " + s for n,s in zip(names,species_genome)]
    
    #plot fancy dendogram 
    fig, ax = plt.subplots(figsize=[16.5,11.7])
    dd=d_analyze.fancy_dendrogram(linkage, names, threshold=0.05, self_thresh=0.05)
    ax.set_title('cluster ' + str(prim_clust))
    
    f.close()
    return(fig)


file_list=os.listdir(pickle_p)
file_list=[f for f in file_list if "primary" not in f]

out_fg="../results/reference_db/figures/secondary_clusters_dendograms/"
if not os.path.exists(out_fg): 
    os.mkdir(out_fg)

for file in file_list:
    filename=file.strip(".pickle")
    dend=plot_clust(pickle_p + file)
    plt.savefig(out_fg+filename+".png", dpi=300, bbox_inches='tight')
    
# create filtered db
filt_db_path="../results/reference_db/dereplicated_genomes_filtered/"
old_db_p="../results/reference_db/dereplicated_genomes/"

if not os.path.exists(filt_db_path): 
    os.mkdir(filt_db_path)
    
gen_to_del=list(to_delete["genome"])


genome_list=os.listdir(old_db_p)
for entry in genome_list:
    if entry not in gen_to_del:
        shutil.copy(old_db_p+entry, filt_db_path+entry)


