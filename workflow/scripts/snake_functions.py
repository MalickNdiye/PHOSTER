import pandas as pd
import numpy as np
import os

def get_direction_r1(wildcards) :
#This  function returns a list of R1 reads when the same sample was sequenced on multiple lanes
    samname=wildcards.sample
    l=config["samples"][samname]

    r_list=[]
    r1_l1=[s for s in l if "L1_R1" in s]
    r1_l2=[s for s in l if "L2_R1" in s]
    r1_l3=[s for s in l if "L3_R1" in s]
    r1_l4=[s for s in l if "L4_R1" in s]

    r_list.extend(r1_l1)
    r_list.extend(r1_l2)
    if len(r1_l3) >0: r_list.extend(r1_l3)
    if len(r1_l4) >0: r_list.extend(r1_l4)
    return(r_list)

def get_direction_r2(wildcards) :
#This  function returns a list of R2 reads when the same sample was sequenced on multiple lanes
    samname=wildcards.sample
    l=config["samples"][samname]

    r_list=[]
    r2_l1=[s for s in l if "L1_R2" in s]
    r2_l2=[s for s in l if "L2_R2" in s]
    r2_l3=[s for s in l if "L3_R2" in s]
    r2_l4=[s for s in l if "L4_R2" in s]

    r_list.extend(r2_l1)
    r_list.extend(r2_l2)
    if len(r2_l3) >0: r_list.extend(r2_l3)
    if len(r2_l4) >0: r_list.extend(r2_l4)
    return(r_list)

def get_files_commas(path, sep=",", remove_hidden=True):
    file_l=os.listdir(path)

    if remove_hidden:
        file_l=[f for f in file_l if not f.startswith(".")]

    file_l2=[]
    for f in file_l:
        file_l2.append(os.path.join(path,f))
    out=sep.join(file_l2)
    return(out)

def get_wildcard_commas(paths, sep=","):
    file_l=[paths]
    out=sep.join(file_l)
    return(out)

def get_genomes(path, ncol=1):
    df=pd.read_csv(path, delimiter="\t")
    genomes=np.array(df.iloc[:,ncol])
    genomes=np.unique(genomes)

    return(genomes)

def aggregate_viral_ids():
    sample_table = checkpoints.split_vamb_contigs.get(**wildcards).output[3]
    IDS = []
    with open(sample_table,'r') as infile:
        for line in infile:
            line = line.rstrip()
            IDS.append(line)
    return(IDS)

def get_vids(wildcards):
    sample_table = checkpoints.split_vamb_contigs.get(**wildcards).output[3]
    IDS = []
    with open(sample_table,'r') as infile:
        for line in infile:
            line = line.rstrip()
            IDS.append(line)

    return(expand(["../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}_dvf/{vSam}.fna_gt2000bp_dvfpred.txt", "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.hmmMiComplete105.tbl", "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.hmmVOG.tbl"], vSam=IDS))

def get_all_SecCluster(path, ignore_sing=True):
    checkpoint_output = path
    df=pd.read_csv(path, delimiter="\t")
    secodary_clusters=list(df["secondary_cluster"])

    if ignore_sing:
        secodary_clusters=set([i for i in secodary_clusters if secodary_clusters.count(i)>1])
    else:
        secodary_clusters=set(secodary_clusters)

    return(secodary_clusters)


def concat_tables(tables_list, skip=0, delimiter="\t", head=0, add={}):
    # this function concatenates a list of tables
    # tables_list: list of tables to concatenate
    # returns: concatenated table
    df_list=[]
    count=0
    for t in tables_list:
        # open file skipping "skip lines" and adding header
        df=pd.read_csv(t,  skiprows=skip, sep=delimiter, header=head)
         
         # "add" is a dictionary with columns to add
         # hey: column name, value: list of values to add, each value correspond to a daframe in tables_list
        for key, value in add.items():
            df[key]=value[count]
        count+=1

        df_list.append(df)

    df=pd.concat(df_list)
    
    return(df)


def get_MAGs_red(wildcards):
    clust_file = checkpoints.parse_dRep.get(**wildcards).output[3]
    df=pd.read_csv(clust_file, delimiter="\t")

    # get column "genome" as list
    genomes=list(df["genome"])

    return(genomes)

def get_bact_annot_prot(wildcards):
    dir=checkpoints.annotate_bacterial_genomes.get(**wildcards).output[1]
    # list fines in dir keeping entire path
    files=os.listdir(dir)
    files=[os.path.join(dir, f) for f in files]

    return(files)

def filter_fasta(fasta, list_fasta, out):
    # this function filters a fasta file keeping only sequences in list
    # fasta: fasta file to filter
    # list: list of sequences to keep
    # out: output file
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    # read fasta file
    records = list(SeqIO.parse(fasta, "fasta"))

    # filter records
    records_filt = [r for r in records if r.id in list_fasta]

    # write filtered records
    SeqIO.write(records_filt, out, "fasta")

