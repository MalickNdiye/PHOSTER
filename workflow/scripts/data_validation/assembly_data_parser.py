from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os

## This scripts takes the contigs files of the assembly of the reads that were
## unmapped by bbsplit and returns a table containing, foir each sample, the
## contig number, its size and the k_mer coverage

outfile1=snakemake.output[0]
outfile2=snakemake.output[1]
file_list=snakemake.input[:]  # list of the contigs files of each sample

# start some list that will be the columns of the final dataframe
samples=[]
contigs=[]
lengths=[]
covs=[]

IDS=[]
seqs=[]



# loop through each file and recover all the info from the headers
for file in file_list:

    sample=file.split("/")[-1]
    sample=sample.split("_")[0]

    print("processing " + sample + " assembly..." )

    i=0
    for record in SeqIO.parse(file, "fasta"):
        if (len(record.seq) >= 1000) and ("P" in sample):
            id_f=sample + "_" + record.id
            IDS.append(id_f)
            seqs.append(record.seq)

        id=record.id
        s_id=id.split("_")

        node_id=s_id[1]
        size=s_id[3]
        cov=s_id[5]

        print("contig: " + node_id + ", size[bp]: " + size + ", kmer coverage: " + cov)

        contigs.append(node_id)
        lengths.append(size)
        covs.append(cov)
        i+=1

    sams= [sample] * i
    samples.extend(sams)

rec_list=[]
for i in range(0, len(IDS)):
    record = SeqRecord(Seq(seqs[i]),id=IDS[i])
    record.description=""
    print("writing " + record.id + " in " + outfile2)
    rec_list.append(record)

SeqIO.write(rec_list, outfile2, "fasta")

d={"sample":samples, "contig":contigs, "size[bp]": lengths, "kmer_cov": covs}
final_df= pd.DataFrame(data=d)


final_df.to_csv(outfile1, index=None, sep='\t')
