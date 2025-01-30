import os
import sys
from Bio import SeqIO

length_threshold = snakemake.params["length_t"]
coverage_threshold = snakemake.params["cov_t"]
f_in = snakemake.input[0]
f_out = snakemake.output[0]
sam_pos=0

print("Filtering contigs from " + f_in + " to " + f_out)

#get sample name
filename=os.path.basename(f_in)
sample=filename.split(".")[sam_pos]


fasta_sequences = SeqIO.parse(open(f_in),'fasta')
filt_conts=open(f_out, "a")

count=1
for fasta in fasta_sequences:
    name, sequence = fasta.description, str(fasta.seq)
    length=len(fasta.seq)

    print("\n--- Contig " + str(count))
    count+=1

    if length > float(length_threshold):
        SeqIO.write(fasta, filt_conts, "fasta")
        
filt_conts.close()
