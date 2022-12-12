import os
import sys
from Bio import SeqIO

length_threshold = snakemake.params["length_t"]
coverage_threshold = snakemake.params["cov_t"]
f_in = snakemake.input[0]
f_out = snakemake.output[0]
tab_out= snakemake.output[1]
sam_pos=0

print("Filtering contigs from " + f_in + " to " + f_out)

#get sample name
filename=os.path.basename(f_in)
sample=filename.split("_")[sam_pos]


def accept_contig(header):
    length = header.split("_")[3]
    cov = header.split("_")[5]
    if float(length) > length_threshold and float(cov) > coverage_threshold:
        return(True)
    else:
        return(False)

def get_head_info(header):
    head_info=header.split("_")
    info={"id": head_info[1], "length": head_info[3], "cov":head_info[5]}
    return(info)

fasta_sequences = SeqIO.parse(open(f_in),'fasta')
filt_conts=open(f_out, "a")
filt_stats=open(tab_out, "a")

count=1
for fasta in fasta_sequences:
    name, sequence = fasta.id, str(fasta.seq)
    head_inf=list(get_head_info(name).values())
    head_inf.insert(0, sample)

    print("\n--- Contig " + str(count) + ":")
    print("\t\n".join(head_inf))
    count+=1

    if accept_contig(name):
        head_inf.append("yes")
        entry="\t".join(head_inf)
        filt_stats.write(entry+"\n")
        SeqIO.write(fasta, filt_conts, "fasta")
    else:
        head_inf.append("no")
        entry="\t".join(head_inf)
        filt_stats.write(entry+"\n")

filt_stats.close()
filt_conts.close()
