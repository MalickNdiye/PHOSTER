import pysam
import os
import pandas as pd

bam_p=snakemake.input["sorted_bam"]

bamfile=pysam.AlignmentFile(bam_p, "rb")

#get sample name
basename = os.path.basename(bam_p)
sam_name= basename.split("_")[0]

# Initialize a dictionary to store the counts
contig_counts = {}
contig_mapqs = {}

# find multireads and store count and average mapQ for each scaffold
for read in bamfile.fetch():
    if read.is_secondary or read.is_supplementary: #skip secondary and supplementary aligment
        continue

    if read.has_tag("XS") and read.has_tag("AS") and read.get_tag("XS") == read.get_tag("AS"): # find multireads
        contig_name = bamfile.get_reference_name(read.reference_id)
        mapq = read.mapping_quality

        contig_counts[contig_name] = contig_counts.get(contig_name, 0) + 1
        contig_mapqs[contig_name] = contig_mapqs.get(contig_name, 0) + mapq



bamfile.close()

for contig in contig_counts:
    if contig in contig_mapqs:
        contig_mapqs[contig] = contig_mapqs[contig] / contig_counts[contig]


df = pd.DataFrame({"scaffold": list(contig_counts.keys()), "multireads": list(contig_counts.values()), "average_mapq_multireads": list(contig_mapqs.values())})
df["sample"]=[sam_name]*df.shape[0]


df.to_csv(snakemake.output["multireads_count"], sep="\t", index=False)
