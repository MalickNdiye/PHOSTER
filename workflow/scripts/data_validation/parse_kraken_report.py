import sys
import os
import pandas as pd

# This script takes in the kraken report and format it in a tabular way that can be more
# easily analyzed

indir=snakemake.params[0]
repdir=indir + "/Reports/"

# write header of summary file
outfile=open(snakemake.output[0], 'a')
outfile.write("sample\tpercentage_of_reads\treads_to_clade\treads_to_taxon\trank_code\ttaxID\ttaxon\tsample_type\n")

# loop through report files
for file in os.listdir(repdir):
    print("\nparsing")
    print(file)
    sample=file.split("_")[0] #take sample name

    if "P" in sample: #define if it comes from the virome or bacteriome
        sam_type="virome"
    else:
        sam_type="bacteriome"

    openfi=open(repdir + file, "r")
    for line in openfi:
        fields=line.split("\t")
        print("\nworking on")
        print(fields)
        fields[5]=fields[5].strip()

        if fields[5]=="Archaea":
            fields[5]="GB_VCs"

        fields.insert(0, sample)
        fields.append(sam_type)

        entry="\t".join(fields)
        entry=entry + "\n"
        print(entry)
        outfile.write(entry)

outfile.close()
