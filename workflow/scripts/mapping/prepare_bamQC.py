import os


bam_list=list(snakemake.input["bam"])
outfile=snakemake.output[0]

f=open(outfile, "w")
for i in bam_list:
    sample=i.split("/")[-1]
    sample=sample.split("_")[0]
    
    rel_p=os.path.relpath(i, outfile)
    rel_p=rel_p.replace('../', '', 1)
    
    print("saving location of bam for sample " + sample)
    print("\tlocation" + rel_p + "\n")
    f.write(sample + "\t" + rel_p + "\n")
    
f.close()