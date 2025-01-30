from zipfile import ZipFile
from zipfile import ZipInfo
import os
import pandas as pd


indirs=list(snakemake.input) # all the fastQC outputs, bith pre and post trimming
outfile=snakemake.output[0] # final summary file
outdir=os.path.dirname(outfile)

def process_zip(file, od=outdir, wanted="fastqc_data.txt"):
# This function opens the zip files outputted by fastQC. Then, it parse
# The info in the fastqc_data.txt file in order to obtain a summary of the QC
# The otput is a list of all the extracted info
    with ZipFile(file,'r') as inzipfile:
        for infile in (name for name in inzipfile.namelist() if name[-1] != '/'):
            if os.path.split(infile)[1] in wanted:
                with inzipfile.open(infile,'r') as txt:

                    for line in txt:
                        line=line.decode() # turn the line from bit to string

                        # Here I extract the basics stats of head file
                        if line.startswith("Filename"):
                            filename=line.split()[-1]
                            print("\nExtracting info from: " + filename)

                            info=filename.split("_")
                            sample=info[0]
                            id=info[-1]

                            if "paired" in filename:
                                trimming="post"
                                lane="NA"
                                direction=info[1]
                            else:
                                trimming="pre"
                                lane=info[1]
                                direction=info[2]




                            print("--- general stats saved")

                        #From here I start to extract the "pass" or "fail" for the following statistics
                        if line.startswith("Total Sequences"):
                            total_reads=line.split()[-1]
                            print("--- total reads saved")
                        if line.startswith("Sequences flagged as poor quality"):
                            poor_Q=line.split()[-1]
                            print("--- Poor quality saved")
                        if line.startswith("Sequence length"):
                            read_length=line.split()[-1]

                        if line.startswith(">>Per tile sequence quality"):
                            tile_Q=line.split()[-1]
                            print("--- Per tile quality saved")
                        if line.startswith(">>Per sequence quality scores"):
                            seq_Q=line.split()[-1]
                            print("--- Per sequence quality scores")
                        if line.startswith(">>Per base sequence content"):
                            base_Q=line.split()[-1]

                        if line.startswith(">>Per sequence GC content"):
                            GC_Q=line.split()[-1]

                        if line.startswith(">>Per base N content"):
                            N_Q=line.split()[-1]

                        if line.startswith(">>Sequence Length Distribution"):
                            length_distribution=line.split()[-1]

                        if line.startswith(">>Overrepresented sequences"):
                            overrep_seqs=line.split()[-1]

                        if line.startswith(">>Adapter Content"):
                            adapters=line.split()[-1]
                            print("--- Adapter content saved")

                summary_list=[sample, trimming, lane, direction, id, total_reads,
                 poor_Q, read_length, tile_Q, seq_Q, base_Q, GC_Q, N_Q,
                  length_distribution, overrep_seqs, adapters]

                txt.close()
                return(summary_list)




out=open(outfile, "w+")
headers=["sample", "trimming", "lane", "direction", "id", "total_reads",
 "poor_Qual_reads", "read_length", "per_tile_Q", "per_seq_Q", "per_base_Q",
  "GC_content", "N_content", "length_distribution", "overrep_seqs", "adapters"] # headers of the summary file

out.write("\t".join(headers))
out.write("\n")
out.close()

for dir in indirs:
    print("\nprocessing " + dir)
    zip_list=[dir + "/" + f for f in os.listdir(dir) if ".zip" in f]
    print("zip files: " + "\n\t".join(zip_list))
    for zip in zip_list:
        info=process_zip(zip)
        with open(outfile, "a") as out:
            out.write("\t".join(info))
            out.write("\n")

out.close()
