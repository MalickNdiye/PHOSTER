import sys
import pandas as pd
import os

cluster=sys.argv[1]
filt_db=sys.argv[2]
clust_info=sys.argv[3]

outfile=sys.argv[4]

df=pd.read_csv(clust_info, delimiter="\t")
df.head()

genomes=df.query("secondary_cluster==@cluster").genome.tolist()
genomes=[os.path.join(filt_db, g) for g in genomes]
print("genomes paths: ")
print(genomes)

with open(outfile, 'w') as f:
    for line in genomes:
        f.write(line + "\n ")

f.close()
