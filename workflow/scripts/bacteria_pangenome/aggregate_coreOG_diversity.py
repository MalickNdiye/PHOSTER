import os
import pandas as pd

infiles=list(snakemake.input[:])
n=len(infiles)
print("aggregateing core genome diversity for the following " + str(n) + " files: ")
print(infiles)

out_tab=snakeamke.output[0]

final_df=pd.DataFrame()
for f in infiles:
    df=pd.read_csv(f, delimiter="\t")
    print(df.head())
    
    final_df=pd.concat([final_df, df])

print("saving summmary table in: " + out_tab)
final_df.to_csv(out_tab, sep="\t", index=False)
