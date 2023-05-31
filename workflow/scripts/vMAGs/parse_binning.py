import pandas as pd
import numpy as np
import os
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re
import copy



parser = argparse.ArgumentParser(
    description="""This script takes vrhyme bin and virus identification stats and filter the assembly for binned or high quality viral contigs. also
    produces nice summary tables""", 
    formatter_class=argparse.RawDescriptionHelpFormatter,
    add_help=False,
)


parser.add_argument("-i", "--infasta", help="Paths to input FASTA files", nargs="+")
parser.add_argument("-b", "--bins", help="path to vRhyme binning directories", nargs="+")
parser.add_argument("-v", "--vi", help="Path to viral identification stats")
parser.add_argument("-f", "--outfilt", help="Path to output filtered assembly")
parser.add_argument("-c", "--outconcat", help="Path to output concatenated assembly")
parser.add_argument("-r", "--outbin", help="Path to output binning data")
parser.add_argument("-t", "--outtab", help="Path to output filtering summary table")
args = parser.parse_args()

print("\n################################################################################")
print("Arguments provided:")
print("Input FASTA file: \n{}".format(args.infasta))
print("\nvRhyme binning directory: \n{}".format(args.bins))
print("\nviral identification stats: {}".format(args.vi))
print("\nOutput filtered assembly: {}".format(args.outfilt))
print("\nOutput concatenated assembly: {}".format(args.outconcat))
print("\nOutput binning data: {}".format(args.outbin))
print("\nOutput filtering summary table: {}".format(args.outtab))
print("################################################################################")



# assemblies
assemblies=args.infasta

# binning directoiries
bins=args.bins

# metadata quality
mtadata_quality=pd.read_csv(args.vi ,sep="\t", low_memory=False)
mtadata_quality=mtadata_quality.round({"quality":0})

# get list of contigs name with quality >0 and length > 2999 according to metadata
quality_contigs_all=mtadata_quality.loc[(mtadata_quality['quality']>0) & (mtadata_quality['len']>2999) ]

# concatenate binning files
binning_files=[os.listdir(x) for x in bins]
binning_files=[j for i in binning_files for j in i]
best_bins=[x for x in binning_files if "membership.tsv" in x]

############################################################################################################3
print("\n\tConcatenating binning files")


def get_binning_data(dir):
    # get semple name from dir name
    sample_name=dir.split("/")[-1]
    if sample_name=="":
        sample_name=dir.split("/")[-2]

    # read file ending with membership.tsv in dir
    print("\n\nReading binning data for sample {}".format(sample_name))
    print(os.listdir(dir))

    mem_list=[x for x in os.listdir(dir) if "membership.tsv" in x]
    
    if len(mem_list)==0:
        print("No binning data for sample {}".format(sample_name))
        #return empt dataframe with columns sample, contig, bin, members, proteins, redundancy 
        return(pd.DataFrame(columns=["sample","contig","bin","members","proteins","redundancy"]))


    mem_file=mem_list[0]
    mem_path=os.path.join(dir,mem_file)
    membership=pd.read_csv(mem_path,sep="\t")
    membership=membership.assign(sample=sample_name)
    
    

    # move sample to first column and rename scaffodl to contig
    membership=membership[["sample","scaffold", "bin"]]
    membership=membership.rename(columns={"scaffold":"contig"})

    # read file ending with summary.tsv in dir
    summary=pd.read_csv(os.path.join(dir,[x for x in os.listdir(dir) if "summary.tsv" in x][0]),sep="\t")

    # merge membership and summary
    joined=pd.merge(membership,summary,on="bin")

    return(joined)

# loop through binning directories and concatenate binning data in one dataframe
binning_data=[get_binning_data(x) for x in bins]
binning_data=pd.concat(binning_data)

###############################################################################################################

# write a function that takes a an object like list(SeqIO.parse(assembly,"fasta")) and a filename and
#  renames all the contigs as the filename + a unique number that's always the same for the same contig
def rename_contigs(assembly, filename, reference):
    # get sample name from assembly file name
    sample_name=filename.split("/")[-1].split("_")[0]
    print("Renaming contigs for sample {}".format(sample_name))

    assembly_obj2=copy.deepcopy(assembly)
    old_names=[x.description for x in assembly_obj2]

    # rename contigs
    for i in range(len(assembly_obj2)):
        sub_name="_".join(assembly_obj2[i].description.split("_")[0:6])

        # check if new name is already in the assembly
        if sum([1 for x in reference if sub_name in x.description])>1:
            # if it is, add a number to the end of the name that is a sum of all the numbers present in the original name
            new_name="{}_vcontig_{}".format(sample_name,"_".join(assembly_obj2[i].description.split("_")[0:6]))
            new_name=new_name+"_{}".format(sum([int(s) for s in re.findall(r'\d+', assembly_obj2[i].description)])+len(assembly_obj2[i].description))
            assembly_obj2[i].id=new_name
            assembly_obj2[i].description=new_name

        elif "bin" in assembly_obj2[i].id:
            new_name=assembly_obj2[i].id
            assembly_obj2[i].id=new_name
            assembly_obj2[i].description=new_name

        else:
            new_name="{}_vcontig_{}".format(sample_name,"_".join(assembly_obj2[i].description.split("_")[0:6]))
            assembly_obj2[i].id=new_name
            assembly_obj2[i].description=new_name

    return(assembly_obj2)




# function that concatenate conigs in same bin
def concat_contigs_bin(sample_bins, sample_binned_data, sample_name, assembly_fasta, outfile):
    for bin in sample_bins:
                
                # get contigs in bin
                bin_contigs=sample_binned_data[sample_binned_data["bin"]==bin]["contig"].tolist()
                print("\n\tConcatenating contigs in bin {} the contigs {} for sample {}".format(bin, " ".join(bin_contigs), sample_name))
        
                # get contigs in bin from assembly
                bin_contigs_assembly=[x for x in assembly_fasta if x.description in bin_contigs]
        
                # concatenate contigs in bin
                bin_contigs_assembly_concatenated=SeqRecord(Seq("".join([str(x.seq) for x in bin_contigs_assembly])),id="{}_bin_{}".format(sample_name,bin),description="")
    
                # write concatenated contigs to file
                with open(outfile,"a+") as f:
                    SeqIO.write(bin_contigs_assembly_concatenated,f,"fasta")


# make a dataframe
# the first column should be the contig name of all contigs in the assembly
# the second column is called "bin" and is the bin number if the contig is binned and the contig name if the contig is not binned
# the third column is called "quality" and is the quality of the contig according to metadata 
# the column is filtered and is yes if the contig is binned or has quality >0 according to metadata and no otherwise

def get_assembly_df(assembly_fasta_org, mtadata_quality, binned_contigs, retained_contigs, binning_data, sample_name):
    # make dataframe
    assembly_df=pd.DataFrame({"old_name":[x.description for x in assembly_fasta_org],
                                "contig":[x.description for x in rename_contigs(assembly_fasta_org, args.outfilt, assembly_fasta_org)],
                                "quality":[mtadata_quality[mtadata_quality["contig"]==x.description]["quality"].tolist()[0] if x.description in mtadata_quality["contig"].tolist() else "unclassified" for x in assembly_fasta_org],
                                "filtered":["yes" if x.description in retained_contigs else "no" for x in assembly_fasta_org]})
    
    # add bin column
    assembly_df["bin"]=[sample_name + "_bin_" + str(binning_data[binning_data["contig"]==x.description]["bin"].tolist()[0]) if x.description in binned_contigs else 
    y.description for x,y in zip(assembly_fasta_org, rename_contigs(assembly_fasta_org, args.outfilt, assembly_fasta_org))]

    # add sample name column
    assembly_df["sample"]=sample_name

    # add length column taken from metadata file for the contig
    assembly_df["length"]=[mtadata_quality[mtadata_quality["contig"]==x]["len"].tolist()[0] if x in mtadata_quality["contig"].tolist() else len(str(y.seq)) for x,y in zip(assembly_df["contig"],assembly_fasta_org)]

    return(assembly_df)

#######################################################################################################

def filter_assembly(assembly):
    binned_contigs=binning_data["contig"].unique().tolist()
    print("binned contigs: {}".format(binned_contigs))

    # get sample name from assembly file name
    sample_name=assembly.split("/")[-1].split("_")[0]
    print("\n\nFiltering assembly for sample {}".format(sample_name))

    ## get list of dereplicated contigs from bin/vRhyme_dereplication/vRhyme_derep_longest_<sample_name>_viral_contigs_trimmed.fa
    bin_dir=[x for x in bins if sample_name in x][0]
    dereplicated_contigs=list(SeqIO.parse(os.path.join(bin_dir,"vRhyme_dereplication","vRhyme_derep_longest_{}_viral_contigs_trimmed.fa".format(sample_name)),"fasta"))
    dereplicated_contigs=[x.description for x in dereplicated_contigs if len(x.seq)>2999 or x.description in binned_contigs]

    # read assembly fasta file
    assembly_fasta_org=list(SeqIO.parse(assembly,"fasta"))

    # get list of contigs name with quality >0 and length > 2999 according to metadata
    quality_contigs=quality_contigs_all.loc[(quality_contigs_all['sample']==sample_name)]["contig"].tolist()

    # add contigs trimmed by checkv
    quality_contigs=quality_contigs + [x.description for x in assembly_fasta_org if x.description not in mtadata_quality["contig"].tolist()]

    # filter assembly for contigs that are either binned
    # or have quality >0
    assembly_fasta=[x for x in assembly_fasta_org if x.description in binned_contigs or x.description in quality_contigs]
    # or have len > 2999 according to metadata 
    assembly_fasta=[x for x in assembly_fasta if len(x.seq)>2999 or x.description in binned_contigs]
    assembly_fasta=[x for x in assembly_fasta if x.description in dereplicated_contigs]
    retained_contigs=[x.description for x in assembly_fasta]

    # write filtered assembly to file
    print("---Writing filtered assembly for sample {}".format(sample_name))
    os.makedirs(os.path.dirname(args.outfilt), exist_ok=True)
    SeqIO.write(rename_contigs(assembly_fasta,args.outfilt, assembly_fasta_org), args.outfilt, "fasta") 
    

    # concatenate contigs in the same bin
    # contigs in the same bin should be concatenated and written as one contig named sample_bin_binnumber
    # contigs that are not binned but have quality >0 and len >2999 according to metadata are written with their original name

    # get binned data for sample
    sample_binned_data=binning_data[binning_data["sample"]==sample_name]

    # get bins for sample
    sample_bins=sample_binned_data["bin"].unique().tolist()

    # loop through bins
    os.makedirs(os.path.dirname(args.outconcat), exist_ok=True)
    concat_contigs_bin(sample_bins, sample_binned_data, sample_name, assembly_fasta, args.outconcat)
    # write contigs that are not binned but have quality >0 according to metadata
    # contigs that are not binned but have quality >0 according to metadata are written with their original name

    # get contigs that are not binned but have quality >0 and length > 2999 according to metadata
    quality_contigs_assembly=[x for x in assembly_fasta if x.description not in binned_contigs]

    # append conntigs that are not binned but have quality >0 according to metadata to file assembly.replace("trimmed","assembly_filtered_concatenated")
    print("\nWriting concat assembly for sample {}".format(sample_name))
    with open(args.outconcat,"a") as f:
        SeqIO.write(rename_contigs(quality_contigs_assembly, args.outconcat, assembly_fasta_org),f,"fasta")


    # make a dataframe
    # the first column should be the contig name of all contigs in the assembly
    # the second column is called "bin" and is the bin number if the contig is binned and the contig name if the contig is not binned
    # the third column is called "quality" and is the quality of the contig according to metadata 
    # the column is filtered and is yes if the contig is binned or has quality >0 according to metadata and no otherwise

    # make dataframe
    assembly_df=get_assembly_df(assembly_fasta_org, mtadata_quality, binned_contigs, retained_contigs, binning_data, sample_name)

    return(assembly_df)

# loop through assemblies and filter them
filtered_assemblies=[filter_assembly(x) for x in assemblies]

# concatenate filtered assemblies
filtered_assemblies=pd.concat(filtered_assemblies)

# write filtered assemblies to file
print("Writing filtered assemblies data to file")
os.makedirs(os.path.dirname(args.outtab), exist_ok=True) 
filtered_assemblies.to_csv(args.outtab,sep="\t",index=False)

# write binning data to file
print("Writing binning data to file")
binning_data.to_csv(args.outbin,sep="\t",index=False)