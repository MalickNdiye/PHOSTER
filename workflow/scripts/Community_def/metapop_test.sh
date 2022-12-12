#!/bin/bash

####--------------------------------------
##SLURM options
####--------------------------------------
#SBATCH --job-name metapop
#SBATCH --account pengel_beemicrophage
#SBATCH --cpus-per-task 25
#SBATCH --mem 50G
#SBATCH --time 7:00:00
#SBATCH --output /work/FAC/FBM/DMF/pengel/beemicrophage/mndiaye1/PHOSTER/workflow/logs/metapop/%x_%j.out
#SBATCH --error /work/FAC/FBM/DMF/pengel/beemicrophage/mndiaye1/PHOSTER/workflow/logs/metapop/%x_%j.err

cd /work/FAC/FBM/DMF/pengel/beemicrophage/mndiaye1/PHOSTER/workflow

bams=../results/mapping/mapdata/
ref_bact=../data/reference_assemblies/hb_bacteria/non_redundant/
cf=../results/mapping/library_count.tsv
threads=25
output=../results/metapop/test/
lib=/work/FAC/FBM/DMF/pengel/beemicrophage/mndiaye1/PHOSTER/workflow/.snakemake/conda/6ee4f263563c6d0db030078388037b81/lib/R/library/

mkdir -p ${output}

metapop --input_samples ${bams} --reference ${ref_bact} --norm ${cf} -l ${lib} -o ${output} --threads ${threads}
