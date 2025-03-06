PHOSTER
==============
**Author:** *Malick ndiaye*

# INTRODUCTION

This is the bioinformatics pipeline of the project PHOSTER. The pre-print is now available!

# PROJECT DESCRIPTION

Bacteriophages (phages) play a crucial role in shaping bacterial communities, yet the relationship between phage and bacterial diversity remains poorly understood, particularly in animal-associated microbiomes. Here, we reconstructed the phage-bacteria interaction network in the gut microbiota of 49 honey bees by leveraging CRISPR spacer matches and genome homology. The resulting interaction network displayed a highly modular structure with nested phage-bacteria interaction within each module. Viral and bacterial alpha and beta diversity were correlated, particularly at the bacterial strain level and when explicitly considering the interaction network in diversity analysis. Overall, our results suggest that the most relevant approach to study phage-bacteria diversity patterns should rely on strain-level resolution and the explicit use of the interaction network. This offer an explanation as to why previous studies have failed to detect phage-bacteria diversity correlations. Finally, we call for further studies building up on these correlation patterns to probe the underlying mechanisms by integrating both bottom-up and top-down regulatory mechanisms in microbiome assembly.

# DIRECTORY STRUCTURE

.
├── README.md --> This file
├── config --> Directory containing config file for snakemake
└── workflow --> COntains Snakefile, conda environment and scripts to run the pipeline

# RAW DATA

From 49 bees collected iin Lausanne(Switzerland) we extracted DNA form the bacterial and viral fraction of the gut microbiota. Libraries were prepared with the Illumina Nextera Flex library kit (Illumina) with unique dual indices (UDI) and sequenced on a Illumina NovaSeq 6000 instrument (PE150) at the Genomic Technologies Facility (GTF) of the University of Lausanne. Raw reads data can be found on SRA, BioProject (PRJNA1232403)[]

# PIPELINE

This snakemake pipeline is divided in several step described below:

## 1. Data Validation

This section is divided in several step:

1. Community profiling using Kraken2
2. Raw reads QC and Trimming
3. Host Filtering

In (1) raw reads are profiled using a kraken database containing all the bacterial, viral and human kmers. Moreover, I added the already charachterized honeybee phages as "archeas" in order to profile those as well. This step is used to have a quick overview of the domain-level composition of our samples. In (2), raw reads are quality checked and trimmed. trimmed reads are then filtered for honeybee genome and human reads in (3), leaving only the bacterial and viral reads (plus contaminants from pollen and other microorganisms). 

At the end of this section, reads are ready for in-depth analysis.

##  2. Bacterial MAGS

In this part we assembled MAGs from the bacterial fraction

1. The host-filtered reads are assembled using metaspades
2. backmapping where the read of each sample are mapped against the assembly of every other sample
3. scaffold are binned into bMAGS in functon of their coverage across samples using metabat2
4. bMAGs are QCed using CheckM, filtered in function of completeness (>75%) and contamination (<10%), and taxonomically classified using GTDB-TK
5. Cluster bMAGsat 95% ANI using dRep to create a reference database

At the end of this section we will a compendium of bMAGs that are taxonomically classified, and divided into bOTUs based on ANI

## 3. Phage Identification

In this section, we identify phages in both the bacterial and virl fraction using virsorter2, VIBRANT and ViralVerify. Then the output of these tools is aggragated. a normalized score is given to all contigs based on how many tools identified it as a phage and the confidence in the prediction. Finally, all contigs that pass a given threshols (more or less equivalent to being identified with low confidence by 3 tools) are retained.

## 4. Polish Viral Contigs

In this section viral contigs are trimmed to remove bacterial contaminations (in the case of prophages). Moreover, the quality of the contigs is determined using checkV. Finally, contigs are filtered again in function of lentght (<10kb) and contamination.

## 5. Dereplicate viral contigs

In this section use dRep to dereplicate the viral contigs at 95% ANI and 85% AF. There are two types of dereplication for 95% ANI and 85 %AF:
- Dereplication using average-linkage clustering. This will yield a dereplicated vMAGs database to map reads.
- Dereplication using single -inkage clustering. This will yield the vOTUS (as suggested by the [dRep developer](https://drep.readthedocs.io/en/latest/choosing_parameters.html))

I do this double dereplication because I noticed that the average-linkage clustering often results in contigs with >98% ANI in different clusters. This is a issue for mapping because reads willl be assigned to multiple vOTUs. On the other hand, The single-linkage create clusters contain genomes that are too divergent >95% ANI and the representative genomes may not be detected in some samples although members of the same vOTU would be detected. I think this issue stems from the high variability of viruses that creates a continuum of divergence in ANI that goes beyond the 95% ANI 85% AF standard threshold. So, to me the single link makes more sense to represent a cohesive genetic unit (vOTU) but I use the average-linkage to avoid mapping issues caused by the divergence within a vOTU (mlti-mapped reads will end up counting for the same vOTU in this case). 

An all vs all alignment is also performed in this section using fastANI.

## 6. Host Assignation

To assign hosts to every viral contigs we will use the CRISPR spacers and genome homology:

1. find CRISPR spacers,DR and cas genes in all MAGs and refernce bacterial genome. filter for good quality spacers
2. Create databases of CRISPR spacers and map them to the redundant viral contigs 
3. Use Fastani to determine genome Homology between phages and bacteria
4. Parse the results

During the results parsing, host X is targeted by phage Y if:
- A CRSIPR from host X maps with maximum 2 mismathces to phage Y (gaps are counted in this case)
- Phage Y has >90% ANI and >50% AF with host X (if AF >80%, it is considered a prophage). usually people use blast to do this, but [Johansen et al. (2022)](https://www.nature.com/articles/s41467-022-28581-5) show that FastANI works just as good

All the host assigned to a given viral contig are assigned to the entire vOTU

## 7. PBIN analyisis

In this section, host-phage linkage info are used to build a phage bacteria interaction network. Then the lpBRIM alogorithm is used to identify modules. Finally nestdness is computed for each modules. These operations are performed 2 times: one for all bacterial genomes (bMAGs and Isolates), and once using only isolates genomes. 

## 8. Bacteria Phylogeny and viral Genetic Relationships

This section is divided in two: Bacteria and Phages.

- Bacteria: all bacteria genomes are annotated. Protein seqcuenced are fed to Orthofinder to find shared Orthogroups. A trimmed alignment of the shared OGs is used to build a phylogeny.

- Phages: All representative vOTU genomes are annotated. Then protein sequences are fed to vContact2 to create a network of protein sharedness. finally, info on the interaction modules from the PBIN are added to the network.

## 9. Map to refrerence database

This section, we map the bacterial fraction's reads to the dereplicated bacteria database and the viral fraction's reads to the dereplicated viral database. 

## 10. Community Analysis

In this section I run Instrain profile on all samples and ummarize the results. Moreover, I count SNVs for each bacterial genome across samples. Finally, I run InStrain compare on all the bacterial genomes.

# DOWNSTREAM ANALYSIS
scripts for downstram analysis and figure generation can be found on [Zenodo]()