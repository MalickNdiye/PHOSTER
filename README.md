# PHOSTER
Pipeline of the anlysis of the honeybee gut bacteriome and virome.

# Data Validation

This section is divided in several step:

  (1) Community profiling using Kraken2
  (2) Raw reads QC and Trimming
  (3) Host Filtering
  (4) Community Profiling using mOTUs
 
In (1) raw reads are profiled using a kraken database containing all the bacterial, viral and human kmers. Moreover, I added the already charachterized honeybee phages as "archeas" in order to profile those as well. This step is used to have a quick overview of the domain-level composition of our samples. In (2), raw reads are quality checked and trimmed. trimmed reads are then filtered for honeybee genome and human reads, leaving only the bacterial and viral reads (plus contaminants from pollen and other microorganis). Finally host-filtered reads are profiled using mOTUs, which is has more taxonomic resolution than Kraken2, allowing us to obtain the genus-level community composition.
