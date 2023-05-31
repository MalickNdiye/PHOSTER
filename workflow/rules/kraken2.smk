################################################################################
######################### Kraken2 DB building ##################################
################################################################################
# Build a costum kraken database to quickly analyse the community of our
# sequencing data. this rule can take more than 5 hrs to run. Honeybee phages
# are classified as "archea" to save me time and avoid to write the whole
# Phylogeny from scratch. Then the classification will be changed when the
# kraken reports ae parsed.
################################################################################
################################################################################

rule build_Krakern2:
    input:
        GB_VC= "../data/reference_assemblies/VCs_kraken/GB_ViralClusters_kraken.fna",
        Amel= "../data/reference_assemblies/A_mellifera/GCF_003254395.2_Amel_HAv3.1_genomic.fna"
    output:
        directory("resources/databases/220131_costum_kraken2db")
    threads: 42
    conda:
        "envs/Kraken2.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "06:00:00",
        disk_mb= 1000000
    log:
        "logs/data_validation/kraken2/database_kraken2.log"
    shell:
        "kraken2-build --download-taxonomy --db {output}; "
        "kraken2-build --download-library bacteria --db {output}; "
        "kraken2-build --download-library viral --db {output}; "
        "kraken2-build --download-library human --db {output}; "
        "kraken2-build --add-to-library {input.GB_VC}  --db {output}; "
        "kraken2-build --add-to-library {input.Amel} --db {output}"