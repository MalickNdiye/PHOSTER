rule concat_viral_assemblies:
    input:
        refs="../data/reference_assemblies/Phages_GB/HBvirDBv1_DB/Single_Samples"
    output:
        "../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna"
    conda:
        "envs/vBinning_env.yaml"
    threads: 1
    log:
        "logs/vMAGs/backmapping/concat_genomes/HBvirDBv1_build_bowtie_index.log"
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 20000,
        runtime= "01:00:00"
    shell:
        "python scripts/vMAGs/concatenate.py {output} {input.ref} --nozip"


rule build_backmapping_index:
    input:
        ref="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna"
    output:
        dir=directory("../results/mapping/HBvirDBv1_DB_index")
    conda:
        "envs/map_env.yaml"
    threads: 15
    params:
        basename="HBvirDBv1_DB",
    log:
        "logs/vMAGs/backmapping/indexing/HBvirDBv1_build_bowtie_index.log"
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 20000,
        runtime= "01:00:00"
    shell:
        "mkdir -p {output.dir}; "
        "bowtie2-build {input.ref} {output.dir}/{params.basename} --threads {threads}"



rule backmapping_viral_db: # risk of running out of buffer memory, run fewer jobs at the time (with --jobs 50 works; maybe one could increase a bit)
    input:
        assembly = "../results/mapping//HBvirDBv1_DB_index", # we bin only bacteria
        R1 = "../data/host_filtered_reads/{sample}_R1_HF.fastq.gz",
        R2 = "../data/host_filtered_reads/{sample}_R2_HF.fastq.gz"
    output:
        sam=temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1.sam")
    resources:
        account="pengel_beemicrophage",
        runtime="10:00:00",
        mem_mb = 10000
    params:
        basename="HBvirDBv1_DB"
    threads: 15
    conda: "envs/map_env.yaml"
    log:
        "logs/vMAGs/backmapping/map/{sample}_mapped_to_HBvirDBv1.log"
    benchmark: 
        "logs/MAGs/backmapping/map/{sample}_mapped_to_HBvirDBv1.benchmark"
    shell:
        "bowtie2 -x {input.assembly}/{params.basename} -1 {input.R1} -2 {input.R2} -S {output.sam} --threads {threads}"


rule sort_virome_backmapping_bam:
    input:
        sam="../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1.sam"
    output:
        bam= temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1.bam"),
        sorted_bam="../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1_sorted.bam",
        bai="../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1_sorted.bam.bai"
    log:
        "logs/vMAGs/backmapping/soringbams/sort_bam_{sample}.log"
    conda:
        "envs/map_env.yaml"
    threads: 5
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "01:30:00"
    shell:
        "samtools view -bh {input.sam} | samtools sort -T {params.tmp} - > {output.bam};"
        "samtools sort {input.bam} -@ {threads} -o {output.sorted_bam}; "
        "samtools index {output.sorted_bam} {output.bai} -@ {threads}"
    

rule run_vamb:
    input:
        ref="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna",
        sorted_bam=expand("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1_sorted.bam", sample=config["samples"])
    output:
        out=temp("../results/MAG_binning/vBins/HBvirDBv1/"),
    log:
        "logs/vMAGs/binning.log"
    conda:
        "envs/vBinning_env.yaml"
    threads: 25
    params:
        min_contig_l=3000
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "4:00:00"
    shell:
        "vamb --outdir {output} --fasta {input.ref} --bamfiles {input.soret_bams} -m {params.min_contig_l} -p {threads} -o C"