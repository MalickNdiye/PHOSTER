####################################################################################
################################## viralMAGS #######################################
####################################################################################
# This part of the pipleine takes takes host contigs identified as viral and
# maps bins them into vMAGS using vRhyme.
# The procedure is similar to the one to create MAGs, just using vRhyme as binning tool.
################################################################################
################################################################################

# This rule builds a mapping index for each sample
rule build_backmapping_index_viruses:
    input:
        ref="../results/assembly/viral_contigs/single_sample/{asmbl}_viral_contigs_trimmed.fasta"
    output:
        dir=directory("../results/assembly/viral_contigs/indexes/{asmbl}_viral_contigs/{asmbl}_viral_contigs_index")
    conda:
        "envs/map_env.yaml"
    threads: 15
    params:
        basename="{asmbl}_viral_contigs_trimmed"
    log:
        "logs/vMAGs/backmapping/indexing/{asmbl}_build_bowtie_index.log"
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 20000,
        runtime= "01:00:00"
    shell:
        "mkdir -p {output.dir}; "
        "bowtie2-build {input.ref} {output.dir}/{params.basename} --threads {threads}"

# This rule uses bowtie2 to map reads to the viral contigs
rule backmapping_viral: 
    input:
        assembly = "../results/assembly/viral_contigs/indexes/{asmbl}_viral_contigs/{asmbl}_viral_contigs_index", # we bin only bacteria
        R1 = "../data/host_filtered_reads/{sample}_R1_HF.fastq.gz",
        R2 = "../data/host_filtered_reads/{sample}_R2_HF.fastq.gz"
    output:
        sam=temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs_trimmed.sam")
    resources:
        account="pengel_beemicrophage",
        runtime="10:00:00",
        mem_mb = 10000
    params:
        basename="{asmbl}_viral_contigs_trimmed"
    threads: 15
    log:
        "../scratch_link/logs/{sample}_to_{asmbl}_back.log"
    conda: "envs/map_env.yaml"
    shell:
        "bowtie2 -x {input.assembly}/{params.basename} -1 {input.R1} -2 {input.R2} -S {output.sam} --threads {threads}"

# This rule takes the bam file and creates a depth file
rule backmapping_depths_viral:
    input:
        sam= "../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs_trimmed.sam"
    output:
        bam= temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs_trimmed.bam"),
        depth= temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs_trimmed.depth")
    resources:
        account="pengel_beemicrophage",
        runtime="01:30:00",
        mem_mb = 10000
    params:
        tmp="../scratch_link/tmp"
    threads: 15
    conda: "envs/sam_env.yaml"
    log:
        "../scratch_link/logs/{sample}_to_{asmbl}_back_depth.log"
    shell:
        "samtools view -bh {input.sam} -@ {threads} | samtools sort -T {params.tmp} - > {output.bam}; "
        "export OMP_NUM_THREADS={threads}; "
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam}"

# This rule merges the depth files
rule merge_depths_viral:
    input:
        depth_files = expand("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{{asmbl}}_viral_contigs_trimmed.depth", sample=config["samples"])
    output:
        depth_file_merged = "../results/vMAG_binning/backmapping/merged_depths/{asmbl}_global_viral_depth.txt"
    resources:
        account="pengel_beemicrophage",
        runtime="1:00:00",
        mem_mb = 10000
    threads: 1
    conda: "envs/mags_env.yaml"
    log:
        "logs/vMAGs/backmapping/depth/{asmbl}_merge_depth.log"
    shell:
        "scripts/MAGs/merge_depths.pl {input.depth_files} > {output.depth_file_merged}"

# this rules convert depth files from METAbat2 format to vRhyme format
rule covert_depths:
    input:
        depth_merged = "../results/vMAG_binning/backmapping/merged_depths/{asmbl}_global_viral_depth.txt"
    output:
        depth_file_merged_vrhyme = "../results/vMAG_binning/backmapping/merged_depths/{asmbl}_global_depth_vrhyme.txt"
    resources:
        account="pengel_beemicrophage",
        runtime="1:00:00",
        mem_mb = 10000
    threads: 1
    params:
        conda="resources/conda_envs/vrhyme"
    log:
        "logs/vMAGs/backmapping/depth/convert/{asmbl}_covert_depth.log"
    shell:
        """
        bash -c '. $HOME/.bashrc
            conda activate {params.conda}
            coverage_table_convert.py -i {input} -o {output}'
        """

# this rule runs vRhyme
rule run_vRhyme:
    input:
        assembly="../results/assembly/viral_contigs/single_sample/{asmbl}_viral_contigs_trimmed.fasta",
        depth_file = "../results/vMAG_binning/backmapping/merged_depths/{asmbl}_global_depth_vrhyme.txt"
    output:
        out= directory("../results/vMAG_binning/vRhyme_bins/{asmbl}/")
    resources:
        account="pengel_beemicrophage",
        runtime="0:30:00",
        mem_mb = 100000
    threads: 15
    params:
        conda="resources/conda_envs/vrhyme",
        tmp="../scratch_link/vMAG_binning/vRhyme_bins/{asmbl}"
    log:
        "logs/vMAGs/backmapping/binning/{asmbl}_virus_binning.log"
    benchmark:
        "logs/vMAGs/backmapping/binning/{asmbl}_virus_binning.benchmark"
    shell:
        """
        bash -c '. $HOME/.bashrc
            conda activate {params.conda}
            vRhyme -i {input.assembly} -c {input.depth_file} -o {params.tmp} --prefix {wildcards.asmbl}_vMAG_ --iter 11 --method longest --derep_id 1  -t {threads}'
            mkdir -p {output}
            mv {params.tmp}/* {output}/
            rm -rf {output}/vRhyme_alternate_bins
        """

# this rule parses the vRhyme output
rule parse_vRhyme:
    input:
        assembly="../results/assembly/viral_contigs/single_sample/{sample}_viral_contigs_trimmed.fasta",
        bins="../results/vMAG_binning/vRhyme_bins/{sample}",
        vi="../results/viral_identification/ViralIdentification_scores.tsv"
    output:
        filtered_cont="../results/assembly/viral_contigs/filtered/{sample}_viral_contigs_filtered.fasta",
        concat_cont="../results/assembly/viral_contigs/concat/{sample}_viral_contigs_concat.fasta",
        binning_data=temp("../results/vMAG_binning/summary_table/{sample}_binning_data.tsv"),
        filtering_data=temp("../results/vMAG_binning/summary_table/{sample}_filtering_data.tsv")
    resources:
        account="pengel_beemicrophage",
        runtime="1:00:00",
        mem_mb = 100000
    threads: 1
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/vMAGs/vRhyme/{sample}_parsing.log"
    benchmark:
        "logs/vMAGs/vRhyme/{sample}_parsing.benchmark"
    shell:
        "python scripts/vMAGs/parse_binning.py -i {input.assembly} -b {input.bins} -v {input.vi} -f {output.filtered_cont} -c {output.concat_cont} -r {output.binning_data} -t {output.filtering_data}"

# this rule aggregates the vRhyme output
rule aggregate_parsed_vrhyme:
    input:
        binning_data=expand("../results/vMAG_binning/summary_table/{sample}_binning_data.tsv", sample=config["samples"]),
        filteing_data=expand("../results/vMAG_binning/summary_table/{sample}_filtering_data.tsv", sample=config["samples"])
    output:
        all_bin="../results/vMAG_binning/summary_table/aggregate_binning_data.tsv",
        all_filt="../results/vMAG_binning/summary_table/aggregate_filtering_data.tsv"
    resources:
        account="pengel_beemicrophage",
        runtime="1:00:00",
        mem_mb = 10000
    threads: 1
    log:
        "logs/vMAGs/vRhyme/aggregate_parsing.log"
    benchmark:
        "logs/vMAGs/vRhyme/aggregate_parsing.benchmark"
    run:
        import pandas as pd

        binnings = [pd.read_csv(i, sep="\t") for i in input.binning_data]
        binnings = pd.concat(binnings)

        filterings = [pd.read_csv(i, sep="\t") for i in input.filteing_data]
        filterings = pd.concat(filterings)

        binnings.to_csv(output.all_bin, sep="\t", index=False)
        filterings.to_csv(output.all_filt, sep="\t", index=False)