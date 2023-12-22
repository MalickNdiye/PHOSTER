# This is the old version of the viral pipeline 

rule parse_VI:
    input:
        vib_score = expand("../results/viral_identification/vibrant/{sample}_vibrant/{sample}_vibrant_score.tsv", sample=config["samples"]),
        vib_prophages=expand("../results/viral_identification/vibrant/{sample}_vibrant/{sample}_vibrant_prophages.tsv", sample=config["samples"]),
        vv_score = expand("../results/viral_identification/viralverify/{sample}_viralverify/{sample}_viralverify_score.csv", sample=config["samples"]),
        vs_score = expand("../results/viral_identification/virsorter/{sample}_virsorter/{sample}_virsorter_score.tsv", sample=config["samples"]),
        vs_prophages= expand("../results/viral_identification/virsorter/{sample}_virsorter/{sample}_virsorter_boundaries.tsv", sample=config["samples"])
    output:
        all_scores ="../results/viral_identification/ViralIdentification_scores.tsv",
        vib_scores="../results/viral_identification/vibrant/vibrant_all_scores.tsv",
        vv_scores ="../results/viral_identification/viralverify/viralverify_all_scores.tsv",
        vs_scores ="../results/viral_identification/virsorter/virsorter_all_scores.tsv",
        tab="../results/viral_identification/identification_tools_tab.tsv"
    log:
        "logs/viral_identification/parse_VI.log"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 10000,
        runtime = "00:10:00"
    threads:1
    conda:
        "envs/base_R_env.yaml"
    script:
        "scripts/viral_identification/parse_viral_identification.R"

rule extract_viral_contigs:
    input:
        all_scores ="../results/viral_identification/ViralIdentification_scores.tsv",
        concat_assembly = "../results/assembly/concat_assembly/{sample}_concat_assembly.fasta"
    output:
        "../results/assembly/viral_contigs/single_sample/{sample}_viral_contigs.fasta"
    log:
        "logs/viral_identification/extract_viruses/{sample}_viral_contigs.log"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 10000,
        runtime = "00:20:00"
    threads:1
    conda:
        "envs/mOTUpan.yaml"
    shell:
        "python scripts/viral_identification/extract_viruses_from_assembly.py {input.concat_assembly} {input.all_scores} {output}"


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
        ref="../results/assembly/viral_contigs/single_sample/{asmbl}_viral_contigs.fasta"
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
        sam=temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs.sam")
    resources:
        account="pengel_beemicrophage",
        runtime="10:00:00",
        mem_mb = 10000
    params:
        basename="{asmbl}_viral_contigs_trimmed"
    threads: 15
    log:
        "../scratch_link/logs/{sample}_to_{asmbl}_back.log"
    conda: 
        "envs/map_env.yaml"
    shell:
        "bowtie2 -x {input.assembly}/{params.basename} -1 {input.R1} -2 {input.R2} -S {output.sam} --threads {threads}"

# This rule takes the bam file and creates a depth file
rule backmapping_depths_viral:
    input:
        sam= "../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs.sam"
    output:
        bam= temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs.bam"),
        depth= temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{asmbl}_viral_contigs.depth")
    resources:
        account="pengel_beemicrophage",
        runtime="01:30:00",
        mem_mb = 10000
    params:
        tmp="../scratch_link/tmp"
    threads: 15
    conda: 
        "envs/sam_env.yaml"
    log:
        "../scratch_link/logs/{sample}_to_{asmbl}_back_depth.log"
    shell:
        "samtools view -bh {input.sam} -@ {threads} | samtools sort -T {params.tmp} - > {output.bam}; "
        "export OMP_NUM_THREADS={threads}; "
        "jgi_summarize_bam_contig_depths --outputDepth {output.depth} {output.bam}"

# This rule merges the depth files
rule merge_depths_viral: 
    input: 
        depth_files = expand("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_{{asmbl}}_viral_contigs.depth", sample=config["samples"])
    output:
        depth_file_merged="../results/vMAG_binning/backmapping/merged_depths/{asmbl}_global_viral_depth.txt"
    resources:
        account="pengel_beemicrophage",
        runtime="1:00:00",
        mem_mb = 10000
    threads: 1
    conda: 
        "envs/mags_env.yaml"
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
        assembly="../results/assembly/viral_contigs/single_sample/{asmbl}_viral_contigs.fasta",
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

####################################################################################
################################## POLISH VIRAL CONTIGS ############################
####################################################################################
# The foilowing rules will:

# 1. remove bacterial contigs (size > 500kb and bacteria according to checkM)
# 2. dereplicate for identical contigs across samples
# 3. run CheckV again to have a last table with all the phages
################################################################################
################################################################################

# this rule bins the viral contigs using the bins from vRhyme
rule bin_viral_contigs:
    input:
        assembly = "../results/assembly/viral_contigs/single_sample/{sample}_viral_contigs.fasta",
        bins="../results/vMAG_binning/vRhyme_bins/{sample}"
    output:
        binned_fasta="../results/assembly/viral_contigs/binning/{sample}_viral_contigs_binned.fasta",
        binning_data="../results/vMAG_binning/summary_tables/{sample}_binning_data.tsv"
    resources:
        account="pengel_beemicrophage",
        runtime="1:00:00",
        mem_mb = 100000
    threads: 1
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/vMAGs/vRhyme/{sample}_binning_fasta.log"
    shell:
        "python scripts/vMAGs/bin_contigs.py -i {input.assembly}  -d {input.bins} -f {output.binned_fasta} -t {output.binning_data}"

# this rule runs checkV on the binned viral contigs
rule run_checkv:
    input:
        assembly = "../results/assembly/viral_contigs/binning/{sample}_viral_contigs_binned.fasta"
    output:
        dir=temp(directory("../results/vMAG_binning/polishing/checkv_viral_contigs/{sample}_checkv")),
        trimmed="../results/assembly/viral_contigs/trimmed/{sample}_viral_contigs_trimmed.fasta",
        tab="../results/vMAG_binning/polishing/checkv/{sample}_checkv_quality_summary.tsv"
    params:
        db="resources/default_DBs/checkv-db-v1.5/"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "02:00:00"
    threads:20
    conda:
        "envs/checkv.yaml"
    log:
        "logs/polish_vMAGs/checkV/{sample}_checkV.log"
    benchmark:
        "logs/polish_vMAGs/checkV/{sample}_checkV.benchmark"
    shell:
        "checkv end_to_end {input.assembly} {output.dir} -t {threads} -d {params.db}; "
        "mv {output.dir}/quality_summary.tsv {output.tab}; "
        "cat {output.dir}/viruses.fna {output.dir}/proviruses.fna > {output.trimmed}"

# this rule filters out contigs that are not viral (found with low confidence by only one tool, too long or too short)
rule viral_contigs_quality_filtering:
    input:
        assembly="../results/assembly/viral_contigs/trimmed/{sample}_viral_contigs_trimmed.fasta",
        binning_data="../results/vMAG_binning/summary_tables/{sample}_binning_data.tsv",
        vi="../results/viral_identification/ViralIdentification_scores.tsv",
        checkv="../results/vMAG_binning/polishing/checkv/{sample}_checkv_quality_summary.tsv"
    output:
        filtered_fasta="../results/assembly/viral_contigs/QC_filtered/{sample}_viral_contigs_filtered.fasta",
        filtering_data=temp("../results/vMAG_binning/summary_tables/{sample}_viral_contigs_metadata.tsv")
    resources:
        account="pengel_beemicrophage",
        runtime="1:00:00",
        mem_mb = 100000
    threads: 1
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/vMAGs/polishing/{sample}_filtering_vMAGs.log"
    benchmark:
        "logs/vMAGs/vRhyme/{sample}_filtering_vMAGs.benchmark"
    shell:
        "python scripts/vMAGs/parse_vMAGs_binning.py -i {input.assembly} -b {input.binning_data} -v {input.vi} -c {input.checkv} -f {output.filtered_fasta}  -t {output.filtering_data}"

rule split_contigs_checkm:
    input:
        "../results/assembly/viral_contigs/QC_filtered/{sample}_viral_contigs_filtered.fasta",
    output:
        tmp_dir=temp(directory("../results/vMAG_binning/polishing/{sample}_single_genomes")),
    resources:
        account="pengel_beemicrophage",
        runtime="0:30:00",
        mem_mb = 8000
    threads: 1
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/polish_vMAGs/checkm/{sample}_split_for_checkm.log"
    shell:
        "python scripts/assembly/split_assembly.py -f {input} -d {output}"

rule checkm_viral_contigs:
    input:
        dir="../results/vMAG_binning/polishing/{sample}_single_genomes"
    output:
        dir=directory("../results/vMAG_binning/polishing/checkm/{sample}_checkm"),
        file="../results/vMAG_binning/polishing/checkm/{sample}_checkm_stats.tsv"
    log:
        "logs/polish_vMAGs/checkm/{sample}_checkm.log"
    benchmark:
        "logs/polish_vMAGs/checkm/{sample}_checkm.benchmark"
    threads: 25
    params:
        db="resources/default_DBs/checkm_db"
    conda:
        "envs/checkm_env.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 200000,
        runtime= "04:00:00"
    shell:
        "export CHECKM_DATA_PATH={params.db}; "
        "checkm lineage_wf {input} {output.dir} -x fasta -t {threads}; "
        "checkm qa {output.dir}/lineage.ms {output.dir} -o 2 -f {output.file} --tab_table"

rule filter_bacteria_from_viral_contigs:
    input:
        fasta="../results/assembly/viral_contigs/QC_filtered/{sample}_viral_contigs_filtered.fasta",
        checkm="../results/vMAG_binning/polishing/checkm/{sample}_checkm_stats.tsv",
        mtdata="../results/vMAG_binning/summary_tables/{sample}_viral_contigs_metadata.tsv"
    output:
        fasta_nobact="../results/assembly/viral_contigs/polished/{sample}_viral_contigs_filtered_nobact.fasta",
        mtdata_nobact="../results/vMAG_binning/summary_tables/{sample}_viral_contigs_metadata_nobact.tsv"
    resources:
        account="pengel_beemicrophage",
        runtime="0:30:00",
        mem_mb = 8000
    threads: 1
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/polish_vMAGs/checkm/{sample}_remove_bacteria.log"
    shell:
        "python scripts/vMAGs/filter_bacteria_from_viral_contigs.py -f {input.fasta} -c {input.checkm} -m {input.mtdata} -o {output.fasta_nobact} -t {output.mtdata_nobact}"

rule dereplicate_viral_contigs:
    input:
        fasta=expand("../results/assembly/viral_contigs/polished/{sample}_viral_contigs_filtered_nobact.fasta", sample=config["samples"]),
    output:
        all_fasta="../results/assembly/viral_contigs/dereplicated/all_viral_contigs.fasta",
        concat_cont_drep="../results/assembly/viral_contigs/dereplicated/all_viral_contigs_drep.fasta",
        concat_cont_clst=temp("../results/assembly/viral_contigs/dereplicated/all_viral_contigs_drep.fasta.clstr")
    resources:
        account="pengel_beemicrophage",
        runtime="04:00:00",
        mem_mb = 200000
    threads: 25
    log:
        "logs/polish_vMAGs/dereplication/drep.log"
    benchmark:
        "logs/polish_vMAGs/dereplication/drep.benchmark"
    conda:
        "envs/cd-hit.yaml"
    shell:
        "cat {input.fasta} > {output.all_fasta}; "
        "cd-hit-est -i {output.all_fasta} -o {output.concat_cont_drep} -c 0.95 -aS 0.85 -M 0 -d 0 -T {threads}"

rule align_viral_contigs:
    input:
        fasta="../results/assembly/viral_contigs/dereplicated/all_viral_contigs.fasta",
    output:
        directory("../results/assembly/viral_contigs/dereplicated/all_genomes_alignment")
    resources:
        account="pengel_beemicrophage",
        runtime="01:00:00",
        mem_mb = 200000
    threads: 16
    log:
        "logs/polish_vMAGs/dereplication/align_viral_genomes.log"
    benchmark:
        "logs/polish_vMAGs/dereplication/align_viral_genomes.benchmark"
    conda:
        "envs/drep_env.yaml"
    shell:
        "python scripts/vMAGs/align_viral_contigs.py -i {input} -o {output} -t {threads}"

rule parse_drep_viral_contigs:
    input:
        derep_fasta="../results/assembly/viral_contigs/dereplicated/all_viral_contigs_drep.fasta.clstr"
    output:
        "../results/vMAG_binning/summary_tables/all_viral_contigs_dereplication.tsv"
    resources:
        account="pengel_beemicrophage",
        runtime="01:00:00",
        mem_mb = 20000
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/polish_vMAGs/dereplication_parsing.log"
    benchmark:
        "logs/polish_vMAGs/dereplication_parsing.benchmark"
    shell:
        "python scripts/vMAGs/parse_cdhit.py -i {input} -o {output}" 

rule aggreagte_polishing_tables:
    input:
        checkv= expand("../results/vMAG_binning/polishing/checkv/{sample}_checkv_quality_summary.tsv", sample=config["samples"]),
        checkm= expand("../results/vMAG_binning/polishing/checkm/{sample}_checkm_stats.tsv", sample=config["samples"]),
        mtdata_nobact=expand("../results/vMAG_binning/summary_tables/{sample}_viral_contigs_metadata_nobact.tsv", sample=config["samples"]),
        drep="../results/vMAG_binning/summary_tables/all_viral_contigs_dereplication.tsv"
    output:
        all_checkv="../results/vMAG_binning/polishing/checkv/all_checkv_quality_summary.tsv",
        all_checkm="../results/vMAG_binning/polishing/checkm/all_checkm_stats.tsv",
        all_mtdata_nobact="../results/vMAG_binning/summary_tables/all_viral_contigs_metadata.tsv"
    resources:
        account="pengel_beemicrophage",
        runtime="01:00:00",
        mem_mb = 20000
    log:
        "logs/polish_vMAGs/aggregrate_polishing_tables.log"
    benchmark:
        "logs/polish_vMAGs/aggregrate_polishing_tables.benchmark"
    run:
        import pandas as pd
        # open drep
        drep=pd.read_csv(input.drep, sep="\t")
        # keep only rows where reference=True
        drep=drep[drep["reference"]==True]
        # store contig_id in a list
        drep_contigs=drep["contig_id"].tolist()
        checkV=concat_tables(input.checkv)
        checkM=concat_tables(input.checkm)
        mtdata_nobact=concat_tables(input.mtdata_nobact)
        checkV.to_csv(output.all_checkv, sep="\t", index=False)
        checkM.to_csv(output.all_checkm, sep="\t", index=False)
        # add column "reference" to mtdata_nobact
        mtdata_nobact["reference"]=mtdata_nobact["genome_id"].isin(drep_contigs)
        mtdata_nobact.to_csv(output.all_mtdata_nobact, sep="\t", index=False)