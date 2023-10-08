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
        "python scripts/vMAGs/bin_contgs.py -i {input.assembly}  -d {input.bins} -f {output.binning_data} -t {output.binned_fasta}"

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
        "python scripts/vMAGs/parse_vMAGs_binning.py -i {input.assembly} -b {input.binning_data} -v {input.vi} -c {input.checkV} -f {output.filtered_fasta}  -t {output.filtering_data}"


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
        "cd-hit-est -i {input} -o {output.concat_cont_drep} -c 0.95 -aS 0.85 -M 0 -d 0 -T {threads}"

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
        mtdata_nobact=expand("../results/vMAG_binning/summary_tables/{sample}_viral_contigs_metadata_nobact.tsv", sample=config["samples"])
    output:
        all_checkv="../results/vMAG_binning/polishing/checkv/all_checkv_quality_summary.tsv",
        all_checkm="../results/vMAG_binning/polishing/checkm/all_checkm_stats.tsv",
        all_mtdata_nobact="../results/vMAG_binning/summary_tables/all_viral_contigs_metadata.tsv"
    resources:
        account="pengel_beemicrophage",
        runtime="01:00:00",
        mem_mb = 20000
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/polish_vMAGs/aggregrate_polishing_tables.log"
    benchmark:
        "logs/polish_vMAGs/aggregrate_polishing_tables.benchmark"
    run:
        checkV=concat_tables(input.checkv)
        checkM=concat_tables(input.checkm)
        mtdata_nobact=concat_tables(input.mtdata_nobact)
        checkV.to_csv(output.all_checkv, sep="\t", index=False)
        checkM.to_csv(output.all_checkm, sep="\t", index=False)
        mtdata_nobact.to_csv(output.all_mtdata_nobact, sep="\t", index=False)




