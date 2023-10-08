################################################################################
######################### MAP TO REFERENCE DB ##################################
#Once our database is well curated, we can TODO
################################################################################
rule concat_genomes:
    input:
        refs="../results/reference_db_filtered/dereplicated_genomes_filtered/"
    output:
        concat="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.fasta",
        stb="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.stb"
    log:
        "logs/ref_db/concat_refs.log"
    threads: 2
    conda:
        "envs/drep_env.yaml"
    params:
        refs=lambda wildcards, input: get_files_commas(input[0], sep=" ")
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1500,
        runtime= "00:30:00"
    shell:
        "cat {input.refs}/*.f* >> {output.concat}; "
        "parse_stb.py --reverse -f {params.refs}  -o {output.stb}"

rule build_ref_index:
    input:
        ref="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.fasta"
    output:
        index=directory("../results/reference_db_filtered/bacteria_index/")
    conda:
        "envs/map_env.yaml"
    threads: 25
    log:
        "logs/mapping/build_bowtie_index.log"
    params:
        basename="all_bacterial_RefGenomes",
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 20000,
        runtime= "01:00:00"
    shell:
        "mkdir -p {output.index}; "
        "bowtie2-build {input.ref} {output.index}/{params.basename} --threads {threads}"

rule MapReads:
    input:
        index="../results/reference_db_filtered/bacteria_index/",
        R1="../data/host_filtered_reads/{sample}{type}_R1_HF.fastq.gz",
        R2="../data/host_filtered_reads/{sample}{type}_R2_HF.fastq.gz",
    output:
        sam=temp("../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.sam"),
    wildcard_constraints:
        sample="\d+"
    conda:
        "envs/map_env.yaml"
    threads: 25
    params:
        basename="all_bacterial_RefGenomes"
    log:
        "logs/mapping/map/{sample}{type}_bbmap_mapping.log"
    benchmark:
        "logs/mapping/map/{sample}{type}_bbmap_mapping.benchmark"
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 200000,
        runtime= "07:00:00"
    shell:
        "bowtie2 -x {input.index}/{params.basename} -1 {input.R1} -2 {input.R2} -S {output.sam} --threads {threads}"

rule generate_bam:
    input:
        sam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.sam"
    output:
        bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.bam",
        cf=temp("../results/mapping/{sample}{type}_library_count.tsv")
    wildcard_constraints:
        sample="\d+"
    conda:
        "envs/map_env.yaml"
    threads: 2
    log:
        "logs/mapping/{sample}{type}_prepare_metapop.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "01:00:00"
    shell:
        "samtools view -S -b {input.sam} > {output.bam}; "
        "touch {output.cf}; "
        "s=$(basename {output.bam}); "
        "lib=${{s%.*}}; "
        "count=$(samtools view -c {output.bam}); "
        'echo -e "${{lib}}\t${{count}}" >> {output.cf}'

##################### Mapping QC ###############################################

rule mappings_stats:
    input:
        bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.bam",
    output:
        temp("../results/mapping/mapout_{sample}{type}/{sample}{type}_mapstats.tsv")
    wildcard_constraints:
        sample="\d+"
    conda:
        "envs/map_env.yaml"
    log:
        "logs/mapping/mapstats/{sample}{type}_mapstats.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "01:00:00"
    threads: 8
    shell:
        "samtools flagstat {input.bam} -@ {threads} > {output}"

acc="B"

rule merge_mapping_stats:
    input:
        filelist=expand("../results/mapping/mapout_{sample}{{type}}/{sample}{{type}}_mapstats.tsv", sample=config["sam_names"])
    output:
        "../results/mapping/all_{type}_mapstats.tsv"
    wildcard_constraints:
        sample="\d+"
    conda:
        "envs/base_R_env.yaml"
    log:
        "logs/mapping/{type}_merge_mapstats.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:15:00"
    threads: 2
    script:
        "scripts/mapping/format_flagstats.R"

rule sort_bam:
    input:
        bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.bam"
    output:
        sorted_bam="../scratch_link/mapping/sorted_bams/{sample}{type}_mapping_sorted.bam",
        bai="../scratch_link/mapping/sorted_bams/{sample}{type}_mapping_sorted.bam.bai"
    log:
        "logs/mapping/QC/sort_bam_{sample}{type}.log"
    conda:
        "envs/map_env.yaml"
    threads: 5
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "01:30:00"
    shell:
        "samtools sort {input.bam} -@ {threads} -o {output.sorted_bam}; "
        "samtools index {output.sorted_bam} {output.bai} -@ {threads}"

rule multireads_count:
    input:
        sorted_bam="../scratch_link/mapping/sorted_bams/{sample}{type}_mapping_sorted.bam"
    output:
        multireads_count="../results/mapping/multireads/{sample}{type}_multireads_count.tsv"
    threads: 1
    log:
        "logs/mapping/multireads/{sample}{type}_multireads.log"
    benchmark:
        "logs/mapping/multireads/{sample}{type}_multireads.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "03:00:00"
    script:
        "scripts/mapping/multireads.py"


rule aggregate_multireads:
    input:
        multireads_count=expand("../results/mapping/multireads/{sample}{{type}}_multireads_count.tsv", sample=config["sam_names"])
    output:
        all_multireads="../results/mapping/multireads/all_{type}_multireads_count.tsv"
    threads: 1
    log:
        "logs/mapping/multireads/all_{type}_multireads.log"
    benchmark:
        "logs/mapping/multireads/all_{type}_multireads.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "00:15:00"
    shell:
        "echo -e 'scaffold\tmultireads\tavg_mapQ\tsample' > {output.all_multireads}; "
        "tail -n +2 {input.multireads_count} >> {output.all_multireads}"

rule prepare_bamQC:
    input:
        bam=expand("../scratch_link/mapping/sorted_bams/{sample}{{type}}_mapping_sorted.bam", sample=config["sam_names"])
    output:
        "../results/mapping/QC_{type}/bamqc_config.txt"
    log:
        "logs/mapping/QC/prepare_{type}_bamQC.log"
    threads: 2
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1000,
        runtime= "00:15:00"
    script:
        "scripts/mapping/prepare_bamQC.py"

rule bamQC:
    input:
        config="../results/mapping/QC_{type}/bamqc_config.txt"
    output:
        directory("../results/mapping/QC_{type}/multi_bamQC/")
    log:
        "logs/mapping/QC/multi_{type}_bamQC.log"
    benchmark:
        "logs/mapping/QC/multi_{type}_bamQC.benchmark"
    threads: 12
    params:
        java_mem="300G"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 300000,
        runtime= "05:30:00"
    conda:
        "envs/qualimap.yaml"
    shell:
        "qualimap multi-bamqc -d {input.config} -outdir {output} -r -c --java-mem-size={params.java_mem}"


################################################################################
################################## INSTRAIN ####################################
################################################################################

rule generate_genelist:
    input:
        ref="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.fasta"
    output:
        "../results/inStrain/all_bacterial_RefGenomes_genes.fna"
    conda:
        "envs/inStrain.yaml"
    threads: 25
    log:
        "logs/instrain/generate_gene_list.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "1:00:00"
    shell:
        "prodigal -i {input.ref} -d {output}"

rule instrain_profile:
    input:
        bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.bam",
        ref="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.fasta",
        genL="../results/inStrain/all_bacterial_RefGenomes_genes.fna",
        stb="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.stb",
    output:
        dir=directory("../scratch_link/inStrain/profile_{type}/{sample}{type}_profile/"),
        sorted_bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.sorted.bam"
    threads: 32
    conda:
        "envs/inStrain.yaml"
    log:
        "logs/instrain/{sample}{type}_profile.log"
    benchmark:
        "logs/instrain/{sample}{type}_profile.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "10:00:00"
    params:
        min_ANI=0.92
    shell:
        "(inStrain profile {input.bam} {input.ref} -o {output.dir} --min_read_ani {params.min_ANI} -p {threads} -g {input.genL} -s {input.stb})2> {log}"


rule aggregate_instain_tabs:
    input:
        IS=expand("../scratch_link/inStrain/profile_{{type}}/{sample}{{type}}_profile", sample=config["sam_names"]),
        tax="../results/reference_db_filtered/summary_data_tables/clust_filtered_winners.tsv"
    output:
        genome_info="../results/inStrain/{type}/data_tables/all_genome_info.tsv",
        gene_info="../results/inStrain/{type}/data_tables/all_gene_info.tsv",
        mapping_info="../results/inStrain/{type}/data_tables/all_mapping_info.tsv",
        scaffold_info="../results/inStrain/{type}/data_tables/all_scaffold_info.tsv",
    conda:
        "envs/base_R_env.yaml"
    params:
        us_func="scripts/useful_func.R"
    log:
        "logs/instrain/{type}/aggreagte_profile_data_tabs.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "01:00:00"
    script:
        "scripts/bacteria_community_analysis/aggregate_instrain_tabs.R"

rule do_rarefaction:
    input:
        IS=expand("../scratch_link/inStrain/profile_{{type}}/{sample}{{type}}_profile", sample=config["sam_names"]),
        clust_info="../results/reference_db_filtered/summary_data_tables/clust_filtered_winners.tsv",
        clust_filt="../results/reference_db_filtered/summary_data_tables/clust_filtered.tsv",
        genome_info="../results/inStrain/{type}/data_tables/all_genome_info.tsv",
        gene_info="../results/inStrain/{type}/data_tables/all_gene_info.tsv",
        stb="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.stb",
        motupan=expand("../results/pangenomes/{{type}}/{genus}/{genus}_mOTUpan.tsv", genus=core_genera)
    output:
        rare="../results/inStrain/{type}/data_tables/snv_counts_mOTUpanCore.tsv"
    conda:
        "envs/rare_R.yaml"
    params:
        us_func="scripts/useful_func.R"
    log:
        "logs/instrain/{type}/do_rarefaction.log"
    threads: 40
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 200000,
        runtime= "01:00:00"
    script:
        "scripts/bacteria_community_analysis/do_rarefaction.R"

rule instrain_compare:
    input:
        IS=expand("../scratch_link/inStrain/profile_{{type}}/{sample}{{type}}_profile/", sample=config["sam_names"]),
        ref="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.fasta",
        stb="../results/reference_db_filtered/all_genomes/all_bacterial_RefGenomes.stb"
    output:
        directory("../results/inStrain/{type}/compare_{type}/compare_{genome}")
    threads: 40
    conda:
        "envs/inStrain.yaml"
    log:
        "logs/instrain/compare/{type}_compare_{genome}.log"
    benchmark:
        "logs/instrain/compare/{type}_compare_{genome}.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 600000,
        runtime= "3-00:00:00"
    params:
        genome=lambda wildcards: wildcards.genome
    shell:
        "inStrain compare -i {input.IS} -o {output} -p {threads} -s {input.stb} --database_mode --genome {params.genome} || mkdir -p {output}/singleton"
        # Since I feed all the genomes of the rference database to this command, if the genom is not present in more than on sample, the command will fail.
        # As a solution, I catch the error and make the output directory anyways to trick snakemake
        # This is not the best solution. TODO make something more elegeant