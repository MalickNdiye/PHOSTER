rule generate_genelist:
    input:
        ref="../results/reference_db/all_genomes/all_bacterial_RefGenomes.fasta"
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
        ref="../results/reference_db/all_genomes/all_bacterial_RefGenomes.fasta",
        genL="../results/inStrain/all_bacterial_RefGenomes_genes.fna",
        stb="../results/reference_db/all_genomes/all_bacterial_RefGenomes.stb",
    output:
        directory("../scratch_link/inStrain/profile_{type}/{sample}{type}_profile/")
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
        runtime= "72:00:00"
    params:
        min_ANI=0.92
    shell:
        "(inStrain profile {input.bam} {input.ref} -o {output} --min_read_ani {params.min_ANI} -p {threads} -g {input.genL} -s {input.stb})2> {log}"

rule aggregate_instain_tabs:
    input:
        IS=expand("../scratch_link/inStrain/profile_{{type}}/{sample}{{type}}_profile", sample=config["sam_names"])
        tax="../results/reference_db_filtered/summary_data_tables/clust_filtered_winners.tsv"
    output:
        genome_info="../results/inStrain/data_tables/all_genome_info.tsv",
        gene_info="../results/inStrain/data_tables/all_gene_info.tsv",
        mapping_info="../results/inStrain/data_tables/all_mapping_info.tsv",
        scaffold_info="../results/inStrain/data_tables/all_scaffold_info.tsv",
    conda:
        "envs/base_R_env"
    params:
        us_func="scripts/useful_func.R"
    log:
        "logs/instrain/aggreagte_profile_data_tabs.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "01:00:00"
    script:
        "scripts/bacteria_community_analysis/aggregate_instrain_tabs.R"

rule instrain_compare:
    input:
        IS=expand("../scratch_link/inStrain/profile_{{type}}/{sample}{{type}}_profile/", sample=config["sam_names"]),
        ref="../results/reference_db/all_genomes/all_bacterial_RefGenomes.fasta",
        stb="../results/reference_db/all_genomes/all_bacterial_RefGenomes.stb"
    output:
        directory("../scratch_link/inStrain/compare_{type}/")
    threads: 46
    conda:
        "envs/inStrain.yaml"
    log:
        "logs/instrain/{type}_compare.log"
    benchmark:
        "logs/instrain/{type}_compare.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "100:00:00"
    shell:
        "(inStrain compare -i {input.IS} -o {output} -p {threads} -s {input.stb})2> {log}"
