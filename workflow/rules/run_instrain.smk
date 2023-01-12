rule generate_genelist:
    input:
        ref="../results/reference_db/all_genomes/all_bacterial_RefGenomes.fasta"
    output:
        "../results/inStrain/all_bacterial_RefGenomes_genes.faa"
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
        genL="../results/inStrain/all_bacterial_RefGenomes_genes.faa",
        stb="../results/reference_db/all_genomes/all_bacterial_RefGenomes.stb",
    output:
        directory("../scratch_link/inStrain/profile_{type}/{sample}{type}_profile/")
    threads: 25
    conda:
        "envs/inStrain.yaml"
    log:
        "logs/instrain/{sample}{type}_profile.log"
    benchmark:
        "logs/instrain/{sample}{type}_profile.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 200000,
        runtime= "72:00:00"
    shell:
        "(inStrain profile {input.bam} {input.ref} -o {output} -p {threads} -g {input.genL} -s {input.stb})2> {log}"

rule instrain_compare:
    input:
        IS=expand("../scratch_link/inStrain/profile_{{type}}/{sample}{{type}}_profile/", sample=config["sam_names"]),
        ref="../results/reference_db/all_genomes/all_bacterial_RefGenomes.fasta",
        genL="../results/inStrain/all_bacterial_RefGenomes_genes.faa",
        stb="../results/reference_db/all_genomes/all_bacterial_RefGenomes.stb",
    output:
        directory("../scratch_link/inStrain/compare_{type}")
    threads: 46
    conda:
        "envs/inStrain.yaml"
    log:
        "logs/instrain/{type}_compare.log"
    benchmark:
        "logs/instrain/{type}_compare.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 200000,
        runtime= "72:00:00"
    shell:
        "(inStrain compare -i {input.IS} -o {output} -p {threads} -g {input.genL} -s {input.stb})2> {log}"
