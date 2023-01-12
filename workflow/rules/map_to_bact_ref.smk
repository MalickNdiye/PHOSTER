rule concat_genomes:
    input:
        refs="../results/reference_db/dereplicated_genomes_filtered/"
    output:
        concat="../results/reference_db/all_genomes/all_bacterial_RefGenomes.fasta",
        stb="../results/reference_db/all_genomes/all_bacterial_RefGenomes.stb"
    log:
        "logs/ref_db/concat_refs.log"
    threads: 2
    conda:
        "envs/drep_env.yaml"
    params:
        refs=get_files_commas("../results/reference_db/dereplicated_genomes_filtered/", sep=" ")
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1500,
        runtime= "00:30:00"
    shell:
        "cat {input.refs}/*.f* >> {output.concat}; "
        "parse_stb.py --reverse -f {params.refs}  -o {output.stb}"

rule build_ref_index:
    input:
        ref="../results/reference_db/all_genomes/all_bacterial_RefGenomes.fasta"
    output:
        index=directory("../results/reference_db/bacteria_index/")
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
        index="../results/reference_db/bacteria_index/",
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
