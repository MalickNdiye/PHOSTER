def get_all_SecCluster(path, ignore_sing=True):
    checkpoint_output = path
    df=pd.read_csv(path, delimiter="\t")
    secodary_clusters=list(df["secondary_cluster"])

    if ignore_sing:
        secodary_clusters=set([i for i in secodary_clusters if secodary_clusters.count(i)>1])
    return(secodary_clusters)

rule find_genomes_cluster:
    input:
        clust_info="../results/reference_db/summary_data_tables/clust_filtered.tsv",
        filtered_mags="../results/MAG_binning/bins/filtered_mags/",
        refs_isolates="../data/reference_assemblies/hb_bacteria/non_redundant/single_genomes/contigs/"
    output:
        paths="../results/pangenomes/{type}/{cluster}_paths.txt",
        tmp_db=temp(directory("../scratch_link/tmp_ref/reference_genomes_{type}_tmp_{cluster}/"))
    log:
        "logs/{type}_pangenomes/{cluster}/find_cluster_paths.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 2000,
        runtime= "00:30:00"
    shell:
        "mkdir -p {output.tmp_db}; "
        "cp {input.filtered_mags}/*.fa {output.tmp_db}; "
        "cp {input.refs_isolates}/*.fna {output.tmp_db}; "
        "python scripts/bacteria_pangenome/find_cluster_genomes.py {wildcards.cluster} {output.tmp_db} {input.clust_info} {output.paths}"

rule annotate_clusters:
    input:
        clust_paths="../results/pangenomes/{type}/{cluster}_paths.txt",
        tmp_db="../scratch_link/tmp_ref/reference_genomes_{type}_tmp_{cluster}/"
    output:
        fnas=directory("../results/pangenomes/{type}/{cluster}/genes_fna"),
        faas=directory("../results/pangenomes/{type}/{cluster}/genes_faa")
    conda:
        "envs/orthofinder.yaml"
    log:
        "logs/{type}_pangenomes/{cluster}/prodiagal_annot.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    shell:
        "mkdir -p {output.faas}; "
        "mkdir -p {output.fnas}; "
        "scripts/bacteria_pangenome/run_prod_annot_clusters.sh {output.fnas} {output.faas} {input.clust_paths}"

rule run_orthofinder:
    input:
        faas="../results/pangenomes/{type}/{cluster}/genes_faa",
    output:
        ortho_out=directory("../scratch_link/pangenomes/{type}/{cluster}/orthofinder_output/")
    conda:
        "envs/orthofinder.yaml"
    log:
        "logs/{type}_pangenomes/{cluster}/run_orthofinder.log"
    params:
        name="results",
        ulim=3000
    threads: 15
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "07:00:00"
    shell:
        "ulimit -n {params.ulim}; "
        "orthofinder  -f {input.faas} -o {output.ortho_out} -n {params.name} -t {threads}"

rule single_OGs_diversity_parse:
    input:
        og="../scratch_link/pangenomes/{type}/{cluster}/orthofinder_output/",
        IS_profile=expand("../scratch_link/inStrain/profile_{{type}}/{sample}{{type}}_profile/", sample=config["sam_names"]),
        stb="../results/reference_db/all_genomes/all_bacterial_RefGenomes.stb"
    output:
        ortho_div="../results/pangenomes/{type}/{cluster}/single_copy_OGs_diversity.tsv"
    log:
        "logs/{type}_pangenomes/{cluster}/run_orthofinder.log"
    params:
        clust=lambda wildcards: wildcards.cluster
    threads: 3
    conda:
     "envs/drep_env.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "00:30:00"
    script:
        "scripts/bacteria_pangenome/FInd_cluster_singleCopy_OGs.py"

rule aggregate_single_OGs_diversity:
    input:
        ortho_div=expand("../results/pangenomes/{{type}}/{cluster}/single_copy_OGs_diversity.tsv", cluster=get_all_SecCluster(rules.parse_dRep.output.clust_final))
    output:
        ortho_div="../results/pangenomes/{type}/all_single_copy_OGs_diversity.tsv"
    log:
        "logs/{type}_pangenomes/aggregate_coreOGs.log"
    threads: 3
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "00:30:00"
    script:
        "scripts/bacteria_pangenome/aggregate_coreOG_diversity.py"
