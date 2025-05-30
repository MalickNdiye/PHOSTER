###############################################################################
######################### BACTERIA PANGENOME ANALYISIS #########################
#Once our database is well curated, we can TODO
###############################################################################

rule find_genomes_cluster:
    input:
        clust_info="../results/reference_db_filtered/summary_data_tables/clust_filtered.tsv",
        filtered_mags="../results/MAG_binning/bins/filtered_mags/",
        refs_isolates="../data/reference_assemblies/hb_bacteria/non_redundant/single_genomes/contigs/"
    output:
        paths="../results/pangenomes/{type}/{genus}_paths.txt",
        tmp_db=temp(directory("../scratch_link/tmp_ref/reference_genomes_{type}_tmp_{genus}/"))
    log:
        "logs/{type}_pangenomes/{genus}/find_cluster_paths.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 2000,
        runtime= "00:30:00"
    shell:
        "mkdir -p {output.tmp_db}; "
        "cp {input.filtered_mags}/*.fa {output.tmp_db}; "
        "cp {input.refs_isolates}/*.fna {output.tmp_db}; "
        "python scripts/bacteria_pangenome/find_cluster_genomes.py {wildcards.genus} {output.tmp_db} {input.clust_info} {output.paths}"

rule annotate_clusters:
    input:
        genus_paths="../results/pangenomes/{type}/{genus}_paths.txt",
        tmp_db="../scratch_link/tmp_ref/reference_genomes_{type}_tmp_{genus}/"
    output:
        fnas=directory("../results/pangenomes/{type}/{genus}/genes_fna"),
        faas=directory("../results/pangenomes/{type}/{genus}/genes_faa")
    conda:
        "envs/orthofinder.yaml"
    log:
        "logs/{type}_pangenomes/{genus}/prodiagal_annot.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    shell:
        "mkdir -p {output.faas}; "
        "mkdir -p {output.fnas}; "
        "scripts/bacteria_pangenome/run_prod_annot_clusters.sh {output.fnas} {output.faas} {input.genus_paths}"

rule run_orthofinder:
    input:
        faas="../results/pangenomes/{type}/{genus}/genes_faa",
    output:
        ortho_out=directory("../scratch_link/pangenomes/{type}/{genus}/orthofinder_output/")
    conda:
        "envs/orthofinder.yaml"
    log:
        "logs/{type}_pangenomes/{genus}/run_orthofinder.log"
    params:
        name="results",
        ulim=40000
    threads: 15
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "24:00:00"
    shell:
        "ulimit -n {params.ulim}; "
        "orthofinder  -f {input.faas} -o {output.ortho_out} -n {params.name} -t {threads}"

rule run_mOTUpan:
    input:
        faas="../results/pangenomes/{type}/{genus}/genes_faa",
        ref_db="../results/reference_db"
    output:
        pan="../results/pangenomes/{type}/{genus}/{genus}_mOTUpan.tsv",
        chemck_out=temp("../results/pangenomes/{type}/{genus}/checkm_out_tab.tsv")
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/{type}_pangenomes/{genus}/run_mOTUpan.log"
    benchmark:
        "logs/{type}_pangenomes/{genus}/run_mOTUpan.benchmark"
    threads: 15
    params:
        boots=100,
        tmp1="../results/pangenomes/{type}/{genus}/tmp1.tsv",
        tmp2="../results/pangenomes/{type}/{genus}/tmp2.tsv"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "03:00:00"
    shell:
        """awk -F',' -v OFS=',' '{{ gsub(/.fa/,"_genes", $1); print }} ' {input.ref_db}/data_tables/Chdb.csv > {params.tmp1}; """
        """awk -F',' -v OFS=',' '{{ gsub(/.fna/,"_genes", $1); print }} ' {params.tmp1} > {params.tmp2}; """
        "sed -e 's/,/\t/g' {params.tmp2} > {output.chemck_out}; "
        "rm {params.tmp1} {params.tmp2}; "
        "mOTUpan.py --faas {input.faas}/* --boots {params.boots} -o {output.pan} --checkm {output.chemck_out} --threads {threads}"


rule single_OGs_diversity_parse:
    input:
        og="../scratch_link/pangenomes/{type}/{genus}/orthofinder_output/",
        clust_info="../results/reference_db_filtered/summary_data_tables/clust_filtered.tsv"
    output:
        genus_OGs="../results/pangenomes/{type}/{genus}/{genus}_single_copy_OGs.tsv",
        isolates_OGs="../results/pangenomes/{type}/{genus}/isolates_{genus}_single_copy_OGs.tsv"
    log:
        "logs/{type}_pangenomes/{genus}/parse_singOG_isolates.log"
    params:
        genus=lambda wildcards: wildcards.genus
    threads: 3
    conda:
        "envs/drep_env.yaml"
    resources:
        account="pengel_beemicrophage",
        mem_mb=50000,
        runtime="00:30:00"
    script:
        "scripts/bacteria_pangenome/Find_genus_singleCopy_OGs.py"

rule single_OGs_diversity_aggregate:
    input:
        genus_OGs=expand("../results/pangenomes/{{type}}/{genus}/{genus}_single_copy_OGs.tsv", genus=core_genera),
        isolates_OGs=expand("../results/pangenomes/{{type}}/{genus}/isolates_{genus}_single_copy_OGs.tsv", genus=core_genera)
    output:
        all_genus_OGs="../results/pangenomes/{type}/all_single_copy_OGs.tsv",
        all_isolates_OGs="../results/pangenomes/{type}/all_isolates_single_copy_OGs.tsv"
    log:
        "logs/{type}_pangenomes/aggregate_singOG_isolates.log"
    threads: 3
    resources:
        account="pengel_beemicrophage",
        mem_mb=5000,
        runtime="00:30:00"
    shell:
        "echo -e 'Orthogroup\tgenome\tgene\tBin_Id\tsecondary_cluster\tspecies\tgenus\tsecondary_cluster_n' > {output.all_genus_OGs}; "
        "echo -e 'Orthogroup\tgenome\tgene\tBin_Id\tsecondary_cluster\tspecies\tgenus\tsecondary_cluster_n' > {output.all_isolates_OGs}; "
        "awk 'FNR>1' {input.genus_OGs} >> {output.all_genus_OGs}; "
        "awk 'FNR>1' {input.isolates_OGs} >> {output.all_isolates_OGs}"