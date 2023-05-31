###############################################################################
########################## Bacterial Reference Database ##################################
################################################################################
# Once we get all our good quality MAGs, it is time to create a reference Database.
# To do so, all filtered MAGs will be aggreagted in a directory with 211 genomes
# from isolates of A. mellifera and A. cerana gut microbiota. Then, Genomes will
# be dereplicated at 95% ANI using dRep. This will yield a species-level database,
# where every entry should correspond to a species.
################################################################################
################################################################################

# This rule aggregates the filtered MAGs and reference genomes in one directory
rule aggregate_refs:
    input:
        filtered_mags="../results/MAG_binning/bins/filtered_mags/",
        refs_isolates="../data/reference_assemblies/hb_bacteria/non_redundant/single_genomes/contigs/"
    output:
        all_refs=temp("../scratch_link/reference_genomes_redundant/")
    log:
        "logs/ref_db/aggregate.log"
    threads: 2
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1500,
        runtime= "00:30:00"
    shell:
        "mkdir -p {output}; "
        "cp {input.filtered_mags}/*.fa {output}; "
        "cp {input.refs_isolates}/*.fna {output}"

# This rule runs dRep for 95% dereplication
rule run_dRep:
    input:
        all_refs="../scratch_link/reference_genomes_redundant/"
    output:
        drep_out=directory("../results/reference_db/")
    log:"logs/ref_db/dRep.log"
    benchmark:"logs/ref_db/dRep.benchmark"
    threads: 25
    conda:
        "envs/drep_env.yaml"
    params:
        compl=75
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "04:00:00"
    shell:
        "dRep dereplicate {output.drep_out} -g {input.all_refs}/* --clusterAlg single --completeness {params.compl} -p {threads}"

# This rule aggreagtes dRep data with checkm and gtdb-TK data to obtain some summary tables
rule parse_dRep:
    input:
        checkm_filt="../results/MAG_binning/checkm_QC/summary/filtered_MAGs_stats.tsv",
        gtdb_dir="../results/MAG_binning/gtdbtk_classification/",
        clust_dir="../results/reference_db",
        mtdata= "../data/metadata/RefGenomes_isolates_mtdata.csv"
    output:
        clust_info="../results/reference_db_filtered/summary_data_tables/clust_info.tsv",
        clust_assign="../results/reference_db_filtered/summary_data_tables/clust_assign.tsv",
        to_delete="../results/reference_db_filtered/summary_data_tables/to_delete.tsv",
        clust_final="../results/reference_db_filtered/summary_data_tables/clust_filtered.tsv",
        clust_win_final="../results/reference_db_filtered/summary_data_tables/clust_filtered_winners.tsv"
    conda:
        "envs/base_R_env.yaml"
    log:
        "logs/ref_db/parse_drep.log"
    params:
        us_func="scripts/useful_func.R"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "00:30:00"
    script:
        "scripts/MAGs/parse_drep.R"

# This rule plots dendograms with taxonomic info for each cluster, also it deletes fro, the dreplicated database all cluster that have no classification at the genus level
rule filter_dRep_db:
    input:
        clust_assign="../results/reference_db_filtered/summary_data_tables/clust_assign.tsv",
        mtdata= "../data/metadata/RefGenomes_isolates_mtdata.csv",
        to_delete="../results/reference_db_filtered/summary_data_tables/to_delete.tsv",
        ref_db="../results/reference_db/"
    output:
        out_fgs=directory("../results/reference_db_filtered/figures/secondary_clusters_dendograms/"),
        filt_db=directory("../results/reference_db_filtered/dereplicated_genomes_filtered/")
    conda:
        "envs/drep_env.yaml"
    log:
        "logs/ref_db/filter_drep_db.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    script:
        "scripts/MAGs/analyze_drep.py"