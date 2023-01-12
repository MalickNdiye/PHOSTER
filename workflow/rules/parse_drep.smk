rule parse_dRep:
    input:
        checkm_filt="../results/MAG_binning/checkm_QC/summary/filtered_MAGs_stats.tsv",
        gtdb_dir="../results/MAG_binning/gtdbtk_classification/",
        clust_dir="../results/reference_db/data_tables/",
        mtdata= "../data/metadata/RefGenomes_isolates_mtdata.csv"
    output:
        clust_info="../results/reference_db/summary_data_tables/clust_info.tsv",
        clust_assign="../results/reference_db/summary_data_tables/clust_assign.tsv",
        to_delete="../results/reference_db/summary_data_tables/to_delete.tsv",
        clust_final="../results/reference_db/summary_data_tables/clust_filtered.tsv",
        clust_win_final="../results/reference_db/summary_data_tables/clust_filtered_winners.tsv"
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

rule filter_dRep_db:
    input:
        clust_assign="../results/reference_db/summary_data_tables/clust_assign.tsv",
        mtdata= "../data/metadata/RefGenomes_isolates_mtdata.csv",
        to_delete="../results/reference_db/summary_data_tables/to_delete.tsv",
        clust_dir="../results/reference_db/data/Clustering_files/",
        old_db="../results/reference_db/dereplicated_genomes/"
    output:
        out_fgs=directory("../results/reference_db/figures/secondary_clusters_dendograms/"),
        filt_db=directory("../results/reference_db/dereplicated_genomes_filtered/")
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
