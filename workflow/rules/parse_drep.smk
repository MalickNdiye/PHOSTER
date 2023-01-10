rule parse_dRep:
    input:
        checkm_filt="../results/MAG_binning/checkm_QC/summary/filtered_MAGs_stats.tsv",
        gtdb_dir="../results/MAG_binning/gtdbtk_classification/",
        clust_dir="../results/reference_db/data_tables/",
        mtdata= "../data/metadata/RefGenomes_isolates_mtdata.csv"
    output:
        clust_info="../results/reference_db/data_tables/clust_info.tsv",
        clust_assign="../results/reference_db/data_tables/clust_assign.tsv",
        to_delete="../results/reference_db/data_tables/clust_assign.tsv",
        clust_final="../results/reference_db/data_tables/clust_filtered.tsv",
        clust_win_final="../results/reference_db/data_tables/clust_filtered_winners.tsv"
    conda:
        "envs/base_R_env"
    log:
        "logs/ref_db/parse_drep.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 5000,
        runtime= "00:30:00"
    script:
        "scripts/MAGs/parse_drep.R"
