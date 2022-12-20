rule aggregate_checkm_QC:
    input:
        bins="../results/MAG_binning/bins",
        stats=expand("../results/MAG_binning/checkm_QC/{sam_name2}B_checkm_QC/{sam_name2}B_checkm_QC_stats.tsv", sam_name2=config["sam_names"])
    output:
        full_stats="../results/MAG_binning/checkm_QC/summary/all_MAGs_stats.tsv",
        filtered_stats="../results/MAG_binning/checkm_QC/summary/filtered_MAGs_stats.tsv",
        filered_mags=directory("../results/MAG_binning/bins/filtered_mags")
    log:
        "logs/MAGs/checkm/aggregate_checkm.log"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    script:
        "scripts/MAGs/aggregate_checkm.py"
