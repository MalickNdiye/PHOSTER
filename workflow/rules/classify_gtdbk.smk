rule classify_gtdbtk:
    input:
        filered_mags="../results/MAG_binning/bins/filtered_mags",
    output:
        class_out=directory("../results/MAG_binning/gtdbtk_classification/")
    log:
        "logs/MAGs/gtdbtk/gtdbtk_classification.log"
    benchmark:
        "logs/MAGs/gtdbtk/gtdbtk_classification.benchmark"
    threads: 8
    conda:
        "envs/gtdbk_env.yaml"
    params:
        db="resources/default_DBs/gtdbtk-2.1.1/db"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 150000,
        runtime= "03:00:00"
    script:
        "export GTDBTK_DATA_PATH={params.db}; "
        "gtdbtk classify_wf --genome_dir {input.filtered_mags} --extension fa --out_dir {output.class_out} --cpus {threads}"
