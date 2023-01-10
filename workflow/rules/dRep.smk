rule aggregate_refs:
    input:
        filtered_mags="../scratch_link/MAG_binning/bins/filtered_mags/",
        refs_isolates="../data/reference_assemblies/hb_bacteria/non_redundant/single_genomes/contigs/"
    output:
        all_refs=directory("../scratch_link/reference_genomes_redundant/")
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
        checkm_db="resources/default_DBs/checkm_db",
        compl=75
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "04:00:00"
    shell:
        "dRep dereplicate {output.drep_out} -g {input.all_refs}/* --clusterAlg single --completeness {params.compl} -p {threads}"
