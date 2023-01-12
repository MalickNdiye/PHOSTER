rule mappings_stats:
    input:
        bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.bam",
    output:
        temp("../results/mapping/mapout_{sample}{type}/{sample}{type}_mapstats.tsv")
    wildcard_constraints:
        sample="\d+"
    conda:
        "envs/map_env.yaml"
    log:
        "logs/mapping/mapstats/{sample}{type}_mapstats.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 20000,
        runtime= "01:00:00"
    threads: 8
    shell:
        "samtools flagstat {input.bam} -@ {threads} > {output}"

rule merge_mapping_stats:
    input:
        filelist=expand("../results/mapping/mapout_{sample}{{type}}/{sample}{{type}}_mapstats.tsv", sample=config["sam_names"])
    output:
        "../results/mapping/all_{type}_mapstats.tsv"
    wildcard_constraints:
        sample="\d+"
    conda:
        "envs/base_R_env.yaml"
    log:
        "logs/mapping/{type}_merge_mapstats.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:15:00"
    threads: 2
    script:
        "scripts/mapping/format_flagstats.R"

rule sort_bam_QC:
    input:
        bam="../scratch_link/mapping/mapdata_{type}/{sample}{type}_mapping.bam"
    output:
        "../scratch_link/mapping/sorted_bams/{sample}{type}_mapping_sorted.bam"
    log:
        "logs/mapping/QC/sort_bamQC.log"
    conda:
        "envs/metapop_env.yaml"
    threads: 5
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "00:30:00
    script:
        "samtools sort {input.bam} -n -@ {threads} -o {output}"

rule prepare_bamQC:
    input:
        bam=expand("../scratch_link/mapping/sorted_bams/{sample}{{type}}_mapping_sorted.bam", sample=config["sam_names"])
    output:
        "../results/mapping/QC/bamqc_config.txt"
    log:
        "logs/mapping/QC/prepare_bamQC.log"
    threads: 2
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 1000,
        runtime= "00:15:00
    script:
        "scripts/mapping/prepre_bamQC.py"

rule bamQC:
    input:
        config="../results/mapping/QC/bamqc_config.txt"
    output:
        directory("../results/mapping/QC/multi_bamQC/")
    log:
        "logs/mapping/QC/multi_bamQC.log"
    benchmark:
        "logs/mapping/QC/multi_bamQC.benchmark"
    threads: 12
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00
    conda:
        "envs/qalimap.yaml"
    shell:
        "qualimap multi-bamqc -d {input.config} -outdir {output} -r -c -nt {threads}"
