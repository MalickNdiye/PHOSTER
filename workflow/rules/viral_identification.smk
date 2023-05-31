################################################################################
############################## Phage identification  ###########################
################################################################################
#This part of the pipeline runs four tools on all the assembly, to identify
#phages, these sequences will then dereplicated and binned. Each phages id tools
#ouput are parsed and merged, a confidence score is attributed to each id
#and a fasta file of the phages sequences is extracted :
#   (1): PHAGES id tools
#   (1.1) concat assembly
#   (1.2) run phage id tools (virstorter, viralverify, deepvirfinder, vibrant)
#   (2) : Parsing outputs
#   (3): Run ChecV
################################################################################

##################################### (1.1) ######################################

rule concat_assembly:
    input:
        assemblies= expand("../results/assembly/HF_assembly/{assembler}/{{sample}}_{assembler}/{{sample}}_contigs_parsed.fasta", assembler=assemblers),
    output:
        concat_assembly = temp("../results/assembly/concat_assembly/{sample}_concat_assembly.fasta")
    log:
        "logs/assembly/concat_assembly/{sample}_concat.log"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 1000,
        runtime = "00:15:00"
    threads:1
    shell:
        "cat  {input} > {output.concat_assembly}"

##################################### (1.2) ######################################

rule run_virsorter:
    input:
        "../results/assembly/concat_assembly/{sample}_concat_assembly.fasta"
    output:
        outdir = temp(directory("../scratch_link/viral_identification/virsorter/{sample}_virsorter/")),
        score = "../results/viral_identification/virsorter/{sample}_virsorter/{sample}_virsorter_score.tsv",
        prophages="../results/viral_identification/virsorter/{sample}_virsorter/{sample}_virsorter_boundaries.tsv"
    params:
        db = directory("resources/default_DBs/virsorter_db"),
        container="resources/containers/virsorter2.sif"
    log:
        "logs/viral_identification/virsorter/{sample}_virsorter.log"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "04:00:00"
    threads: 20
    conda:
        "envs/singularity.yaml"
    shell:
        "org_dir=$PWD; "
        "infile=$(basename {input}); "
        "mkdir -p {output.outdir}; "
        "cp {params.container} {output.outdir}; cp {input} {output.outdir}; "
        "cd {output.outdir}; "
        "singularity run -B $PWD virsorter2.sif run -w ./output -i ${{infile}} --keep-original-seq -j {threads};"
        "rm virsorter2.sif ${{infile}}; "
        "cd ${{org_dir}}; "
        "dir=$(dirname {output.score}); mkdir -p ${{dir}}; "
        "mv {output.outdir}/output/final-viral-score.tsv {output.score}; "
        "mv {output.outdir}/output/final-viral-boundary.tsv {output.prophages}"

rule run_viralverify:
    input:
        "../results/assembly/concat_assembly/{sample}_concat_assembly.fasta"
    output:
        outdir = temp(directory("../scratch_link/viral_identification/viralverify/{sample}_viralverify/")),
        score = "../results/viral_identification/viralverify/{sample}_viralverify/{sample}_viralverify_score.csv"
    params:
        hmm = "resources/default_DBs/nbc_hmms.hmm"
    log:
        "logs/viral_identification/viralverify/{sample}_viralverify.log"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "04:0:00"
    threads:20
    conda:
        "envs/viralverify.yaml"
    shell:
        "viralverify -f {input} -o {output.outdir} --hmm {params.hmm} -t {threads};"
        "dir=$(dirname {output.score}); mkdir -p ${{dir}}; "
        "mv {output.outdir}/{wildcards.sample}_concat_assembly_result_table.csv {output.score}"

rule run_vibrant:
    input:
        assembly = "../results/assembly/concat_assembly/{sample}_concat_assembly.fasta"
    output:
        outdir = temp(directory("../scratch_link/viral_identification/vibrant/{sample}/")),
        score = "../results/viral_identification/vibrant/{sample}_vibrant/{sample}_vibrant_score.tsv",
        prophages="../results/viral_identification/vibrant/{sample}_vibrant/{sample}_vibrant_prophages.tsv"
    log:
        "logs/viral_identification/vibrant/{sample}_vibrant.log"
    params:
        db="resources/default_DBs/databases"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "04:00:00"
    threads:20
    conda:
        "envs/vibrant.yaml"
    shell:
        "VIBRANT_run.py -i {input.assembly} -folder {output.outdir} -d {params.db} -t {threads};"
        "dir=$(dirname {output.score}); mkdir -p ${{dir}}; "
        "mv {output.outdir}/VIBRANT_{wildcards.sample}_concat_assembly/VIBRANT_results_{wildcards.sample}_concat_assembly/VIBRANT_genome_quality_{wildcards.sample}_concat_assembly.tsv {output.score}; "
        "mv {output.outdir}/VIBRANT_{wildcards.sample}_concat_assembly/VIBRANT_results_{wildcards.sample}_concat_assembly/VIBRANT_integrated_prophage_coordinates_{wildcards.sample}_concat_assembly.tsv {output.prophages};"
        "cp {output.outdir}/VIBRANT_{wildcards.sample}_concat_assembly/VIBRANT_results_{wildcards.sample}_concat_assembly/* ${{dir}}"

##################################### (2) ######################################
rule parse_VI:
    input:
        vib_score = expand("../results/viral_identification/vibrant/{sample}_vibrant/{sample}_vibrant_score.tsv", sample=config["samples"]),
        vib_prophages=expand("../results/viral_identification/vibrant/{sample}_vibrant/{sample}_vibrant_prophages.tsv", sample=config["samples"]),
        vv_score = expand("../results/viral_identification/viralverify/{sample}_viralverify/{sample}_viralverify_score.csv", sample=config["samples"]),
        vs_score = expand("../results/viral_identification/virsorter/{sample}_virsorter/{sample}_virsorter_score.tsv", sample=config["samples"]),
        vs_prophages= expand("../results/viral_identification/virsorter/{sample}_virsorter/{sample}_virsorter_boundaries.tsv", sample=config["samples"])
    output:
        all_scores ="../results/viral_identification/ViralIdentification_scores.tsv",
        vib_scores="../results/viral_identification/vibrant/vibrant_all_scores.tsv",
        vv_scores ="../results/viral_identification/viralverify/viralverify_all_scores.tsv",
        vs_scores ="../results/viral_identification/virsorter/virsorter_all_scores.tsv",
        tab="../results/viral_identification/identification_tools_tab.tsv"
    log:
        "logs/viral_identification/parse_VI.log"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 10000,
        runtime = "00:10:00"
    threads:1
    conda:
        "envs/base_R_env.yaml"
    script:
        "scripts/viral_identification/parse_viral_identification.R"

rule extract_viral_contigs:
    input:
        all_scores ="../results/viral_identification/ViralIdentification_scores.tsv",
        concat_assembly = "../results/assembly/concat_assembly/{sample}_concat_assembly.fasta"
    output:
       temp("../results/assembly/viral_contigs/single_sample/{sample}_viral_contigs.fasta")
    log:
        "logs/viral_identification/extract_viruses/{sample}_viral_contigs.log"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 10000,
        runtime = "00:20:00"
    threads:1
    conda:
        "envs/mOTUpan.yaml"
    shell:
        "python scripts/viral_identification/extract_viruses_from_assembly.py {input.concat_assembly} {input.all_scores} {output}"

##################################### (3) ######################################
rule run_checkv:
    input:
        assembly = "../results/assembly/viral_contigs/single_sample/{sample}_viral_contigs.fasta"
    output:
        trimmed="../results/assembly/viral_contigs/single_sample/{sample}_viral_contigs_trimmed.fasta",
        tab="../results/viral_identification/checkv/{sample}_checkv_quality_summary.tsv"
    log:
        "logs/viral_identification/checkv/{sample}_checkv.log"
    params:
        outdir = directory("../scratch_link/viral_identification/checkv/{sample}/"),
        db="resources/default_DBs/checkv-db-v1.5/"
    resources:
        account="pengel_beemicrophage",
        mem_mb= 100000,
        runtime = "04:00:00"
    threads:20
    conda:
        "envs/checkv.yaml"
    shell:
        "mkdir -p {params.outdir}; "
        "checkv end_to_end {input.assembly} {params.outdir} -t {threads} -d {params.db}; "
        "mv {params.outdir}/quality_summary.tsv {output.tab}; "
        "cat {params.outdir}/viruses.fna {params.outdir}/proviruses.fna > {output.trimmed}"