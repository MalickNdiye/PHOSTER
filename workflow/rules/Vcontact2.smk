####################################################################################
################################## vContact2 #######################################
####################################################################################
# The foilowing rules are used to prepare the viral contigs/bins file to run vContact2
################################################################################
################################################################################

# this rule split the contigs into single genomes fasta files
rule split_contig:
    input:
        concat=expand("../results/assembly/viral_contigs/concat/{sample}_viral_contigs_concat.fasta", sample=config["samples"]),
        filt=expand("../results/assembly/viral_contigs/filtered/{sample}_viral_contigs_filtered.fasta", sample=config["samples"])
    output:
        concat_dir=temp(directory("../scratch_link/viral_contigs/concat_single_genomes/")),
        filt_dir=temp(directory("../scratch_link/viral_contigs/filtred_single_genomes/"))
    resources:
        account="pengel_beemicrophage",
        runtime="0:30:00",
        mem_mb = 8000
    threads: 1
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/vMAGs/split_contigs.log"
    shell:
        """
        for file in {input.concat} ; do
            python scripts/assembly/split_assembly.py -f $file -d {output.concat_dir}
        done

        for file in {input.filt} ; do
            python scripts/assembly/split_assembly.py -f $file -d {output.filt_dir}
        done
        """

rule prodigal_concat_bins:
    input:
        fasta=expand("../results/assembly/viral_contigs/concat/{sample}_viral_contigs_concat.fasta", sample=config["samples"])
    output:
        agg_fasta=temp("../results/assembly/viral_contigs/concat/all_viral_contigs_concat.fasta"),
        proteins = "../results/assembly/viral_contigs/annotations/all_viral_contigs_concat.faa",
        genes = "../results/assembly/viral_contigs/annotations/all_viral_contigs_concat.fna"
    threads: 1
    conda:
        "envs/prodigal.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "01:30:00"
    log:
        "logs/vMAGs/annotations/prodigal.log"
    shell:
        "cat {input.fasta} > {output.agg_fasta}; "
        "prodigal -i {output.agg_fasta} -d {output.genes} -a {output.proteins} -p meta -g 11"

rule gene_2_genome:
    input:
        all_vprot = "../results/assembly/viral_contigs/annotations/all_viral_contigs_concat.faa"
    output:
        gene_2_genome = "../results/Vcontact2/gene_to_genome.csv"
    threads: 1
    conda:
        "envs/mags_env.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "00:10:00"
    log:
        "logs/vcontact/gene_to_genome.log"
    shell:
        "python scripts/Viral_classification/gene2genome.py -p {input.all_vprot} -o {output.gene_2_genome} -s 'Prodigal-FAA'"

rule run_vcontact:
    input:
        all_vprot = "../results/assembly/viral_contigs/annotations/all_viral_contigs_concat.faa",
        gene_2_genome = "../results/Vcontact2/gene_to_genome.csv"
    output:
        directory("../results/Vcontact2/vCONTACT_results")
    threads: 48
    params:
        condaenv="resources/conda_envs/vcontact2"
    log:
        "logs/vcontact/run_vcontact2.log"
    benchmark:
        "logs/vcontact/run_vcontact2.benchmark"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 512000,
        runtime= "15:00:00"
    shell:
        """
        bash -c '. $HOME/.bashrc
            conda activate {params.condaenv}
            vcontact2 -t {threads} --raw-proteins {input.all_vprot} --rel-mode 'Diamond' --proteins-fp {input.gene_2_genome} --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {params.condaenv}/bin/cluster_one-1.0.jar --output-dir {output} -e 'cytoscape' -e 'csv''
        """
