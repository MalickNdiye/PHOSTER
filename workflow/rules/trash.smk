###############################################################################
################################## METAPOP #####################################
################################################################################

rule megre_libs:
    input:
        expand("../results/mapping/{sample}{type}_library_count.tsv", sample=config["sam_names"], type=acc)
    output:
        "../results/mapping/{type}_library_count.tsv"
    threads: 2
    log:
        "logs/mapping/concat_lib_{type}.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 2000,
        runtime= "01:00:00"
    shell:
        "cat {input} >> {output}"



# To make metapop run I modified the environment by using the script that are in the etapop github page instead of the ones in conda
rule run_metapop_virus:
    input:
        bams="../results/mapping/mapdata_virus/",
        cf="../results/mapping/virus_library_count.tsv",
        ref_bact="../data/reference_assemblies/viral_clusters_GB/"
    output:
        directory("../results/metapop/viruses/")
    conda:
        "envs/metapop_env.yaml"
    threads: 25
    params:
        lib="/work/FAC/FBM/DMF/pengel/beemicrophage/mndiaye1/PHOSTER/workflow/.snakemake/conda/6ee4f263563c6d0db030078388037b81/lib/R/library/"
    log:
        "logs/metapop/run_metapop.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 70000,
        runtime= "07:00:00"
    shell:
        "(metapop --input_samples {input.bams} --reference {input.ref_bact} --norm {input.cf} -l {params.lib} -o {output} --threads {threads})2> {log}"

# rule run_metapop_bacteria:
#     input:
#         bams="../results/mapping/mapdata_{type}/",
#         cf="../results/mapping/{type}_library_count.tsv",
#         ref_bact="../scratch_link/Database_play/dRep_out/dereplicated_genomes/"
#     output:
#         directory("../results/metapop/{type}/")
#     conda:
#         "envs/metapop_env.yaml"
#     threads: 46
#     params:
#         lib="/work/FAC/FBM/DMF/pengel/beemicrophage/mndiaye1/PHOSTER/workflow/.snakemake/conda/6ee4f263563c6d0db030078388037b81/lib/R/library/"
#     log:drep_db.fasta
#         "logs/metapop/run_metapop_{type}.log"
#     resources:
#         account = "pengel_beemicrophage",
#         mem_mb = 200000,
#         disk_mb= 10000,
#         runtime= "72:00:00"
#     shell:
#         "(metapop --input_samples {input.bams} --reference {input.ref_bact} --norm {input.cf} -l {params.lib} --whole_genomes -o {output} --threads {threads})2> {log}"

rule run_metapop_bacteria_test:
    input:
        bams="../results/mapping/mapdata_{type}_test/",
        cf="../results/mapping/{type}test_library_count.tsv",
        ref_bact="../scratch_link/Database_play/dRep_out/dereplicated_genomes/"
    output:
        directory("../results/metapop/{type}_test/")
    conda:
        "envs/metapop_env.yaml"
    wildcard_constraints:
        sample="\d+"
    threads: 25
    log:
        "logs/metapop/run_metapop_{type}test.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 200000,
        runtime= "72:00:00"
    shell:
        "(metapop --input_samples {input.bams} --reference {input.ref_bact} --norm {input.cf} --whole_genomes -o {output} --threads {threads})2> {log}"


################################ Unmapped reads assembly #######################
rule unmapped_metaspades_assembly:
    input:
        unmapped_R1="../results/mapping/mapout_{sample}{type}/{sample}{type}_R1_unmapped.fastq",
        unmapped_R2="../results/mapping/mapout_{sample}{type}/{sample}{type}_R2_unmapped.fastq"
    output:
        dir=directory("../results/assembly/unmapped_assembly/{sample}{type}_metaspades_unampped/"),
        file="../results/assembly/unmapped_assembly/{sample}{type}_metaspades_unampped/{sample}{type}_unmapped_contigs.fasta"
    threads: 25
    wildcard_constraints:
        sample="\d+"
    log:
        "logs/assembly/unmapped_assembly/{sample}{type}_unmapped_assembly.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 120000,
        runtime= "04:00:00"
    conda:
        "envs/assembly.yaml"
    shell:
        "spades.py --meta --pe1-1 {input.unmapped_R1} --pe1-2 {input.unmapped_R2} \
        -o {output.dir} \
        -k 21,33,55,77,99,127 -m 120 -t {threads}; "
        "mv {output.dir}/contigs.fasta {output.file}"

sample_nam="francesco"
# This rule creates both a file containing genereal info for all the assemblies and
# a fasta with all the contigs
rule parse_unmapped_assemblies:
    input:
        expand("../results/data_validation/reads_mapping/unmapped_asssembly/{sample}_metaspades_unampped/{sample}_unmapped_contigs.fasta", sample=sample_nam)
    output:
        "../results/data_validation/reads_mapping/unmapped_asssembly/assemblies_summary.txt",
        "../results/data_validation/reads_mapping/unmapped_asssembly/P_all_filtered_contigs.fasta"
    threads: 2
    log:
        "logs/data_validation/reads_mapping/parse_assemblies.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 4000,
        runtime= "00:30:00"
    script:
        "scripts/data_validation/assembly_data_parser.py"


# blast the unmapped contigs and return the best hit
rule blast_P_unmapped_contigs:
    input:
        "../results/data_validation/reads_mapping/unmapped_asssembly/P_all_filtered_contigs.fasta"
    output:
        "../results/data_validation/reads_mapping/unmapped_asssembly/blast/P_filtered_contigs_blastout.txt"
    log:
        "logs/data_validation/unmapped_assembly/blast_long_unmapped.log"
    conda:
        "envs/blast.yaml"
    threads: 25
    resources:
        mem_mb = 4000,
        runtime= "10:00:00",
        account = "pengel_beemicrophage"
    shell:
        "./scripts/data_validation/unmapped_blast.sh {input} {output} {threads}"


rule cd_hit_unmapped_assembly:
    input:
        "../results/data_validation/reads_mapping/unmapped_asssembly/P_all_filtered_contigs.fasta"
    output:
        clst_80="../results/data_validation/reads_mapping/unmapped_asssembly/clustering/unmapped_contigs_clstr_80.out",
        clst_100="../results/data_validation/reads_mapping/unmapped_asssembly/clustering/unmapped_contigs_clstr_100.out"
    log:
        "logs/data_validation/unmapped_assembly/unmapped_clustering.log"
    conda:
        "envs/cd-hit.yaml"
    threads: 4
    resources:
        account = "pengel_beemicrophage",
    shell:
        "cd-hit-est -d 0 -i {input} -o {output.clst_80} -c 0.8; "
        "cd-hit-est -d 0 -i {input} -o {output.clst_100} -c 1"

rule parse_unmapped_clusters:
    input:
        "../results/data_validation/reads_mapping/unmapped_asssembly/clustering/unmapped_contigs_clstr_80.out",
        "../results/data_validation/reads_mapping/unmapped_asssembly/clustering/unmapped_contigs_clstr_100.out"
    output:
        "../results/data_validation/reads_mapping/unmapped_asssembly/clustering/unmapped_contigs_clstr_parsed_80.txt",
        "../results/data_validation/reads_mapping/unmapped_asssembly/clustering/unmapped_contigs_clstr_parsed_100.txt"
    threads: 2
    conda:
        "envs/base_R_env.yaml"
    log:
        "logs/data_validation/reads_mapping/parse_assemblies.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 4000,
        runtime= "00:30:00"
    script:
        "scripts/data_validation/parse_cd_hit.R"


################################################################################
########################### ASSEMBLY ###########################################
################################################################################
rule merge_reads:
    input:
        R1="../data/trimmed_reads/{sample}_R1_paired.fastq.gz",
        R2="../data/trimmed_reads/{sample}_R2_paired.fastq.gz"
    output:
        dir=directory("../data/merged_reads/{sample}_merging_output/"),
        fw="../data/merged_reads/{sample}_merging_output/{sample}.unassembled.forward.fastq",
        rv="../data/merged_reads/{sample}_merging_output/{sample}.unassembled.reverse.fastq",
        merged="../data/merged_reads/{sample}_merging_output/{sample}.assembled.fastq"
    threads: 20
    log:
        "logs/assembly/{sample}_reads_merging.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 6000,
        runtime= "00:20:00"
    conda:
        "envs/trimmomatic.yaml"
    shell:
        "mkdir -p {output.dir}; "
        "pear -f {input.R1} -r {input.R2} -o {output.dir}/{wildcards.sample} -n 40 -j {threads} -v 1; "# is -v=1 good?
        "touch {output.fw} {output.rv} {output.merged}"


rule metaspades_assembly:
    input:
        R1="../data/merged_reads/{sample}_merging_output/{sample}.unassembled.forward.fastq",
        R2="../data/merged_reads/{sample}_merging_output/{sample}.unassembled.reverse.fastq",
        merged="../data/merged_reads/{sample}_merging_output/{sample}.assembled.fastq"
    output:
        dir=directory("../data/assemblies/metaspades/{sample}_metaspades_assembly/"),
        file="../data/assemblies/metaspades/{sample}_metaspades_assembly/{sample}_contigs.fasta"
    threads: 40
    log:
        "logs/assembly/metaspades/{sample}_metaspades_assembly.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 120000,
        runtime= "04:00:00"
    conda:
        "envs/assembly.yaml"
    shell:
        "spades.py --meta --merged {input.merged} \
        -s {input.R1} -s {input.R2} -o {output.dir} \
        -k 21,33,55,77,99,127 -m 120 -t {threads}; "
        "mv {output.dir}/contigs.fasta {output.file}"

rule quast_QC:
    input:
        contigs=expand("../data/assemblies/{assembler}/{sample}_{assembler}_assembly/{sample}_contigs.fasta", assembler=["metaspades"], sample=sample_nam),
        reads=expand("../data/merged_reads/{sample}_merging_output/{sample}.assembled.fastq", sample=sample_nam)
    output:
        directory("../results/assembly/quast/")
    threads: 25
    log:
        "logs/assembly/quast/quast_QC.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 50000,
        runtime= "03:00:00"
    conda:
        "envs/quast_env.yaml"
    shell:
        "metaquast.py {input.contigs} -o {output} \
         --contig-thresholds 1000,2000,2500,3000,5000,7500,10000,15000,20000,25000,30000,35000,40000,45000,50000,55000,60000,70000,80000,90000,100000 \
         -t {threads} --unique-mapping"
         #--mp12 {input.reads}

# rule get_unmapped:
#     input:
#         bam="../results/mapping/mapdata_{type}/{sample}{type}_mapping.bam"
#     output:
#         sorted_unmapped=".../results/mapping/mapout_{sample}{type}/{sample}{type}_unmapped.sorted.bam",
#         unmapped="../results/mapping/mapout_{sample}{type}/{sample}{type}_unmapped.bam",
#         unmapped_R1="../results/mapping/mapout_{sample}{type}/{sample}{type}_R1_unmapped.fastq",
#         unmapped_R2="../results/mapping/mapout_{sample}{type}/{sample}{type}_R2_unmapped.fastq"
#     wildcard_constraints:
#         sample="\d+"
#     conda:
#         "envs/manip_bap_env.yaml"
#     threads: 8
#     log:
#         "logs/mapping/unmapped/{sample}{type}_extract_unmapped.log"
#     resources:
#         account= "pengel_beemicrophage",
#         mem_mb= 20000,
#         runtime= "01:00:00"
#     shell:
#         "samtools view -b -f 4 -@ {threads} {input.bam} > {output.unmapped}; "
#         "samtools sort {output.unmapped} -n -@ {threads} -o {output.sorted_unmapped}; "
#         "bedtools bamtofastq -i {output.sorted_unmapped} -fq {output.unmapped_R1} -fq2 {output.unmapped_R2}"

# rule parse_HBdVir:
#     input:
#         assemblies= "../data/reference_assemblies/Phages_GB/HBvirDBv1_DB/Single_Samples/{vGB}.fasta",
#     output:
#         concat_assembly = "../results/assembly/concat_assembly/{vGB}_concat_assembly.fasta"
#     params:
#         length_t = 1000,
#         cov_t = 1
#     conda: "envs/mags_env.yaml"
#     resources:
#         account="pengel_beemicrophage",
#         mem_mb= 2000,
#         runtime= "00:30:00"
#     log:
#         "logs/assembly/HF/vGB/{vGB}_parse_assembly.log"
#     script:
#         "scripts/assembly/parse_reference_metagenome.py"

# rule run_deepvirfinder: discarded beacause it takes ages and desn't add many contigs, and those that it adds are weird
#     input:
#          "../results/assembly/concat_assembly/{sample}_concat_assembly.fasta"
#     output:
#         outdir =temp(directory("../scratch_link/viral_identification/deepvirfinder/{sample}_deepvirfinder")),
#         score = "../results/viral_identification/deepvirfinder/{sample}_deepvirfinder/{sample}_deepvirfinder_score.tsv"
#     log:
#         "logs/viral_identification/deepvirfinder/{sample}_dvfinder.log"
#     resources:
#         account="pengel_beemicrophage",
#         mem_mb= 100000,
#         runtime = "12:00:00"
#     threads: 20
#     params:
#         dvfdir="resources/softwares/DeepVirFinder"
#     conda:
#         "envs/dvf.yaml"
#     shell:
#         "export THEANO_FLAGS='base_compiledir={output.outdir}'; "
#         "python {params.dvfdir}/dvf.py -i {input} -o {output.outdir} -c {threads}; "
#         "dir=$(dirname {output.score}); mkdir -p ${{dir}}; "
#         "mv {output.outdir}/{wildcards.sample}_concat_assembly.fasta_gt1bp_dvfpred.txt {output.score}"

# rule dereplicate_viral_contigs_dRep:
#     input:
#         concat_cont="../scratch_link/viral_contigs/concat_single_genomes"
#     output:
#         derep_dir=directory("../results/assembly/viral_contigs/dereplicated/dRep"),
#         drep_fasta="../results/assembly/viral_contigs/dereplicated/all_viral_contigs_concat_drep.fasta"
#     resources:
#         account="pengel_beemicrophage",
#         runtime="24:00:00",
#         mem_mb = 500000
#     threads: 25
#     conda:
#         "envs/drep_env.yaml"
#     log:
#         "logs/polish_vMAGs/dereplication_dep.log"
#     benchmark:
#         "logs/polish_vMAGs/dereplication_dep.benchmark"
#     shell:
#         """
#         mkdir -p {output.derep_dir}/data_tables
#         touch {output.derep_dir}/data_tables/Bdb.csv
#         echo "genome,location" > {output.derep_dir}/data_tables/Bdb.csv

#         for file in {input}/*; do
#             bs=$(basename ${{file}})
#             abs=$(realpath ${{file}})

#             if [[ "${{bs}}" != .* ]]; then
#                 echo "${{bs}},${{abs}}" >> {output.derep_dir}/data_tables/Bdb.csv
#             fi
#         done

#         dRep dereplicate {output.derep_dir} -p {threads} --multiround_primary_clustering --greedy_secondary_clustering --S_algorithm fastANI -nc .85 -sa .97 -l 3000 -N50W 0 -sizeW 1 --ignoreGenomeQuality --clusterAlg single
#         find {output.derep_dir}/dereplicated_genomes/ -type f -exec cat {} \; > {output.drep_fasta}
#         """

rule parse_viral_contigs_dRep: # a few contigs are ignored by drep for some weird rson, let's put them back in the file
    input:
        derep_dir="../results/assembly/viral_contigs/dereplicated/dRep",
        concat_cont=expand("../results/assembly/viral_contigs/concat/{sample}_viral_contigs_concat.fasta", sample=config["samples"]),
        drep_fasta="../results/assembly/viral_contigs/dereplicated/all_viral_contigs_concat_drep.fasta"
    output:
        concat_cont_all=temp("../results/assembly/viral_contigs/concat/all_viral_contigs_concat_tmp.fasta"),
        drep_fasta_final="../results/assembly/viral_contigs/dereplicated/all_viral_contigs_concat_drep2.fasta",
        outtab="../results/assembly/viral_contigs/dereplicated/drep_name_changes.tsv"
    resources:
        account="pengel_beemicrophage",
        runtime="01:00:00",
        mem_mb = 20000
    conda:
        "envs/mOTUpan.yaml"
    log:
        "logs/polish_vMAGs/dereplication_dep_parsing.log"
    benchmark:
        "logs/polish_vMAGs/dereplication_dep_parsing.benchmark"
    shell:
        "cat {input.concat_cont} > {output.concat_cont_all}; "
        "python scripts/vMAGs/parse_drep_output.py -i {input.drep_fasta} -d {input.derep_dir} -c {output.concat_cont_all} -o {output.drep_fasta_final} -t "