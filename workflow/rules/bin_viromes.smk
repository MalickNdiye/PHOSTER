rule concat_viral_assemblies:
    input:
        refs="../data/reference_assemblies/Phages_GB/HBvirDBv1_DB/Single_Samples"
    output:
        concat="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna",
        tab="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_tab.ref"
    conda:
        "envs/mOTUpan.yaml"
    threads: 1
    log:
        "logs/vMAGs/backmapping/concat_genomes/HBvirDBv1_build_bowtie_index.log"
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 20000,
        runtime= "01:00:00"
    shell:
        "python scripts/vMAGs/concat_viral_assemblies.py {output.concat} {output.tab} {input.refs}/*.fasta"


rule build_backmapping_viruses_index:
    input:
        ref="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna"
    output:
        dir=directory("../results/mapping/HBvirDBv1_DB_index")
    conda:
        "envs/map_env.yaml"
    threads: 15
    params:
        basename="HBvirDBv1_DB",
    log:
        "logs/vMAGs/backmapping/indexing/HBvirDBv1_build_bowtie_index.log"
    resources:
        account= "pengel_beemicrophage",
        mem_mb= 20000,
        runtime= "01:00:00"
    shell:
        "mkdir -p {output.dir}; "
        "bowtie2-build {input.ref} {output.dir}/{params.basename} --threads {threads}"


rule backmapping_viral_db: # risk of running out of buffer memory, run fewer jobs at the time (with --jobs 50 works; maybe one could increase a bit)
    input:
        assembly = "../results/mapping/HBvirDBv1_DB_index", # we bin only bacteria
        R1 = "../data/host_filtered_reads/{sample}_R1_HF.fastq.gz",
        R2 = "../data/host_filtered_reads/{sample}_R2_HF.fastq.gz"
    output:
        sam=temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1.sam")
    resources:
        account="pengel_beemicrophage",
        runtime="10:00:00",
        mem_mb = 10000
    params:
        basename="HBvirDBv1_DB"
    threads: 15
    conda: "envs/map_env.yaml"
    log:
        "logs/vMAGs/backmapping/map/{sample}_mapped_to_HBvirDBv1.log"
    benchmark:
        "logs/MAGs/backmapping/map/{sample}_mapped_to_HBvirDBv1.benchmark"
    shell:
        "bowtie2 -x {input.assembly}/{params.basename} -1 {input.R1} -2 {input.R2} -S {output.sam} --threads {threads}"


rule sort_virome_backmapping_bam:
    input:
        sam="../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1.sam"
    output:
        bam= temp("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1.bam"),
        sorted_bam="../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1_sorted.bam",
        bai="../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1_sorted.bam.bai"
    log:
        "logs/vMAGs/backmapping/soringbams/sort_bam_{sample}.log"
    conda:
        "envs/map_env.yaml"
    threads: 5
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "01:30:00"
    shell:
        "samtools view -bh {input.sam} > {output.bam};"
        "samtools sort {output.bam} -@ {threads} -o {output.sorted_bam}; "
        "samtools index {output.sorted_bam} {output.bai} -@ {threads}"


rule run_vamb:
    input:
        ref="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna",
        sorted_bam=expand("../scratch_link/vMAG_binning/backmapping/{sample}/{sample}_mapped_to_HBvirDBv1_sorted.bam", sample=config["samples"])
    output:
        directory("../results/MAG_binning/vBins/HBvirDBv1/")
    log:
        "logs/vMAGs/binning.log"
    conda:
        "envs/vBinning_env.yaml"
    threads: 25
    params:
        min_contig_l=2000
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 500000,
        runtime= "4:00:00"
    shell:
        "vamb --outdir {output} --fasta {input.ref} --bamfiles {input.sorted_bam} -m {params.min_contig_l} -p {threads} -o C"



#######################################PHAMB############################################################################################
MICOMPLETEDB = "resources/default_DBs/phamb/Bact105.hmm"
VOGDB= "resources/default_DBs/phamb/AllVOG.hmm"
DVFDIR = "resources/softwares/DeepVirFinder"
TMP_DIR = "../scratch_link/mag_annotation_tmp"
CONTIGSUFFIX = '.fna'
ASSEMBLY = "../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna"


checkpoint split_vamb_contigs:
    input:
        contigs = "../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna"
    output:
        asmbl= directory("../results/MAG_binning/vBins/phamb/assembly/"),
        contigs=temp("../results/MAG_binning/vBins/phamb/contigs.npz"),
        length=temp("../results/MAG_binning/vBins/phamb/contig_lengths.npz"),
        sam_tab="../results/MAG_binning/vBins/phamb/sample_table.txt"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "00:30:00"
    params:
        out="../results/MAG_binning/vBins/phamb/",
        asmbl= "assembly",
        contigs="contigs.npz",
        length="contig_lengths.npz",
        sam_tab="sample_table.txt"
    threads: 1
    conda:
        "envs/vBinning_env.yaml"
    log:
        "logs/magannotation_log/split_contigs.log"
    shell:
        "split_contigs.py -c {input.contigs}; "
        "mkdir -p {params.out}; "
        "mv {params.asmbl} {params.contigs} {params.length} {params.sam_tab} {params.out}"

rule PRODIGAL:
    input:
        sam_tab= "../results/MAG_binning/vBins/phamb/sample_table.txt"
    output:
        proteins = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.predicted_proteins.faa",
        genes = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.predicted_proteins.fna"
    params:
        tmp_contigs = "sample_annotation/{vSam}.unzipped_contigs.fna",
        contigs = "../results/MAG_binning/vBins/phamb/assembly/{vSam}/{vSam}" + CONTIGSUFFIX
    threads: 1
    conda:
        "envs/prodigal.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "00:30:00"
    log:
        "logs/magannotation_log/prodigal/{vSam}.prodigal.log"
    shell:
        """
        if [[ {params.contigs} = *.gz ]]; then
            gunzip -c {params.contigs} > {params.tmp_contigs}
            prodigal -i {params.tmp_contigs} -d {output.genes} -a {output.proteins} -p meta -g 11 -q 2>{log}
            rm {params.tmp_contigs}
        else
            prodigal -i {params.contigs} -d {output.genes} -a {output.proteins} -p meta -g 11 -q 2>{log}
        fi
        """

rule miComplete:
    input:
        proteins = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.predicted_proteins.faa",
        DB = MICOMPLETEDB
    output:
        hmmfile = temp("../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.hmmMiComplete105.tbl")
    params:
        tmpoutput = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.micomplete.tmp"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "00:30:00"
    threads: 5
    conda:
        "envs/hmmer.yaml"
    log:
        "logs/magannotation_log/hmm/{vSam}.micomplete.log"
    shell:
        "hmmsearch --cpu {threads} -E 1.0e-05 -o {params.tmpoutput} --tblout {output.hmmfile} {input.DB} {input.proteins} 2>{log}"

rule VOG:
    input:
        proteins = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.predicted_proteins.faa",
        DB = VOGDB
    output:
        hmmfile = temp("../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.hmmVOG.tbl")
    params:
        tmpoutput = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.hmmVOG.tmp"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "00:30:00"
    threads: 5
    conda:
        "envs/hmmer.yaml"
    log:
        "logs/magannotation_log/hmm/{vSam}.VOG.log"
    shell:
        "hmmsearch --cpu {threads} -E 1.0e-05 -o {params.tmpoutput} --tblout {output.hmmfile} {input.DB} {input.proteins} 2>{log}"


rule DeepVirFinder:
    input:
        sam_tab= "../results/MAG_binning/vBins/phamb/sample_table.txt"
    output:
        dvffile = temp("../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}_dvf/{vSam}.fna_gt2000bp_dvfpred.txt")
    params:
        dvfscript = os.path.join(DVFDIR,'dvf.py'),
        contigs = "../results/MAG_binning/vBins/phamb/assembly/{vSam}/{vSam}" + CONTIGSUFFIX,
        dvfdir = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}_dvf"
    threads: 5
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    conda:
        "envs/dvf.yaml"
    log:
        "logs/magannotation_log/DVF/{vSam}.dvf.log"
    shell:
        """
        export THEANO_FLAGS="base_compiledir={params.dvfdir}"
        python {params.dvfscript} -i {params.contigs} -o {params.dvfdir} -l 2000 -c {threads} 2>{log}
        rm -r {params.dvfdir}/compiledir_Linux*
        """

rule aggergate_anno:
    input:
        anno=get_vids
    output:
        out_anno = directory("../results/MAG_binning/vBins/phamb/all_annotations")
    params:
        anno="../results/MAG_binning/vBins/phamb/sample_annotation"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    conda:
        "envs/vBinning_env.yaml"
    log:
        "logs/magannotation_log/QC/RF.log"
    shell:
        "mkdir -p {output.out_anno}; "
        "cat {params.anno}/*/*hmmMiComplete105.tbl > {output.out_anno}/all.hmmMiComplete105.tbl; "
        "cat {params.anno}/*/*hmmVOG.tbl > {output.out_anno}/all.hmmVOG.tbl; "
        "cat {params.anno}/*/*_dvf/*dvfpred.txt > {output.out_anno}/DVF.predictions.txt; "
        "head -n1 {output.out_anno}/DVF.predictions.txt > {output.out_anno}/DVF.header; "
        "grep -v 'pvalue' {output.out_anno}/DVF.predictions.txt > {output.out_anno}/DVF.predictions; "
        "cat {output.out_anno}/DVF.header {output.out_anno}/DVF.predictions > {output.out_anno}/all.DVF.predictions.txt"


rule RF_QC: # TODO I had to copy the phamb git (https://github.com/RasmussenLab/phamb.git) dbs directory in .snakemake/conda/6f50c953011b54affdd78f02f83f847d_/bin/ to make it work
    input:
        contigs="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna",
        vamb_dir="../results/MAG_binning/vBins/HBvirDBv1",
        anno="../results/MAG_binning/vBins/phamb/all_annotations"
    output:
        out_rf = directory("../results/MAG_binning/vBins/phamb/QC/")
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    conda:
        "envs/vBinning_env.yaml"
    log:
        "logs/magannotation_log/QC/RF.log"
    shell:
        "run_RF.py {input.contigs} {input.vamb_dir}/vae_clusters.tsv {input.anno} {output.out_rf}"


rule vamb_parser:
    input:
        contig_tab="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_tab.ref",
        all_anno = "../results/MAG_binning/vBins/phamb/all_annotations",
        qc="../results/MAG_binning/vBins/phamb/QC",
        vae="../results/MAG_binning/vBins/HBvirDBv1"
    output:
        vamb_clst = "../results/MAG_binning/vBins/vamb_clst_metadata.tsv"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    conda:
        "envs/base_R_env.yaml"
    log:
        "logs/vMAGs/parse_vamb.log"
    script:
        "scripts/vMAGs/vamb_parser.R"


rule create_bins:
    input:
        tab="../results/MAG_binning/vBins/vamb_clst_metadata.tsv",
        fasta="../results/MAG_binning/vRef/HBvirDBv1_DB_concat_vamb.fna"
    output:
        directory("../results/MAG_binning/vBins/bins")
    threads: 1
    conda: "envs/mOTUpan.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    log:
        "logs/vMAGs/create_bins.log"
    shell:
        "python scripts/vMAGs/vambinning.py {input.tab} {input.fasta} {output}; "
        "cat {output}/single_bins/* > {output}/all_viral_bins.fasta"

########################################### DeReplicate Viruses ################################################################################################

rule vRhyme_derep:
    input:
        bins="../results/MAG_binning/vBins/bins"
    output:
        directory("../results/MAG_binning/vBins/dereplicated_bins")
    threads: 20
    params:
        conda="resources/conda_envs/vrhyme"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "01:30:00"
    log:
        "logs/vMAGs/dereplication/derep.log"
    shell:
        """
        bash -c '. $HOME/.bashrc
            conda activate {params.conda}
            vRhyme -i {input}/all_viral_bins.fasta --derep_only --method composite --derep_id 0.97 --frac 0.50  --model NN -t {threads} -o {output}'
        """

rule filter_derep:
    input:
        bin="../results/MAG_binning/vBins/dereplicated_bins"
    output:
        fasta="../results/MAG_binning/vRef/filtered_ref/filtered_derep_bins.fasta",
        filter_tab="../results/MAG_binning/vRef/filtered_ref/filtering_tab.tsv",
        derep_tab="../results/MAG_binning/vRef/filtered_ref/derep_tab.tsv",
        singles=directory("../results/MAG_binning/vRef/filtered_ref/single_genomes")
    threads: 1
    conda: "envs/mOTUpan.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:15:00"
    log:
        "logs/vMAGs/dereplication/filter_derep.log"
    shell:
        "python scripts/vMAGs/filter_vRhyme_derep.py {input.bin}/vRhyme_dereplication/vRhyme_derep_composite_all_viral_bins.fa {input.bin}/vRhyme_dereplication/vRhyme_derep_composited-seqs_all_viral_bins.tsv {output.fasta} {output.filter_tab} {output.derep_tab} {output.singles}"
