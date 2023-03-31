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
        contigs="../results/MAG_binning/vBins/phamb/contigs.npz",
        length="../results/MAG_binning/vBins/phamb/contig_lengths.npz",
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
        hmmfile = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.hmmMiComplete105.tbl"
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
        hmmfile = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.hmmVOG.tbl"
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
        dvffile = "../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}_dvf/{vSam}.fna_gt2000bp_dvfpred.txt"
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