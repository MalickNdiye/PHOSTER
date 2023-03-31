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

def get_all_vProt(wildcards):
    sample_table = checkpoints.split_vamb_contigs.get(**wildcards).output[3]
    IDS = []
    with open(sample_table,'r') as infile:
        for line in infile:
            line = line.rstrip()
            IDS.append(line)
    return(expand("../results/MAG_binning/vBins/phamb/sample_annotation/{vSam}/{vSam}.predicted_proteins.faa", vSam=IDS))


rule gene_2_genome:
    input:
        proteins = get_all_vProt
    output:
        all_vprot = "../results/Vcontact2/proteins/all_viral_prot.faa",
        gene_2_genome = "../results/Vcontact2/proteins/gene_to_genome.csv"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 8000,
        runtime= "00:10:00"
    shell:
        "cat {input} > {output.all_vprot}; "
        "python scripts/Viral_classification/gener2genome.py -p {output.all_vprot} -o {output.gene_2_genome} -s 'Prodigal-FAA'"


rule run_vcontact:
    input:
        all_vprot = "../results/Vcontact2/proteins/all_viral_prot.faa",
        gene_2_genome = "../results/Vcontact2/proteins/gene_to_genome.csv"
    output:
        directory("../results/Vcontact2/vCONTACT_results")
    threads: 48
    params: 
        condaenv="resources/conda_envs/vcontact2"
    log:
        "logs/vcontact/run_vcontact2.log"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "24:00:00"
    shell:
        "conda activate {params.condaenv}; "
        "vcontact2 -t {threads} --raw-proteins {input.all_vprot} --rel-mode 'Diamond' --proteins-fp {input.gene_2_genome} --db 'ProkaryoticViralRefSeq211-Merged' --pcs-mode MCL --vcs-mode ClusterONE --c1-bin {params.condaenv}/bin/cluster_one-1.0.jar --output-dir {output} -e 'cytoscape' -e 'csv'"