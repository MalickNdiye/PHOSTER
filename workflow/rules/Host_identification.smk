rule find_CRISPR_spacers:
    input:
        all_refs="../scratch_link/reference_genomes_redundant/"
    output:
        spacers_dir=temp(directory("../results/spacers_db/hb_spacers/"))
    log:
        "logs/spacers_db/find_spacers.log"
    benchmark:
        "logs/spacers_db/find_spacers.benchmark"
    threads: 1
    conda:
        "envs/singularity.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "04:00:00"
    shell:
        "for i in {input}/*; "
        "do bs=$(basename ${i}); "
        "gen=${{i}}-spacers; "
        "singularity exec -B $PWD resources/containers/CrisprCasFinder.simg perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl \
        -so /usr/local/CRISPRCasFinder/sel392v2.so -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 \
        -drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv \
        -rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv \
        -out {output.spacers_dir}/${{gen}} -in ${{i}}; done"

rule parse_CCF:
    input:
        spacers_dir="../results/spacers_db/hb_spacers/"
    output:
        parsed_crispr=directory("../results/spacers_db/hb_spacers_parsed/")
    log:
        "logs/spacers_db/parse_CFFF.log"
    threads: 1
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "00:30:00"
    script:
        "scripts/assign_host/CFF_Parser.py"

rule create_spacersDB:
    input:
        openDB_spacers="resources/default_DBs/CrisprOpenDB/CrisprOpenDB/SpacersDB/SpacersDB.fasta",
        my_spacers="../results/spacers_db/hb_spacers_parsed"
    output:
        spacers_db=directory("../results/spacers_db/spacersDB/")
    log:
        "logs/spacers_db/create_spacersDB.log"
    benchmark:
        "logs/spacers_db/create_spacersDB.benchmark"
    threads: 1
    conda:
        "envs/CrisprOpenDB.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 10000,
        runtime= "01:00:00"
    shell:
        "cat {input.my_spacers}/CRISPR_v1.fasta {openDB_spacers} > {output.spacers_db}/myspacersDB.fasta",
        "makeblastdb -in {output.spacers_db}/myspacersDB.fasta -dbtype nucl -out {output.spacers_db}/mySpacersDB"

rule assign_host:
    input:
        phageDB="../data/reference_assemblies/viral_clusters_GB/HBvirDBv1_DB_derep_c99_aS95.fasta",
        spacers_db="../results/spacers_db/spacersDB/
    output:
        blastout="../results/host_assigniation/spacers_hit_blastout.txt",
        spacers_report="../results/host_assigniation/CRISPRopenDB_report.txt"
    log:
        "logs/assign_host/assign_host.log"
    benchmark:
        "logs/assign_host/assign_host.benchmark"
    threads: 30
    conda:
        "envs/CrisprOpenDB.yaml"
    resources:
        account = "pengel_beemicrophage",
        mem_mb = 100000,
        runtime= "07:00:00"
    shell:
        "scripts/assign_host/run_assign_host.sh $PWD {input.phageDB} {input.spacers_db} {output.blastout} {output.spacers_report} {threads}"
