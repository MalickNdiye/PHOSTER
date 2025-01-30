#!/bin/bash

configfile=./popcogent_config.sh
source ${configfile}
source ~/miniconda3/bin/activate PopCOGenT
source ${mugsy_env}

echo "######################################################################"
echo "######################################################################"
echo "######################################################################"
echo -e "\n"


rm -rf ${alignment_dir}
rm  ${genome_dir}${genus}/*.renamed.*
export MUGSY_INSTALL=/users/mndiaye1/beemicrophage_dir/PHOSTER/workflow/resources/db2/mugsy/mugsy_x86-64-v1r2.3
echo ${MUGSY_INSTALL}
echo ${alignment_dir}

final_output_dir=$1
mkdir -p ${final_output_dir}
mkdir -p ${alignment_dir}${genus}
genus=$2

if [ "${slurm_str}" = "" ]
	then
        echo "Running PopCOGenT get_alignment_and_length_bias."
		echo "get_alignment_and_length_bias.py --genome_dir ${genome_dir}${genus} --genome_ext ${genome_ext} --alignment_dir ${alignment_dir}${genus} --mugsy_path ${mugsy_path} --mugsy_env ${mugsy_env} --base_name ${base_name} --final_output_dir ${final_output_dir} --num_threads ${num_threads} ${keep_alignments}"
		python get_alignment_and_length_bias.py --genome_dir ${genome_dir}${genus}/ --genome_ext ${genome_ext} --alignment_dir ${alignment_dir}${genus}/ --mugsy_path ${mugsy_path} --mugsy_env ${mugsy_env} --base_name ${base_name} --final_output_dir ${final_output_dir} --num_threads ${num_threads} ${keep_alignments}
		echo "Running PopCOGenT cluster."
        python cluster.py --base_name ${base_name} --length_bias_file ${final_output_dir}/${base_name}.length_bias.txt --output_directory ${final_output_dir} --infomap_path ${infomap_path} ${single_cell}
fi