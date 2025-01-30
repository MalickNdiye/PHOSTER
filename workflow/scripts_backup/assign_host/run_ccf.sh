# This script takes a filtered MAG and refrence genome and uses CrisprCasFinder to identiry CRISR spacers in the genomes

# set global variables
input=$1
output=$2

genome=$3

new_wd=${output%/*}/${genome}_run # CCF is a terrible software, it spits all intermediate files in the working directory, so I will launch the command in a directory for each single genome

# Create output directoris
mkdir -p ${output%/*}
mkdir -p ${new_wd}

cp  ${input}/${genome} ${new_wd} # The genome file must be moved in the working directory because the container cannot follow sym links

echo "input:"
echo ${input}/${genome}

echo "output:"
echo ${output}

echo "temporary WD:"
echo ${new_wd}

# create a copy of the container in the
cp resources/containers/CrisprCasFinder.simg ${new_wd}/CrisprCasFinder_${genome}.simg

#move to working directory
cd ${output%/*}/${genome}_run

# launch CCF
singularity exec -B $PWD --no-home CrisprCasFinder_${genome}.simg perl /usr/local/CRISPRCasFinder/CRISPRCasFinder.pl \
-so /usr/local/CRISPRCasFinder/sel392v2.so -cf /usr/local/CRISPRCasFinder/CasFinder-2.0.3 -cas \
-drpt /usr/local/CRISPRCasFinder/supplementary_files/repeatDirection.tsv \
-rpts /usr/local/CRISPRCasFinder/supplementary_files/Repeat_List.csv \
-out ./temp_${genome} -in ${genome}

# move results to final directory
mv ./temp_${genome} ${genome}-spacers
mv ${genome}-spacers ..

# Remove temporary Files
rm ${genome}
rm CrisprCasFinder_${genome}.simg
cd ..
rm -rf ${genome}_run