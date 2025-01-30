# working directory
wd=$1
echo "working directory: "
echo ${wd}

#inputs
phages=$2
spacersDB=$3
sample=$7

#outputs
blastout=$4
report=$5

#n
mkdir -p ${blastout%/*}

#CrisprOpenDB_path
OpenDB=${wd}/resources/default_DBs/CrisprOpenDB
OpenDB_tmp=${wd}/../scratch_link/${sample}_DB_tmp/

# copy content of CrisprOpenDB to OpenDB_tmp 
cp -r ${OpenDB} ${OpenDB_tmp}



cd ${OpenDB_tmp}
echo "CrisprOpenDB Directory:"
pwd -P

#echo "-------------- Installing CrisprOpenDB"
#python setup.py install
echo "-------------- identifying Hosts for viral contigs in ${phages} using ${spacersDB} as spacers database"
python CL_Interface.py -i ${wd}/${phages} -u -b ${wd}/${spacersDB}/mySpacersDB -m 2 -n $6 -t > CRISPRopenDB_${sample}.tmp

sed "s/['()]//g" CRISPRopenDB_${sample}.tmp | grep ',' | tr "," "\t" | sed 's/No hits found. Sorry/NoHitsFoundSorry/g' | sed 's/None/0/g' > ${wd}/${report} 

echo -e "Hit_nr,SPACER_ID,Query,identity,alignement_length,mismatch,gap,q_start,q_end,s_start,s_end,e_value,score,GENEBANK_ID,ORGANISM_NAME,SPECIES,GENUS,FAMILY,ORDER,SPACER,SPACER_LENGTH,COUNT_SPACER,POSITION_INSIDE_LOCUS,true_num_mismatch" > ${wd}/${blastout}
tail -n +2 -q *.csv >> ${wd}/${blastout} 

rm *.csv
rm CRISPRopenDB_${sample}.tmp

cd ${wd}/${blastout%/*}/
rm -rf ${OpenDB_tmp}

echo "DONE!"
