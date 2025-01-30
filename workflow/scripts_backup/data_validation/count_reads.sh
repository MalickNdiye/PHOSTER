#!/bin/bash

## This script counts the reads and total bp before and after trimming and
## writes them into a file

# input files
pre_T=$1
post_T=$2
post_F=$3

# output file
outfile=$4
echo -e "sample\treads_preQC\tbases_preQC\treads_postQC\tbases_postQC\treads_postFilt\tbases_postFilt" > $outfile


echo ${pre_T}
fname_preT="$(basename -- ${pre_T})"
fname_preT=$(echo $fname_preT | cut -d'_' -f1)
echo $fname_preT

echo ${post_T}
fname_postT="$(basename -- ${post_T})"
fname_postT=$(echo $fname_postT | cut -d'_' -f1)
echo $fname_postT

echo ${post_F}
fname_postF="$(basename -- ${post_F})"
fname_postF=$(echo $fname_postF | cut -d'_' -f1)
echo $fname_postF

##counting reads before trimming and writing it togheter with sample name to a file
reads_bases_preT=$(zcat ${pre_T} |paste - - - -|cut -f2|wc -c -l|awk -v OFS="\n" '{print $0}')

##counting reads after trimming and appending it to the row corresponding to the sample
reads_bases_postT=$(zcat ${post_T} |paste - - - -|cut -f2|wc -c -l|awk -v OFS="\n" '{print $0}')

##counting reads after Filtering out the host and human and appending it to the row corresponding to the sample
reads_bases_postF=$(zcat ${post_F} |paste - - - -|cut -f2|wc -c -l|awk -v OFS="\n" '{print $0}')

#write into final file
echo -e "${fname_preT}\t${reads_bases_preT}\t${reads_bases_postT}\t${reads_bases_postF}" >> $outfile


# Trasform all spaces into tab (doesn't really work)
sed -E 's/[[:space:]]+/\t/g' $outfile > testfile.tmp && mv testfile.tmp $outfile
