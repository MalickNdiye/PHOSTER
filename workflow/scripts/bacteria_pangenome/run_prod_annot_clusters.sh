outdir_fna=$1
outdir_faa=$2
clust_genome=$3



cat ${clust_genome} | while read fa
do
  filename=$(basename ${fa})
  genome=${filename%.fa}
  genome=${genome%.fna}
  echo -e "\n\t-----annotating ${fa}..."
  prodigal -i ${fa} -d ${outdir_fna}/${genome}_genes.fna -a ${outdir_faa}/${genome}_genes.faa
done
