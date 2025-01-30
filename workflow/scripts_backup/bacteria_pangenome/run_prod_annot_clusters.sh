outdir_fna=$1
outdir_faa=$2
clust_genome=$3

# this

mkdir -p $1
mkdir -p $2

echo ${clust_genome}
echo ${clust_genome}/*


for fa in ${clust_genome}/*; do
  echo ${fa}
  filename=$(basename ${fa})
  genome=${filename%.fa}
  genome=${genome%.fna}
  echo -e "\n\t-----annotating ${fa} ..."
  prodigal -i ${fa} -d ${outdir_fna}/${genome}_genes.fna -a ${outdir_faa}/${genome}_genes.faa -p meta
done
