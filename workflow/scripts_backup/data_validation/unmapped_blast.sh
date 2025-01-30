infile=$1
indir=$(dirname ${1})
outdir=$(dirname ${2})

mkdir -p $outdir
echo $outdir

awk -F "|" '/^>/ {close(F); ID=$1; gsub("^>", "", ID); F=ID".fna"} {print >> F}' $infile
mv *.fna $outdir

for i in $outdir/*.fna; do
  echo "$i"
  blastn -db nt -query $i -max_target_seqs 1 -outfmt 6 -out ${i}.txt -remote
  cat ${i}.txt >> $2
done

rm ${outdir}/*.fna
#-num_threads ${3}
#-mt_mode 1
