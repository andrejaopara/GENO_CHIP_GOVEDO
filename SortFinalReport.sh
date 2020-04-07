pedfile=$1
codedir=$2
echo $pedfile
echo $codedir
head -n9 $pedfile  > Header.txt
Rscript ${codedir}/OrderPedFile.R $pedfile
cat Header.txt OrderedPedFile.txt > PEDFILE.txt
mv $pedfile ${pedfile}_unsorted
mv PEDFILE.txt $pedfile
