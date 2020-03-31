pedfile=$1
head -n9 $pedfile  > Header.txt
Rscript OrderPedFile.R $pedfile
cat Header.txt OrderedPedFile.txt > PEDFILE.txt
mv $pedfile ${pedfile}_unsorted
mv PEDFILE $pedfile
