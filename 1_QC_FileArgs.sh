#!/bin/bash

#first argument is the PLINK ped and map file name, second is the chip
#Thirs argument is the path to plink
#Fourth argument is the temporary file for genotypes
#Edit path to script! 
#Run within the directory of the files
plink=$3
export S=$4

echo "##########################################################################"
echo "Check individual call rate"
echo "##########################################################################"

echo "##########################################################################"
echo "Check for genotype call, heterozygosity"
echo "##########################################################################"
$3 --file $1 --missing --het --cow --recode --out $1_$2


echo "##########################################################################"
echo "Extrac the individuals that have genotype call less than 0.9"
echo "##########################################################################"
tail -n +2 $1_$2.imiss | awk '$6 >= 0.10 { print $1,$2 }' > IDsWithCallRateLessThan0.90$1_$2.txt


echo "##########################################################################"
echo "Draw plot for percentage of genotype call rate and heterozygosity rate for individuals"
echo "##########################################################################"
## Draw plot for percentage of genotype call rate and heterozygosity rate for individuals
## You need geneplotter package to run next script. Install geneplotter with the following command
# source("https://bioconductor.org/biocLite.R")
# biocLite("geneplotter")
${S}/DrawHetMissPlotInd.R $1_$2
echo "##########################################################################"
echo "Created plot:$1_$2_imiss-vs-het.pdf"
echo "##########################################################################"


## Remove individuals with low genotype call rate and heterozygosity rate out of range +-6SD of heterozygosity
FILE=./IndividualsToExlcudeMissHet.txt
if [ -f "$FILE" ]; then
      echo "##########################################################################"
      echo "File IndividualsToExclude.txt exists"
      echo "##########################################################################"
      $3 --file $1_$2 --remove IndividualsToExlcudeMissHet.txt --missing --recode --cow --out $1_$2_CleanInds
      mv IndividualsToExlcudeMissHet.txt IndividualsToExlcudeMissHet_$1_$2.txt
else 
      echo "##########################################################################"
      echo "File IndividualsToExclude.txt DOESNT exists"
      echo "##########################################################################"
      mv $1_$2.ped $1_$2_CleanInds.ped
      mv $1_$2.map $1_$2_CleanInds.map
      mv $1_$2.lmiss  $1_$2_CleanInds.lmiss
fi

FILE2=$1_$2_CleanInds.ped
if [ -f $FILE2 ]; then 
        echo "##########################################################################"
        echo "File CleanIndsMarkers.txt exists"
        echo "##########################################################################"
	echo "##########################################################################"
	echo "Check for SNP call rate"
	echo "##########################################################################"
	## Extract SNP with missing information on more than 10% of genotypes
	tail -n +2 $1_$2_CleanInds.lmiss | awk '$5 > 0.10 { print $1,$2 }' > MarkersWithCallRateLessThan0.90Plink_$1_$2.txt

	## Draw plot for percentage of genotype call rate and heterozygosity rate for markers
	${S}/DrawMissPlotMarkers.R $1_$2_CleanInds
	echo "Created plot:$1_$2_CleanInds_lmiss.pdf"

	## Remove SNP with missing information on more than 10% of genotypes
	$3 --file $1_$2_CleanInds --recode --geno 0.10 --cow --out $1_$2_CleanIndsMarkers
fi
