
#!bin/bash


## Laura E Cook, Unversity of Melbourne
## 1 Feb 2019

## This script takes the orthologous dna sequences for the TWARs and peforms a maf alignment
## This will be used to generate phylogenetic trees for comparison

TRA=($(for file in *_nogaps.fa; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

mafft ${tr}.fa > ${tr}_plusTammar_aligned.fa

done
