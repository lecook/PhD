
#!bin/bash


## Laura E Cook, Unversity of Melbourne
## 4 Feb 2019
## Remove alignment gaps from a fasta file, output with new name

TRA=($(for file in *.fa; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

sed 's/[-]//g' ${tr}.fa > ${tr}_nogaps.fa


done
