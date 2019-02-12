
#!bin/bash


## Laura E Cook, Unversity of Melbourne
## 1 Feb 2019

## This script takes a previously extracted sequence from a maf alignment (maf_extract_from_WGA.sh)
## Keeps the sequence only (column 7)
## Removes any new line
## And adds a fasta header for each TWAR to the file
## The ouput is the sequence in fasta format

TRA=($(for file in *.maf; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

cat ${tr}.maf | awk '{print $7; }' | tr -d '\n' | awk -v var="${tr}" '{print ">" var "\n" $0}' > ${tr}_gaps.fa

done
