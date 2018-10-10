#!/bin/bash

TRA=($(for file in Tcyn_chr*_TWAR_*.maf; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

less ${tr}.maf | tr -d '\n' | awk -v var="${tr}" '{print ">" var "\n" $0}' > ~/laura_analyses/tfbs/fasta_TWAR_files/thylacine/${tr}.fa

done
