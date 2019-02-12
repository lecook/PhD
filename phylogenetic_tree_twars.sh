
#!bin/bash


## Laura E Cook, Unversity of Melbourne
## 4 Feb 2019
## This script generates a phylogenetic tree for a given MAFFT alignment
## Input is an alignment as a FASTA file and output is a newick tree that can be plotted in R

TRA=($(for file in *_plusTammar_aligned.fa; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

mafft --retree 0 --treeout --reorder --localpair --weighti 0 --averagelinkage ${tr}.fa


done
