#!bin/bash

## Laura E Cook, University of Melbourne
## 5 Oct 2018
## Take a mutiple sequence fasta file and outputs the length of each of the sequences in that fasta

cat file.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'

done