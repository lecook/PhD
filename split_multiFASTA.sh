#!bin/bash

## Laura E Cook, University of Melbourne
## 6 Feb 2019
## Split a mutiFASTA file into single fasta files

while read line
do
    if [[ ${line:0:1} == '>' ]]
    then
        outfile=${line#>}.fa
        echo $line > $outfile
    else
        echo $line >> $outfile
    fi
done < macEug3_missing_twars.fa