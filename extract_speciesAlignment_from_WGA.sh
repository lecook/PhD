#!/bin/bash

TRA=($(for file in chr*.maf; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

grep Tcyn ${tr}.maf | sed 's/ \+/\t/g'| cut -f 7 > Tcyn_${tr}.maf

done
