#!/bin/bash

TRA=($(for file in chr*.maf; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}

sed -i '/#eof/d' ${tr}.maf | sed -i '/##maf version/d' ${tr}.maf | sed -i '/maf_parse/d' ${tr}.maf | sed -i '1s/^/##maf version=1 scoring=roast.v3.3\n/' ${tr}.maf

done
