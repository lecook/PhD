#!/bin/bash

## Laura E Cook, Unversity of Melbourne
## 1 Feb 2019

## This script extracts the aligned sequence for each TWAR in every species in the alignment
## Ideally it would be great to learn how to do this is a loop
## I could generate a txt file with the species genome names and then loop over those names?
## Look at make this more efficient and soft code when I come back to it


TRA=($(for file in chr*.maf; do echo $file |cut -d "." -f 1;done))

echo ${TRA[@]}

for tr in ${TRA[@]};

do

echo ${tr}
grep felCat5 ${tr}.maf  > species_mafs/felCat5_${tr}.maf
grep mm10 ${tr}.maf  > species_mafs/mm10_${tr}.maf
grep cavPor3 ${tr}.maf  > species_mafs/cavPor3_${tr}.maf
grep dipOrd1 ${tr}.maf  > species_mafs/dipOrd1_${tr}.maf
grep hetGla2 ${tr}.maf  > species_mafs/hetGla2_${tr}.maf
grep ochPri2 ${tr}.maf  > species_mafs/ochPri2_${tr}.maf
grep oryCun2 ${tr}.maf  > species_mafs/oryCun2_${tr}.maf
grep rn5 ${tr}.maf  > species_mafs/rn5_${tr}.maf
grep speTri2 ${tr}.maf  > species_mafs/speTri2_${tr}.maf
grep tupBel1 ${tr}.maf  > species_mafs/tupBel1_${tr}.maf
grep calJac3 ${tr}.maf  > species_mafs/calJac3_${tr}.maf
grep gorGor3 ${tr}.maf  > species_mafs/gorGor3_${tr}.maf
grep hg19 ${tr}.maf  > species_mafs/hg19_${tr}.maf
grep micMur1 ${tr}.maf  > species_mafs/micMur1_${tr}.maf
grep nomLeu2 ${tr}.maf  > species_mafs/nomLeu2_${tr}.maf
grep otoGar3 ${tr}.maf  > species_mafs/otoGar3_${tr}.maf
grep panTro4 ${tr}.maf  > species_mafs/panTro4_${tr}.maf
grep papHam1 ${tr}.maf  > species_mafs/papHam1_${tr}.maf
grep ponAbe2 ${tr}.maf  > species_mafs/ponAbe2_${tr}.maf
grep rheMac3 ${tr}.maf  > species_mafs/rheMac3_${tr}.maf
grep saiBol1 ${tr}.maf  > species_mafs/saiBol1_${tr}.maf
grep tarSyr1 ${tr}.maf  > species_mafs/tarSyr1_${tr}.maf
grep ailMel1 ${tr}.maf  > species_mafs/ailMel1_${tr}.maf
grep bosTau7 ${tr}.maf  > species_mafs/bosTau7_${tr}.maf
grep Clup ${tr}.maf  > species_mafs/clup_${tr}.maf
grep choHof1 ${tr}.maf  > species_mafs/choHof1_${tr}.maf
grep dasNov3 ${tr}.maf  > species_mafs/dasNov3_${tr}.maf
grep Tcyn ${tr}.maf  > species_mafs/tcyn_${tr}.maf
grep sarHar1 ${tr}.maf  > species_mafs/sarHar1_${tr}.maf
grep echTel1 ${tr}.maf  > species_mafs/echTel1_${tr}.maf
grep equCab2 ${tr}.maf  > species_mafs/equCab2_${tr}.maf
grep eriEur1 ${tr}.maf  > species_mafs/eriEur1_${tr}.maf
grep loxAfr3 ${tr}.maf  > species_mafs/loxAfr3_${tr}.maf
grep myoLuc2 ${tr}.maf  > species_mafs/myoLuc2_${tr}.maf
grep oviAri1 ${tr}.maf  > species_mafs/oviAri1_${tr}.maf
grep proCap1 ${tr}.maf  > species_mafs/proCap1_${tr}.maf
grep pteVam1 ${tr}.maf  > species_mafs/pteVam1_${tr}.maf
grep sorAra1 ${tr}.maf  > species_mafs/sorAra1_${tr}.maf
grep susScr3 ${tr}.maf  > species_mafs/susScr3_${tr}.maf
grep triMan1 ${tr}.maf  > species_mafs/triMan1_${tr}.maf
grep turTru2 ${tr}.maf  > species_mafs/turTru2_${tr}.maf
grep vicPac1 ${tr}.maf  > species_mafs/vicPac1_${tr}.maf
grep anoCar2 ${tr}.maf  > species_mafs/anoCar2_${tr}.maf
grep chrPic1 ${tr}.maf  > species_mafs/chrPic1_${tr}.maf
grep danRer7 ${tr}.maf  > species_mafs/danRer7_${tr}.maf
grep fr3 ${tr}.maf  > species_mafs/fr3_${tr}.maf
grep gadMor1 ${tr}.maf  > species_mafs/gadMor1_${tr}.maf
grep galGal4 ${tr}.maf  > species_mafs/galGal4_${tr}.maf
grep gasAcu1 ${tr}.maf  > species_mafs/gasAcu1_${tr}.maf
grep latCha1 ${tr}.maf  > species_mafs/latCha1_${tr}.maf
grep macEug2 ${tr}.maf  > species_mafs/macEug2_${tr}.maf
grep melGal1 ${tr}.maf  > species_mafs/melGal1_${tr}.maf
grep melUnd1 ${tr}.maf  > species_mafs/melUnd1_${tr}.maf
grep monDom5 ${tr}.maf  > species_mafs/monDom5_${tr}.maf
grep oreNil2 ${tr}.maf  > species_mafs/oreNil2_${tr}.maf
grep ornAna1 ${tr}.maf  > species_mafs/ornAna1_${tr}.maf
grep oryLat2 ${tr}.maf  > species_mafs/oryLat2_${tr}.maf
grep petMar1 ${tr}.maf  > species_mafs/petMar1_${tr}.maf
grep taeGut1 ${tr}.maf  > species_mafs/taeGut1_${tr}.maf
grep tetNig2 ${tr}.maf  > species_mafs/tetNig2_${tr}.maf
grep xenTro3 ${tr}.maf > species_mafs/xenTro3_${tr}.maf

done
