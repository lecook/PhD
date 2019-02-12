#!usr/bin/env python3

## Laura E Cook
## Nov 2018
## This script counts the number of gaps ("-") in a fasta alignment file

import sys
import re

# get filename from the command line
filename = sys.argv[1]

orig_stdout = sys.stdout
output = sys.argv[2]
f=open(output, 'w')
sys.stdout = f 


#open fasta file
with open(filename) as file:

	#create an empty variable
  datalist = []

  #for each line in a file
  for line in file:

  	# if the line start with a fasta header
    if line.startswith('>'):

    	# append the dataList with the fasta header and sequence as key and value
      datalist.append([line.strip()[1:], ''])
    else:
      datalist[-1][1] += line.strip()
     

  # then for the data in datalist print the header and count the number of times "-" occurs
  for data in datalist:
    matches1 = re.findall(r'[A-Za-z]-[A-Za-z]', data[1])
    matches2 = re.findall(r'[A-Za-z]--[A-Za-z]', data[1])
    matches3 = re.findall(r'[A-Za-z]---[A-Za-z]', data[1])
    matches4 = re.findall(r'[A-Za-z]----[A-Za-z]', data[1])
    matches5 = re.findall(r'[A-Za-z]-----[A-Za-z]', data[1])
    matches6 = re.findall(r'[A-Za-z]------[A-Za-z]', data[1])
    matches7 = re.findall(r'[A-Za-z]-------[A-Za-z]', data[1])
    matches8 = re.findall(r'[A-Za-z]--------[A-Za-z]', data[1])
    matches9 = re.findall(r'[A-Za-z]---------[A-Za-z]', data[1])
    matches10 = re.findall(r'[A-Za-z]----------[A-Za-z]', data[1])
    matches11 = re.findall(r'[A-Za-z]-----------[A-Za-z]', data[1])
    matches12 = re.findall(r'[A-Za-z]------------[A-Za-z]', data[1])
    matches13 = re.findall(r'[A-Za-z]-------------[A-Za-z]', data[1])
    matches14 = re.findall(r'[A-Za-z]--------------[A-Za-z]', data[1])
    matches15 = re.findall(r'[A-Za-z]---------------[A-Za-z]', data[1])
    matches16 = re.findall(r'[A-Za-z]----------------[A-Za-z]', data[1])
    matches17 = re.findall(r'[A-Za-z]-----------------[A-Za-z]', data[1])
    matches18 = re.findall(r'[A-Za-z]------------------[A-Za-z]', data[1])
    matches19 = re.findall(r'[A-Za-z]-------------------[A-Za-z]', data[1])
    matches20 = re.findall(r'[A-Za-z]--------------------[A-Za-z]', data[1])
    print(data[0], '    ', len(matches1), len(matches2), len(matches3), len(matches4), len(matches5), len(matches6), len(matches7), len(matches8), len(matches9), len(matches10), len(matches11), len(matches12), len(matches13), len(matches14), len(matches15), len(matches16), len(matches17), len(matches18), len(matches19), len(matches20))

sys.sdout = orig_stdout
f.close() 