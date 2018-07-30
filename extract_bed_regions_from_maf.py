#!usr/bin/env python3

import sys
import os


# Collect command line arguments. This script takes a bed file
# of coordinates in the maf reference genome and a maf whole
# genome alignment as inputs

bed_regions = sys.argv[1]
mafWGA = sys.argv[2]

# Open the bed file
with open(bed_regions, 'r') as inputfile:

# Iterate through bed file line-by-line
# without loading the whole file into the program.
	for line in inputfile:

  # split the bed region at whitespace (i.e. tabs)
        line = line.strip("\n")
        region = line.split("\t")
        print(region)
        # collect start position
        start_coordinate = region[1]

        # collect end position
        end_coordinate = region[2]

        # collect TWAR name
        TWAR = region[3]
        # format output file name
        outputFile = TWAR + ".maf"

        # format the maf_parse command
        bashCommand = "maf_parse --start {} --end {} {} > {}".format(start_coordinate, end_coordinate, mafWGA, outputFile)

        # run the command
        os.system(bashCommand)
