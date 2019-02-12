#!usr/bin/env python3

#Note, maf_parse assumes that each input maf alignment is against
#a single reference scaffold. In other words, you don't have any
#way to tell it to search against mouse chromosome 10 because
#the file should only contain alignments against mouse chromosome 10

#After running, use these sed commands to remove unwanted lines

#sed -i '/#eof/d' [myfiles].maf
#sed -i '/##maf version/d' [myfiles].maf
#sed -i '/maf_parse/d' [myfiles].maf

#And run this to add a header line to the file
#sed -i '1s/^/##maf version=1 scoring=roast.v3.3\n/' [myfiles].maf

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



