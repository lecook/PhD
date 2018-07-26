#!usr/bin/env python3

# This is the final script I used to split the BAM files into smaller chunks for analyses
# This python script allows you to take a BAM file and split it into X chunks (chunks) where each chunk is a separate SAM file
# containing all alignments against Y scaffolds (scaffoldsPerChunk).
# This sam file can then be converted into a BAM file.

# import python modules used in the script
import sys
import os
import math
import random

# get arguments from the command line (list of scaffold names, number of chunks, and prefix for the output files)
# scaffold names is a text file with each name of the scaffold on a separate line
scaffoldNames = sys.argv[1]
chunks = int(sys.argv[2])
prefix = str(sys.argv[3])

# hard-coded files - change this later
mergedBAM = "/home/laura/laura_analyses/refassembly/Thylacine_trimmed_qual20_minlen50_DCed_OCed.marked_duplicates.bam"


### ---------------------- STORE SCAFFOLD NAMES IN A LIST ---------------------- ###

# create empty list variable to store each scaffold name
scaffoldList = []

# open the scaffold.txt file and save each scaffold name into a list
# in a for loop read each line (scaffold names are all on an individual line) and store each name in a list.
# append the list with the next scaffold name

# open txt file which contains the scaffolds names each on an individual line (from command arguments)
with open(scaffoldNames, 'r') as inputfile:

    #in a for loop read each line and store each name in a list
    for line in inputfile:

        # append the list with the next scaffold name
        scaffoldList.append(line.strip())


# randomised the list because the genome is sorted by largest scaffold first so the first BAM file will be huge compared to the next etc.
random.shuffle(scaffoldList)

# count the number of scaffolds, divide by the number of chunks you want (from command arguments) and round to the nearest integer
scaffoldsPerChunk = math.ceil(len(scaffoldList)/chunks)


### ------------- EXTRACTING ALIGNMENTS INTO SEPARATE CHUNKS ------------- ###

# this will be used to keep track of how many scaffolds have been put into each file
scaffoldCounter = 0

# this will be used to keep track of which output chunk file the scaffold is being extracted to.
chunkCounter = 1

# make the output file name for chunk1 once  outside the loop just to output the sam header
outputFile = prefix + str(chunkCounter) + '.sam'

## make a bash command to run samtools to extract the sam header
bashCommand = "samtools-1.8 view -H /home/laura/laura_analyses/refassembly/Thylacine_trimmed_qual20_minlen50_DCed_OCed.marked_duplicates.bam > {}".format(outputFile)

## run the command to add the sam header
os.system(bashCommand)

# for every scaffold in the list
for scaffold in scaffoldList:

    # create an output file name for the scaffolds.
    # concatenate the prefix (from the command line arguments), with the chunk number for the current output file and the .bam suffix
    # chunkCounter must be a string to concatenate use str()
    outputFile = prefix + str(chunkCounter) + '.sam'
    
    # create a string for the samtools command
    # using .format() at the end will replace the  first {} with the the large bam file, the second {} with the scaffold names 
    # and the third {} with the outputFile name specified above
    bashCommand = "samtools-1.8 view {} {} >> {}".format(mergedBAM, scaffold, outputFile)

    # run the string above as a command line command.
    os.system(bashCommand)

    # increment your scaffold counter by 1 to indicate 1 scaffold has been added to this output file
    scaffoldCounter += 1;

    # check to see if the number of scaffolds in the current output chunk file has passed the threshold set by (scaffoldsPerChunk)
    if (scaffoldCounter >= scaffoldsPerChunk):



        # increment the chunkCounter to move onto a new chunk file
        int(chunkCounter)
        chunkCounter += 1

        # set the scaffoldCounter back to 0 so that the for loop will start counting up to scaffoldTotal again
        scaffoldCounter = 0

        # make the output file name for chunk1 once  outside the loop just to output the sam header
        outputFile = prefix + str(chunkCounter) + '.sam'

        ## make a bash command to run samtools to extract the sam header
        bashCommand = "samtools-1.8 view -H /home/laura/laura_analyses/refassembly/Thylacine_trimmed_qual20_minlen50_DCed_OCed.marked_duplicates.bam > {}".format(outputFile)

        ## run the command to add the sam header
        os.system(bashCommand)
