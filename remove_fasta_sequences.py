#!/usr/bin/env python3

## Laura E Cook
## 9 Jan 2019

## Takes a multisequence fasta file and removes sequences based on a list of fasta headers


from Bio import SeqIO
import sys

ffile = SeqIO.parse(sys.argv[1], "fasta")
header_set = set(line.strip() for line in open(sys.argv[2]))

orig_stdout = sys.stdout
f=open('ailMel1_TWARs_keepGaps_80percent_overlap.fasta', 'w')
sys.stdout = f 

for seq_record in ffile:
    try:
        header_set.remove(seq_record.name)
    except KeyError:
        print(seq_record.format("fasta"))
        continue
if len(header_set) != 0:
    print(len(header_set),'of the headers from list were not identified in the input fasta file.', file=sys.stderr)

sys.sdout = orig_stdout
f.close()
