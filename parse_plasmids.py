#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 11, 2017
# Purpose: Parse fasta files and make new copies excluding all
# sequences that include the word plasmid in their description

import sys
import re
from Bio import SeqIO

inFilename = sys.argv[1]



outFilename = inFilename[0:inFilename.find(".")] + "_no_plasmids.fasta"

# Alt. method
#outFile = open(outFilename,'w')

#with open(inFilename,'rU') as file:
#    for record in SeqIO.parse(file,'fasta'):
#        if re.search("plasmid", record.description) == None:
#            print(record.id)


input_seq_iterator = SeqIO.parse(inFilename,"fasta")
not_plasmid_iterator = (record for record in input_seq_iterator if (re.search('plasmid',record.description) == None))
SeqIO.write(not_plasmid_iterator, outFilename, "fasta")
