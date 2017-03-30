#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 13, 2017
# Purpose: Parse fasta files and make new copies excluding all
# sequences that include have an N._* identitfier

import sys
import re
from Bio import SeqIO

inFilename = sys.argv[1]


outFilename = inFilename[0:inFilename.find(".")] + "_no_duplicates.fasta"


input_seq_iterator = SeqIO.parse(inFilename,"fasta")
not_plasmid_iterator = (record for record in input_seq_iterator if (re.search('N._',record.description) == None))
SeqIO.write(not_plasmid_iterator, outFilename, "fasta")
