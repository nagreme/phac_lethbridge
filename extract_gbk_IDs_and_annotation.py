#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 16, 2017
# Purpose: Parse genbank files and extract gene IDs and
#          annotation (product name) and print these out

import sys
import re
from Bio import SeqIO

inFilename = sys.argv[1]


outFilename = inFilename[0:inFilename.find(".")] + "_annotated_gene_list.csv"
outFile = open(outFilename, 'w')

#print header for file
outFile.write("id\tname\tdesc\ttype\tlocation\tproduct\n")


for record in SeqIO.parse(inFilename, "genbank"):
    for feature in record.features:
        if ('product' in feature.qualifiers):
            out_str = '{id}\t{name}\t{desc}\t{type}\t{location}\t{product}\n'.format(id = record.id, name = record.name, desc = record.description, type = feature.type, location = feature.location, product = feature.qualifiers['product'])
            outFile.write(out_str)

outFile.close()
