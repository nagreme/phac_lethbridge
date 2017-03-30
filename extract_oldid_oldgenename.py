#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# February 22, 2017
# Purpose: Parse the pan_genome_reference.fa (sys.argv[1]) file and
#          print gene_id\tgene_name\n to a file for later use.
#----------------------------------------------------------------------

import sys
import re
from Bio import SeqIO

#---------------------------------Main---------------------------------

def main():

    #open output file
    refFile = open("{}/oldgenename2seqid.txt".format(sys.argv[1][0:sys.argv[1].rfind("/")]), 'w')

    #for each entry extract the gene id and description/annotation/gene name
    for record in SeqIO.parse(sys.argv[1], "fasta"):
        desc = record.description
        #get the individual pieces
        match = re.search("(.*?) (.*)",desc)

        #write them to the output file in the wanted format
        refFile.write("{}\t{}\n".format(match.group(1),match.group(2)))

    #close the file
    refFile.close()


if __name__ == '__main__':
    main()
