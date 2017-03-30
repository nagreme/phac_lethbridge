#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# March 6, 2017
# Purpose: Parse a multifasta pangenome and extract only the accessory
#          (given by a list of accessory gene file)
#   sys.argv[1] = multifasta
#   sys.argv[2] = list of acc genes
#----------------------------------------------------------------------

import sys
import re
from Bio import SeqIO
#----------------------------------------------------------------------


def main():

    #open the gene list
    genes = open(sys.argv[2], 'rU')

    #create a list of the accessory genes
    acc_genes = []
    for gene in genes:
        acc_genes.append(gene.strip("\n\t '\""))

    #once done close that file
    genes.close()


    #Then move on to the multifasta file
    #Parse through it and extract only the accessory genes


    # outFilename = "{}_accessory.fasta".format(sys.argv[1][0:sys.argv[1].find(".")])

    # Alt. method
    #outFile = open(outFilename,'w')

    #with open(inFilename,'rU') as file:
    #    for record in SeqIO.parse(file,'fasta'):
    #        if re.search("plasmid", record.description) == None:
    #            print(record.id)


    input_seq_iterator = SeqIO.parse(inFilename,"fasta")
    not_plasmid_iterator = (record for record in input_seq_iterator if (re.search('plasmid',record.description) == None))
    SeqIO.write(not_plasmid_iterator, outFilename, "fasta")



if __name__ == '__main__':
    main()
