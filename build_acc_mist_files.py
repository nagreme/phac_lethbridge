#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# March 6, 2017
# Purpose: Given a multifasta containing all the sequences to include
#          build the files required for a MIST run (.fasta files for
#          each gene/allele and a markers file)
#   sys.argv[1] = multifasta
#   sys.argv[2] = accessory gene list (.txt)
#----------------------------------------------------------------------

import sys
import re
from Bio import SeqIO

#---------------------------------Main---------------------------------

def main():

    if (len(sys.argv) < 3):
        print("Please include a multifasta file and a accessory gene list as input.\n")

    else:
        #open the gene list
        genes = open(sys.argv[2], 'rU')

        #create a list of the accessory genes
        acc_genes = []
        for gene in genes:
            gene = gene.strip("\n\t\"")
            gene = re.sub("/","_",gene)
            gene = re.sub("\(","",gene)
            gene = re.sub("\)","",gene)
            gene = re.sub(" ","_",gene)
            gene = re.sub("'","",gene)
            gene = re.sub("\.","_",gene)
            gene = re.sub(",","_",gene)
            acc_genes.append(gene)

        #once done close that file
        genes.close()


        #ready the markers file
        markersFile = open("/home/nadege/mist/needed_files/acc.markers", 'w')

        #print the header
        markersFile.write("Marker Name\tTest Name\tTest Type\tForward Primer\tReverse Primer\tAmplicon Size (bp)\tAmplicon Range Factor (e.g. 0.1)\tAllelic Database Filename\tRepeat Size\n")

        for record in SeqIO.parse(sys.argv[1], "fasta"):
            #extract just the gene name
            match = re.search(".*? (.*)",record.description)
            curr_gene = match.group(1)

            #Take care of weirdness that could mess up the file name
            curr_gene = curr_gene.strip("\n\t\"")
            curr_gene = re.sub("/","_",curr_gene)
            curr_gene = re.sub("\(","",curr_gene)
            curr_gene = re.sub("\)","",curr_gene)
            curr_gene = re.sub(" ","_",curr_gene)
            curr_gene = re.sub("'","",curr_gene)
            curr_gene = re.sub("\.","_",curr_gene)
            curr_gene = re.sub(",","_",curr_gene)


            #if the gene is an accessory gene (i.e. in the list) then include it
            if (curr_gene in acc_genes):
                #write each gene to its own file
                outFile = open("/home/nadege/mist/needed_files/acc_alleles/{}.fasta".format(curr_gene), 'w')
                outFile.write(">1\n{}\n".format(record.seq))
                outFile.close()

                #and add an entry to the markers file for it
                markersFile.write("{gene}\tacc\t1\t\t\t-1\t0\t{gene}.fasta\t0\n".format(gene=curr_gene))





if __name__ == '__main__':
    main()
