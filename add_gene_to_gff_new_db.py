#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 25, 2017
# Purpose: Parse .gff files made using Prokka with my new Campylobacter
#          DB to add gene=______ tags so Roary registers it and uses it
#          to annotate. The product name will be used as the gene name.
#----------------------------------------------------------------------

import sys
import re

#------------------------------Constants-------------------------------



#---------------------------------Main---------------------------------

def main():
    #open the input file
    inFile = open(sys.argv[1], 'rU')

    #open the output file
    outFilename = "{}_mod_gene.gff".format(sys.argv[1][0:sys.argv[1].find(".")])
    outFile = open(outFilename, 'w')

    #we haven't reached the ##FASTA section yet
    fasta_flag = False

    #read line by line (but break when you reach ##FASTA)
    for line in inFile:
        #just copy over the line as is if it's a comment or something else
        if ("#" in line or "repeat_region" in line or fasta_flag):
            outFile.write(line)
            if ("##FASTA" in line):
                fasta_flag = True

        else: #valid line
            #if it's not a comment parse it with a regex to extract the piece before
            match = re.search("(.*ID=.*?;)(.*;product=)(.*)", line)

            # print(line)
            # print(match.group(1),"\n", match.group(2), "\n", match.group(3), "\n\n")

            #gene=_____ should go, the part between, the product, and whatever may come after that
            #use the product to make the gene=_____ and reconstruct the entry using the parsed pieces
            new_line = "{first}gene={product};{mid}{product}\n".format(first=match.group(1), product=match.group(3), mid=match.group(2))

            #then print the reconstructed line to the output file
            outFile.write(new_line)

    #close files
    inFile.close()
    outFile.close()


if __name__ == '__main__':
    main()
