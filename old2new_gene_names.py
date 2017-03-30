#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# February 22, 2017
# Purpose: Parse gene_presence_absence.Rtab and change the old gene_id
#          gene names to the new gene names in a copy of the input file.
#          sys.argv[1] => oldgenename2seqid.txt
#          sys.argv[2] => seqid2newgenename.txt
#          sys.argv[3] => gene_presence_absence.Rtab
#----------------------------------------------------------------------

import sys
import re

#---------------------------------Main---------------------------------

def main():

    #open the input files
    oldnamesFile = open(sys.argv[1], 'rU')
    newnamesFile = open(sys.argv[2], 'rU')
    inFile = open(sys.argv[3], 'rU')

    #open output file
    outFile = open("{}/gene_pres_abs_new_genenames.Rtab".format(sys.argv[3][0:sys.argv[3].rfind("/")]), 'w')

    #build the old gene names dictionary
    old_names = {}
    for line in oldnamesFile:
        tokens = line.strip("\n").split("\t")
        old_names[tokens[1]] = tokens[0] #dnaA => 30628_00001

    #build the new gene names dictionary
    new_names = {}
    for line in newnamesFile:
        tokens = line.strip("\n").split("\t")
        new_names[tokens[0]] = tokens[1] #30628_00001 => dnaA

    #print the header as is
    outFile.write(inFile.readline())

    #for each line of the input file write it to the output file replacing the old gene name with the new one
    for line in inFile:
        old_gene_name = line[0:line.find("\t")].strip('"')
        new_gene_name = new_names[old_names[old_gene_name]]
        outFile.write("{new_name}{rest}".format(new_name=new_gene_name, rest=line[line.find("\t"):]))



    #close files
    oldnamesFile.close()
    newnamesFile.close()
    inFile.close()
    outFile.close()

    print("Done.\n")


if __name__ == '__main__':
    main()
