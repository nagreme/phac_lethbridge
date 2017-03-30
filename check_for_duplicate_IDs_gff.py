#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# February 3, 2017
# Purpose: Parse Prokka .gff files chekcing IDs and pointing out duplicates
#----------------------------------------------------------------------

import sys
import re

#------------------------------Constants-------------------------------



#---------------------------------Main---------------------------------

def main():

    id_list = []

    #collect all the IDs from all the files and then check if there are duplicates
    for index in range(1,len(sys.argv)):
        #open the input file
        inFile = open(sys.argv[index], 'rU')

        #we haven't reached the ##FASTA section yet
        fasta_flag = False

        #read line by line (but break when you reach ##FASTA)
        for line in inFile:
            #just copy over the line as is if it's a comment or something else
            if ("#" in line or "repeat_region" in line or fasta_flag):
                if ("##FASTA" in line):
                    fasta_flag = True

            else: #valid line
                #if it's not a comment parse it with a regex to extract the piece before
                match = re.search(".*ID=(.*?);.*", line)
                #add the ID part to our list of IDs
                id_list.append(match.group(1))
                print(match.group(1))

        #close files
        inFile.close()

    #check if there are duplicates first
    id_set = set(id_list)

    if (len(id_set) != len(id_list)):
        print("There are duplicate IDs")


if __name__ == '__main__':
    main()
