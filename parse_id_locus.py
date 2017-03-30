#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 23, 2017
# Purpose: Read a list of gene id/locus and parse out and format the
#          different parts nicely
#----------------------------------------------------------------------

import sys
import re


#------------------------------Constants-------------------------------
ID_LOCUS_INDEX = 0

#---------------------------------Main---------------------------------

def main():
    #open the input file
    inFile = open(sys.argv[1], 'rU')

    #open the output file
    outFilename = "{}_parsed.csv".format(sys.argv[1][0:sys.argv[1].find(".")])
    outFile = open(outFilename, 'w')

    #print header row in output file
    outFile.write("id\tstart\tend\tstrand\n")

    #assume the first line of the input file contains headers
    inFile.readline()

    #read in the input line by line and put into a list
    for line in inFile:

        match = re.search('(^\w+[0-9]*\.[0-9])(\w?\w?\w?\w?\{?\[)([0-9]+)(:)([0-9]+)(\])(\([+-]\))(, )?(\[[0-9]+:[0-9]+\]\([-+]\), )*(\[)?([0-9]+)?(:)?([0-9]+)?(\]\([-+]\)\})?',line)

        id_locus = match.group(1)
        start = match.group(3)
        end = match.group(5)
        strand = match.group(7)

        if (match.group(13)):
            end = match.group(13)

        #build the output line and write it to the output file
        output_line = "{}\t{}\t{}\t{}\n".format(id_locus, start, end, strand)

        outFile.write(output_line)


    #once we've gone through all the lines close your files
    inFile.close()
    outFile.close()




if __name__ == '__main__':
    main()
