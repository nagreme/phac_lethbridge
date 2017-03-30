#!/usr/bin/env python

import sys
import re


def main():
    inFile = open(sys.argv[1],'r')
    outFile = open(sys.argv[1]+'.extracted_IDs.txt', 'w')

    for line in inFile:

        match = re.search('(NZ_)?(CP[0-9]+\.[0-9])',line)

        if match != None:
            outFile.write(match.group())
            outFile.write('\n')


    inFile.close()
    outFile.close()



if __name__ == '__main__':
    main()
