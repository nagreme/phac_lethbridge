#!/usr/bin/env python

import sys
import re


def main():
    shortList = open(sys.argv[1], 'r')
    longList = open(sys.argv[2],'r')
    outFile = open(differing_IDs.txt', 'w')

    for line in shortList:

        match = re.search('(NZ_)?(CP[0-9]+\.[0-9])',line)

        if match != None:
            outFile.write(match.group())
            outFile.write('\n')


    inFile.close()
    outFile.close()



if __name__ == '__main__':
    main()
