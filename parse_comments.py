#!usr/bin/env python

import sys
import re


def main():
    inFile = open(sys.argv[1],'r')
    outFile = open(sys.argv[1]+'.parsed.csv', 'w')

    for line in inFile:
        if "#" not in line:
            outFile.write(line)


    inFile.close()
    outFile.close()



if __name__ == '__main__':
    main()
