#!/usr/bin/env python

#Give a txt file with a list of ncbi ids to retrieve as fasta files (.fa)
#(One id per line)

import sys
import os

def main():
    file = open(sys.argv[1], 'rU')

    for line in file:
        cmd = "echo " + line.strip()
        os.system(cmd)


        outfile = "%s.fasta" %(line.strip())

        cmd = "esearch -db assembly -query \"Campylobacter jejuni [ORGN] AND %s [SQID]\" | elink -target nuccore | efetch -format fasta > %s" %(line.strip(),outfile)
        #print cmd
        os.system(cmd)


    file.close()

if __name__ == '__main__':
    main()
