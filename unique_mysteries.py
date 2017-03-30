#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 25, 2017
# Purpose: For every product=mystery found add a number to make it
#          unique. Eventually change it to  process hypothetical
#          protein to mystery_#
#      >>> Feb 15, 2017 (continued Feb 21, 2017)
#          Modifications to make this script a metagenome prep/cleanup
#          parser. Change "hypothetical protein" --> mystery_#, add
#          Name=_; and gene=_; fields for those, and change all the
#          sequence IDs to match their ID=_; field to avoid confusion
#----------------------------------------------------------------------

import sys
import re

#---------------------------------Main---------------------------------

def main():
    #start at 1 for biologists
    myst_count = 1
    seq_count = 1
    rna_count = 1
    seq_id_dict = {}

    #open the input file
    inFile = open(sys.argv[1], 'rU')

    #open the output file
    outFilename = "{}_mystified.gff".format(sys.argv[1][0:sys.argv[1].find(".")])
    outFile = open(outFilename, 'w')
    refFile = open("{}/seqid2newgenename.txt".format(sys.argv[1][0:sys.argv[1].rfind("/")]), 'w')

    fasta_flag = False

    #read line by line
    for line in inFile:
        #assume we'll just copy over the line as is
        new_line = line

        #if it's a comment check for sequence-regions and ##FASTA
        if ("#" in line or "repeat_region" in line):
            #if it's a sequence make a hash so we can replace all other occurences
            if ("##sequence-region" in line):
                match = re.search("(##sequence-region )(.*_[0-9]*)( .*)",line)
                seq_id_dict[match.group(2)] = "0_metagenome_{:0>5}".format(seq_count) #left-pad with 0 up to 5 digits
                seq_count += 1
                new_line = "{first}{id}{last}\n".format(first=match.group(1), id=seq_id_dict[match.group(2)], last=match.group(3))

            #if we reach the fasta part toggle the flag
            elif ("##FASTA" in line):
                fasta_flag = True


        #if it's not a comment but we've passed the fasta flag
        elif (fasta_flag):
            #we're either dealing with a header, which we have to replace using our dicitonary
            if (">" in line):
                match = re.search(">(.*)",line)
                new_line = ">{}\n".format(seq_id_dict[match.group(1)])
            #or a line of nucleotides which we just copy over

        #Otherwise, it'll be a line with an entry
        else:
            #look for the sequence name and get the right replacement from our dictionary
            match = re.search("(.*?)\t.*", line)
            new_entry_title = seq_id_dict[match.group(1)]

            #dissect the line, split the parts
            match = re.search("(.*?)(\t.*ID=.*?;)(eC_num.*?;)?(Name=.*?;)?(gene=.*?;)?(inference=.*product=)(.*)", line)

            id_part = match.group(2)
            name_part = match.group(4)
            gene_part = match.group(5)
            middle = match.group(6)
            product = match.group(7)

            #if we're dealing with a hypothetical protein line we need to add Name and gene fields and modify product
            if ("product=hypothetical protein" in line):
                name_part = "Name=mystery_{num};".format(num=myst_count)
                gene_part = "gene=myst_{num};".format(num=myst_count)
                product = "mystery_{}".format(myst_count)
                myst_count += 1

            #for some reason some proteins fell through the cracks...
            if (gene_part == None):
                name_part = "Name={};".format(product)
                gene_part = "gene={};".format(product)

            #build the new line
            new_line = "{new_title}{id}{name_part}{gene_part}{middle}{product}\n".format(new_title=new_entry_title,
            id=id_part, name_part=name_part, gene_part=gene_part, middle=middle, product=product )

            #if we're dealing with rRNA it won't have a gene_part, but we need
            #one for gene_presence_absence.Rtab so make one up
            if ("rRNA\t" in line):
                gene_part = "rRNA_{:0>3};".format(rna_count)
                rna_count += 1
            #print a line of output to the seqid2newgenename file
            refFile.write("{old_seq_id}\t{new_gene_name}\n".format(old_seq_id=match.group(1), new_gene_name=gene_part[5:len(gene_part)-1]))


        #then print the new line to the output file (whatever it currently is)
        outFile.write(new_line)

    #end for loop

    #close files
    inFile.close()
    outFile.close()
    refFile.close()


if __name__ == '__main__':
    main()
