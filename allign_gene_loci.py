#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 18, 2017
# Purpose: Parse a tab-delimited file of gene loci (ID[start:end](strand)) (columns correspond to dif lists)
#          build lists for each column, sort them by key (ID and strand) then start position, and print
#          them to an output file aligning them while maintaining columns.

#Note: This program is functional but it places '' instead of blanks in the output
#      I have not yet found a way to supress this

#sys.argv[1] is the input filename (Required)
#sys.argv[2] is the fuzzy range to use (Optopnal)

import sys
import re
from operator import attrgetter



#------------------------------Constants-------------------------------
DELIM = "\t"
FUZZY_RANGE = None
FUZZY_FLAG = False

if (len(sys.argv) > 2):
    FUZZY_RANGE = int(sys.argv[2])
    FUZZY_FLAG = True


#----------------------------Class Gene--------------------------------
# A Class to hold a gene object => a genome accession and a locus
#----------------------------------------------------------------------
class Gene:
    def __init__(self, gene_locus_str):
        #build the object from a string in the input format
        #example: AL111168.1[1011761:1012148](+)
        #other: AL111168.1join{[1325667:1326117](+), [1326119:1326567](+), [1326565:1327107](+)}
        match = re.search('(^\w+[0-9]*\.[0-9])(\w?\w?\w?\w?\{?\[)([0-9]+)(:)([0-9]+)(\])(\([+-]\))(, )?(\[[0-9]+:[0-9]+\]\([-+]\), )*(\[)?([0-9]+)?(:)?([0-9]+)?(\]\([-+]\)\})?',gene_locus_str)

        self.id = match.group(1)
        self.strand = match.group(7)
        self.key = '{id}{strand}'.format(id = match.group(1), strand = match.group(7)) #AL111168.1(+)
        self.start = int(match.group(3)) #1011761

        if (match.group(13)):
            self.end = int(match.group(13))
        else:
            self.end = int(match.group(5)) #1012148

        self.len = self.end - self.start #387
        self.str_repr = gene_locus_str

    def __repr__(self):
        return self.str_repr

    #comparison method based on key then start locus
    #return -1 if self is before other
    #return 0 if self and othe have the same key and start
    #return 1 if self is after other
    #when fuzzy is True allow for a range to be equal
    def comparator(self, other):
        var = 0

        if (not FUZZY_FLAG):
            if (other == "" or self.key < other.key or (self.key == other.key and self.start < other.start)):
                var = -1
            elif (self.key > other.key or (self.key == other.key and self.start > other.start)):
                var = 1
        else: #fuzzy
            if (other == "" or self.key < other.key or (self.key == other.key and self.start < other.start - FUZZY_RANGE)):
                var = -1
            elif (self.key > other.key or (self.key == other.key and self.start > other.start + FUZZY_RANGE)):
                var = 1

        return var


#---------------------------------Main---------------------------------

def main():
    #open the input file
    inFile = open(sys.argv[1], 'rU')

    #build the output filename
    logic_type = "crisp"
    if (FUZZY_FLAG):
        logic_type = "fuzzy_{}".format(FUZZY_RANGE)

    outFilename = '{}_{}_aligned_lists.csv'.format(sys.argv[1][0:sys.argv[1].find(".")], logic_type)

    #Assume header in the first row so read the first line and do nothing wiht it
    inFile.readline()

    #get the number of columns from the number of tokens in a line
    tokens = inFile.readline().strip('\n').split(DELIM)
    num_cols = len(tokens)

    #build a list that will hold each column list
    columns = []

    #don't forget to add the current line to the right columns
    for i in range(num_cols):
        columns.append([])

        #only add a gene if the token has a valid length
        if (len(tokens[i]) > 0):
            columns[i].append(Gene(tokens[i]))

    #then read the rest of the lines and add the tokens to the right columns
    for line in inFile:
        tokens = line.strip('\n').split(DELIM)

        for i in range(num_cols):
            #only add a gene if the token has a valid length
            if (len(tokens[i]) > 0):
                columns[i].append(Gene(tokens[i]))

    #once done close the input file and sort the lists by key then start
    inFile.close()

    gene_lists = []

    for i in range(num_cols):
        gene_lists.append([])
        gene_lists[i] = sorted(columns[i], key = attrgetter('key','start'))


    # for i in range(num_cols):
    #     print("column",i,gene_lists[i],"\n")
    #
    # out_line = build_next_output_line(gene_lists)
    #
    # print(out_line)


    #open the output file
    outFile = open(outFilename, 'w')

    #keep building and printing output lines until gene_listsis empty
    while (not empty(gene_lists)):
        out_line = build_next_output_line(gene_lists)
        outFile.write(out_line)


    outFile.close()


#------------------------build_next_output_line------------------------
# PURPOSE: Build a line of output aligning the genes by selecting the
#          smallest in the current candidates at the top of the lists
#          and inserting delimiters as spacing. Remove the used items
#          from gene_lists.
# PARAMETERS: A list of columns containing Gene objects
# RETURNS: The next line of output tot be printed
#----------------------------------------------------------------------

def build_next_output_line(gene_lists):
    #finding the smallest linearly is probably fine (must consider ID in addition to start)

    temp_line = []
    for i in range(len(gene_lists)):
        if (not empty(gene_lists[i])):
            temp_line.append(gene_lists[i][0])
        else:
            temp_line.append("")

    #setup
    curr_smallest = temp_line[0]
    smallest_index = 0
    eq_found = False
    eq_overwrite = False

    #look through the rest of the list looking for the smallest element
    for i in range(1,len(temp_line)):
        #if we found a new smallest
        if (curr_smallest == "" or curr_smallest.comparator(temp_line[i]) == 1):
            #overwrite the previously smallest value with delimiter
            temp_line[smallest_index] = ""
            #then update our smallest info
            curr_smallest = temp_line[i]
            smallest_index = i

            #if we're overwriting a smallest that ahd several copies we need
            #to overwrite at least one more so set the overwrite flag
            if (eq_found):
                eq_found = False #and reset this for the new smallest
                eq_overwrite = True

        #if the curr item is equal toggle the equal flag
        elif (curr_smallest.comparator(temp_line[i]) == 0):
            eq_found = True

        #if the curr item is larger, then replace it with the delimiter
        else:
            temp_line[i] = ""

    #if at some point we did encounter an equality that was later replaced we need to overwrite those
    #at the same time, make sure to remove the used objects from gene_lists

    for i in range(len(temp_line)):
        #just go over the temp_line again and overwrite anything that doesn't
        #match the curr smallest if eq_overqrite was toggled
        if (temp_line[i] != ""):
            if (curr_smallest.comparator(temp_line[i]) != 0):
                if (eq_overwrite):
                    temp_line[i] = ""

            #pop off the used object so that it's not reused ad infinitum
            #note if it != 0 wasn't triggered it must be == since the previous loop should have caught < and >
            else:
                gene_lists[i].pop(0)

    #format output line and return it
    out_line = "{}\n".format("\t".join(repr(gene) for gene in temp_line))  #the output contains '' probably because of the repr function on this line

    print(out_line)

    return out_line

#-------------------------------empty----------------------------------
# PURPOSE: Recursively check that a multidimensional list is empty
# PARAMETERS: List to be checked
# RETURNS: True if list empty on all levels, otherwise False
# Obtained from stack overflow
#----------------------------------------------------------------------
def empty(seq):
    try:
        return all(empty(x) for x in seq)
    except TypeError:
        return False



if __name__ == '__main__':
    main()
