#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 25, 2017
# Purpose: Take a tab-delimited input file containing start, end,
#          strand, and annotation and build a .gbk file minus the
#          header information.
#          Note: the files I'm using as input seem to have shifted the
#          start values to a 0 index rather than the original 1 index
#          so I will be doing +1 to correct that
#
#          Give the input file as the first cmd line argument and then
#          give all the extracted gbk info files (_genes_info.csv)
#          after that. (Added build_dicitonary Feb 2nd, 2017)
#----------------------------------------------------------------------

import sys
import re

#------------------------------Constants-------------------------------
#Column indices of input data
START_INDEX = 0
END_INDEX = 1
STRAND_INDEX = 2
ANNOTATION_INDEX = 3

#for build dictionary
PID = "Pid"
GENE_NAME = "GN"

#---------------------------------Main---------------------------------

def main():

    #Get that dictionary built
    gene_name_dict = build_dictionary()

    #open the input file
    inFile = open(sys.argv[1], 'rU')

    #Assume the first line of the input file contains a header row
    inFile.readline()

    #open the output file
    outFilename = "{}gene_gbk_body.txt".format(sys.argv[1][sys.argv[1].rfind("/")+1:sys.argv[1].rfind("merged")])
    outFile = open(outFilename, 'w')

    #set up a locus tag counter
    id_tag = 1 #biologists like to start numbering things at 1

    #for each line of input
    for line in inFile:
        #parse out the relevant pieces
        tokens = line.split("\t")

        start = int(tokens[START_INDEX]) + 1
        end = int(tokens[END_INDEX])
        strand = tokens[STRAND_INDEX]
        annotation = tokens[ANNOTATION_INDEX].strip("\n\"[]\"' ")

        #build a formatted CDS entry for that line
        locus_line = "{start}..{end}".format(start=start,end=end)

        if (strand == "(-)"):
            locus_line = "complement({})".format(locus_line)


        #Determine what to use for gene_name using the dictionary that was built
        gene_name = ""

        #check if there is an entry for this key in the dictionary
        gene_key = "{start}{end}{strand}".format(start=start-1,end=end,strand=strand)
        if (gene_key in gene_name_dict):
            gene_name = gene_name_dict[gene_key] #use the entry if present
        else:
            gene_name = annotation #default if there is no gene name in the dictionary

        #and tack that onto the end of locus line if there is a gene name
        if (gene_name):
            locus_line = "{}\n                     /gene=\"{}\"".format(locus_line,gene_name[gene_name.find(":")+1:])


        #the section's label will be CDS in most cases
        label = "     CDS             "

        #but if we're dealing with RNA it'll be different
        if (annotation.find("tRNA-",0) > -1 and len(annotation) < 14): #these length tests are meant to filter out synthases and the like
            label = "     tRNA            "

        elif (annotation.find("RNA",0) > -1 and len(annotation) < 20): #but they might unfortunately not catch some rRNA and tRNA => change those manually
            label = "     rRNA            "

        #If the label is CDS, we need these in the section
        codon_and_transl = "\n                     /codon_start=1\n                     /transl_table=11"

        #If the annotation is pseudogene (____), we need to include the /pseudo label
        if (annotation.find("pseudogene") > -1):
            codon_and_transl = "\n/pseudo{}".format(codon_and_transl)

        #But if the label is rRNA or tRNA then we can omit those
        if (label.find("CDS") == -1):
            codon_and_transl = ""


        #Finally build the section from the assembled pieces #Add /gene before the /locus_tag *************************************************#################################
        cds = '''{label}{locus}
                     /locus_tag="{filename}_{id_tag}"{codon_and_transl}
                     /product="{product}"\n'''.format(label=label,locus=locus_line, gene_name=gene_name, filename=sys.argv[1][sys.argv[1].rfind("/")+1:sys.argv[1].rfind("_merged")],
                     id_tag=id_tag,codon_and_transl=codon_and_transl, product=annotation)

        #update the counter
        id_tag += 1

        #print it to the output file
        outFile.write(cds)

    #close files
    inFile.close()
    outFile.close()



#---------------------------build_dictionary---------------------------
# PURPOSE: Fix it function to add gene_name annotation to already made
#          .gbk files. If I ever need to do this again I should modify
#          main to take a more complete input file and add those as I'm
#          building the rest (and therefore not use this fuinction at all)
# RETURNS: A dicitonary containing key-value pairs of keys being a conca-
#          tenation of start+end+strand and values being the gene names
#----------------------------------------------------------------------
def build_dictionary():
    #set up and empty dictionary
    gene_name_dict = {}

    #each file after the main input file is an input file for this function
    #for each one
    for index in range(2,len(sys.argv)):
        #open it
        dictFile = open(sys.argv[index], 'rU')

        #for each line
        for line in dictFile:
            #split it up for access
            tokens = line.split("\t")

            #only bother adding a key-value pair to the dictionary if there is a value (selected gene name)
            if (tokens[7]):
                #convert strand annotation for easier lookup later
                strand = "(+)"
                if (tokens[4] == "-1"):
                    strand = "(-)"

                #build a key
                gene_key = "{start}{end}{strand}".format(start=tokens[2], end=tokens[3], strand=strand)

                #insert into dicitonary with corresponding value
                insert_into_gene_dict(gene_name_dict, gene_key, tokens[7].strip('"\''))


        #close that file so we can move on to the next
        dictFile.close()

    #once you've gone through all the files and the dictionary is complete, return it
    return gene_name_dict



#----------------------------------------------------------------------
# PURPOSE: Given a gene name dictionary and a key value pair insert it
#          into the given dictionary. This means checking if it is
#          already there or not, and if there is a collision handle it.
# PARAMETERS: A gene name dicitonary, a key, and it's corresponding value
#----------------------------------------------------------------------
def insert_into_gene_dict(gene_name_dict, key, new_value):
    #first check if the key is already in the dictionary
    #if it isn't, simply add the entry
    if (key not in gene_name_dict):
        gene_name_dict[key] = new_value

    #if it is resolve the collision (GN:gene_name > Pid: protein_id)
    else:
        #first check if it really is a collision: are the values the same?
        old_value = gene_name_dict[key]
        if (old_value != new_value):
            #then handle the replacement cases
            #if curr is Pid and new is GN
            if (PID in old_value and GENE_NAME in new_value):
                #the new one replaces the old one
                gene_name_dict[key] = new_value

            #if both are the same type but have different values, keep both
            elif ((PID in old_value and PID in new_value) or (GENE_NAME in old_value and GENE_NAME in new_value)):
                #check if the keys are a substring of one another (keep the longest one)
                if (old_value in new_value or new_value in old_value):
                    if (len(old_value) < len(new_value)):
                        gene_name_dict[key] = new_value
                    #else old_value is the longest so no replacement

                else:
                    gene_name_dict[key] = "{}_or_{}".format(old_value,new_value[new_value.find(":")+1:])

            #if the old if GN and the new is Pid there is no replacement







if __name__ == '__main__':
    main()
    print("Done.\n")
