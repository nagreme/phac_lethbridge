#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 20, 2017
# Purpose: Read in a tab-delimited file with 5 columns: id/locus,
#          prokka_annotation, rast_annotation, ncbi_annotation without
#          hypotheticals, and ncbi with hypotheticals. For each
#          gene (id/locus) choose one annotation. Print out this list
#          of genes with annotations to a new file.
#          This program REQUIRES a specific format.
#----------------------------------------------------------------------

import sys
import re


#------------------------------Constants-------------------------------
ID_LOCUS_INDEX = 0
PROKKA_INDEX = 1
RAST_INDEX = 2
NCBI_NH_INDEX = 3
NCBI_H_INDEX = 4

#---------------------------------Main---------------------------------

def main():
    #open the input file
    inFile = open(sys.argv[1], 'rU')

    #open the output file
    outFilename = "{}_selected.csv".format(sys.argv[1][0:sys.argv[1].find(".")])
    outFile = open(outFilename, 'w')

    #print header row in output file
    outFile.write("id\tstart\tend\tstrand\tannotation\n")

    #assume the first line of the input file contains headers
    inFile.readline()

    #read in the input line by line and put into a list
    for line in inFile:
        row = line.split('\t')

        #strip off the [] and the '' for the products/annotations if there is one
        for i in range(1,len(row)):
            if (row[i]):
                row[i] = row[i].strip("[]' \n")

        #no need to store everything because we have enough information
        #in one line to make a selection, build out output line and print it
        chosen_annotation = choose_annotation(row[PROKKA_INDEX], row[RAST_INDEX], row[NCBI_NH_INDEX], row[NCBI_H_INDEX])


        #extract some usefel information for printing (for easier sorting and manipulation of the output)
        match = re.search("(\w{2}[0-9]*\.[0-9])(\[)([0-9]*)(:)([0-9]*)(\])(\([-+]\))",row[ID_LOCUS_INDEX])

        #build the output line and write it to the output file
        output_line = "{id}\t{start}\t{end}\t{strand}\t{annotation}\n".format(id=match.group(1),annotation=chosen_annotation,
                                                                                start=match.group(3), end=match.group(5), strand=match.group(7))

        outFile.write(output_line)


    #once we've gone through all the lines close your files

    inFile.close()
    outFile.close()




#---------------------------choose_annotation--------------------------
# PURPOSE: Wrapper for selecting the right annotation
# PARAMETERS: Prokka, RAST, and both NCBI annotations as strings
# RETURNS: The selected annotation
#----------------------------------------------------------------------
def choose_annotation(prokka_annotation, rast_annotation, ncbi_nh_annotation, ncbi_h_annotation):
    #get a case number
    case_num = eval_case_value_pa(prokka_annotation, rast_annotation, ncbi_nh_annotation)

    #evaluate the corresponding function to get an annotation
    func = annotation_function(case_num)
    annotation = func(prokka_annotation, rast_annotation, ncbi_h_annotation)

    return annotation


#--------------------------eval_case_value_pa--------------------------
# PURPOSE: Assign a case number based on which annotations are present
# PARAMETERS: The Prokka, RAST, and NCBI annotation strings
# RETURNS: An int representing the case number
#----------------------------------------------------------------------
def eval_case_value_pa(prokka_annotation, rast_annotation, ncbi_nh_annotation): #pa: presence-abscence
    case_num = 0

    if (prokka_annotation):
        case_num += 4

    if (rast_annotation):
        case_num += 2

    if (ncbi_nh_annotation):
        case_num += 1

    return case_num


#--------------------------eval_case_value_kw--------------------------
# PURPOSE: Assign a case number based on which annotations contain possible/putative
# PARAMETERS: The Prokka, RAST, and NCBI annotation strings
# RETURNS: An int representing the case number
#----------------------------------------------------------------------
def eval_case_value_kw(prokka_annotation, rast_annotation, ncbi_nh_annotation): #kw: keyword
    case_num = 0

    #Only give weight to an annotation if they do not contain the keywords

    if (prokka_annotation.find("possible") == -1 and prokka_annotation.find("putative") == -1):
        case_num += 4

    if (rast_annotation.find("possible") == -1 and rast_annotation.find("putative") == -1):
        case_num += 2

    if (ncbi_nh_annotation.find("possible") == -1 and ncbi_nh_annotation.find("putative") == -1):
        case_num += 1

    return case_num



#--------------------------annotation_function-------------------------
# PURPOSE: A function case/switch using dictionary mappings to selecting
#          which function should be used to select the correct annotation
#          based on which case is encountered.
# PARAMETERS: An int indicating which case has occured, corresponds to
#             which function should be used to selecct annotation
# RETURNS: The function that should be used to select annotation
#----------------------------------------------------------------------
def annotation_function(case_num):
    switcher = {
        0: zero_one,
        1: zero_one,
        2: two,
        3: three_five,
        4: four,
        5: three_five,
        6: six,
        7: seven
    }
    # Get the function from switcher dictionary
    func = switcher.get(case_num, lambda: "")

    return func


#----------------------------------------------------------------------
# Switch functions to select annotation
#----------------------------------------------------------------------

#-------------------------------zero_one-------------------------------
# 0 0 0
# 0 0 1
# NCBI wins by default
#----------------------------------------------------------------------
def zero_one(prokka_annotation, rast_annotation, ncbi_h_annotation):
    return ncbi_h_annotation


#---------------------------------two----------------------------------
# 0 1 0
# RAST wins by default
#----------------------------------------------------------------------
def two(prokka_annotation, rast_annotation, ncbi_h_annotation):
    return rast_annotation


#------------------------------three_five------------------------------
# 0 1 1
# 1 0 1
# If NCBI has a length greater than 5 it wins by default
# Otherwise the other one does (Prokka or RAST)
#----------------------------------------------------------------------
def three_five(prokka_annotation, rast_annotation, ncbi_h_annotation):
    annotation = ncbi_h_annotation

    if (len(ncbi_h_annotation) < 6):
        case_num = eval_case_value_pa(prokka_annotation, rast_annotation, ncbi_h_annotation) - 1
        #go to the prokka or rast case if ncbi is no good by subtracting one
        func = annotation_function(case_num)
        annotation = func(prokka_annotation, rast_annotation, "")

    return annotation


#---------------------------------four---------------------------------
# 1 0 0
# Prokka wins by default
#----------------------------------------------------------------------
def four(prokka_annotation, rast_annotation, ncbi_h_annotation):
    return prokka_annotation


#---------------------------------six----------------------------------
# 1 1 0
# Unless they're disfavoured by containing "possible/putative"
# keywords, Prokka wins by default
#----------------------------------------------------------------------
def six(prokka_annotation, rast_annotation, ncbi_h_annotation):
    annotation = prokka_annotation

    #Prokka wins unless it contains a keyword and RAST doesn't
    if ((prokka_annotation.find("possible") > -1 or prokka_annotation.find("putative") > -1)
        and rast_annotation.find("possible") == -1 and rast_annotation.find("putative") == -1):
        annotation = rast_annotation

    return annotation


#---------------------------------seven--------------------------------
# 1 1 1
# If NCBI has a length greater than 5: keyword search, NCBI > Prokka > RAST
# Else return six
#----------------------------------------------------------------------
def seven(prokka_annotation, rast_annotation, ncbi_h_annotation):
    annotation = ncbi_h_annotation

    #If NCBI is short choose one of Prokka and RAST
    if (len(ncbi_h_annotation) < 6):
        annotation = six(prokka_annotation, rast_annotation, "")

    #Otherwise, keyword search between all three to determine a new case number
    else:
        case_num = eval_case_value_kw(prokka_annotation, rast_annotation, ncbi_h_annotation)

        #if one of them was disfavoured the case_number will be smaller so execute that
        if (case_num < 7):
            func = annotation_function(case_num)
            annotation = func(prokka_annotation, rast_annotation, ncbi_h_annotation)

    return annotation



if __name__ == '__main__':
    main()
