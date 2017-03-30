#!/usr/bin/env python

# Written by Nadège Pulgar-Vidal
# January 16, 2017
# Purpose: Create a dictionary for each file and compare the key sets

import sys
import re
from Bio import SeqIO


file1 = open(sys.argv[1], 'rU') #prokka
file2 = open(sys.argv[2], 'rU') #rast
file3 = open(sys.argv[3], 'rU') #ncbi

dict1 = {}
dict2 = {}
dict3 = {}

for line in file1:
    tokens = line.split('\t')
    key = tokens[0] + tokens[4]
    value = tokens[5]

    dict1[key] = value

for line in file2:
    tokens = line.split('\t')
    key = tokens[0] + tokens[4]
    value = tokens[5]

    dict2[key] = value

for line in file3:
    tokens = line.split('\t')
    key = tokens[0] + tokens[4]
    value = tokens[5]

    dict3[key] = value


set1 = set(dict1.keys())
set2 = set(dict2.keys())
set3 = set(dict3.keys())


# Write intersection list then difference lists
# outFile = open("output.txt",'w')
# outFile.write("Intersecting keys:\n")
# outFile.write(repr(set1 & set2))
# outFile.write("\n\nKeys unique to set 1:\n")
# outFile.write(repr(set1 - set2))
# outFile.write("\n\nKeys unique to set 2:\n")
# outFile.write(repr(set2 - set1))
# outFile.close()
# file1.close()
# file2.close()



outFile = open("intersection_prokka_rast_ncbi.txt",'w')
inter123 = (set1 & set2) & set3
for item in inter123:
    outFile.write(item)
    outFile.write("\n")
outFile.close()


outFile = open("intersection_prokka_rast.txt",'w')
inter12 = (set1 & set2) - inter123
for item in inter12:
    outFile.write(item)
    outFile.write("\n")
outFile.close()

outFile = open("intersection_prokka_ncbi.txt",'w')
inter13 = (set1 & set3) - inter123
for item in inter13:
    outFile.write(item)
    outFile.write("\n")
outFile.close()

outFile = open("intersection_rast_ncbi.txt",'w')
inter23 = (set2 & set3) - inter123
for item in inter23:
    outFile.write(item)
    outFile.write("\n")
outFile.close()


outFile = open("unique_prokka.txt",'w')
unique1 = (set1 - set2) - set3
for item in unique1:
    outFile.write(item)
    outFile.write("\n")
outFile.close()


outFile = open("unique_rast.txt",'w')
unique2 = (set2 - set1) - set3
for item in unique2:
    outFile.write(item)
    outFile.write("\n")
outFile.close()

outFile = open("unique_ncbi.txt",'w')
unique3 = (set3 - set1) - set2
for item in unique3:
    outFile.write(item)
    outFile.write("\n")
outFile.close()
