#!/bin/bash

# ./extract_file_headers.sh > output.txt  (use in curr dir)

for f in *.fa
do
    echo $f
    grep ">" $f
done 
