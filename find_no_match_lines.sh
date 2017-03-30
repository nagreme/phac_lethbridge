#!/bin/bash

#Usage: ./find_no_match_lines.sh longList_filename shortList_filename > output.txt

while IFS= read -r line
do
    grep -q "$line" $2 || echo $line
        
done < $1
