#!/bin/bash

count=0

while IFS= read -r var
do
    echo "$var"
    let count=count+1
done < $1

echo $count   
