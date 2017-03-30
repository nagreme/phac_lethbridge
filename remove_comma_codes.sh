#!/bin/bash

#replace all %2C with a nothing in the given file

sed -e 's/%2C//g' $1 >$1.out

echo "Done."
