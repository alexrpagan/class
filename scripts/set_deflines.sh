#!/bin/bash

for f in `ls *.fa`
do
    sed -i '1s/.*/>'$f'/' $f && echo "Done with $f"
done
