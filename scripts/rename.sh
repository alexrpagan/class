#!/bin/bash

for f in `ls *.fasta`
    do
        genome=`echo $f | sed 's/_22:/-/' | sed 's/\.fasta//'`
        mv $f $genome || exit 1
done
