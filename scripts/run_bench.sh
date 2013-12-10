#!/bin/bash

CABLAST='/home/apagan/code/cast_v0.9/cablast/cablast'
CLASS='/home/apagan/code/class/src/cpp/search/search'
NCBI='/home/apagan/misc/ncbi_cxx--12_0_0/GCC447-Release64'
BLAST=$NCBI'/bin/blastn'

QUERIES=/datadisk1/apagan/queries/queries_final.fa

IN='/datadisk1/apagan/outbench'
OUT='/datadisk1/apagan/benchresults'

mkdir -p $OUT

for size in 1 5 10 15 20 25
    do

    echo "Getting hits from blast cmd line"
    blast_dir=$OUT/blast_$size
    mkdir -p $blast_dir || exit
    time -p ( \
        $BLAST -db $IN/$size/fulldb/all \
               -outfmt 7                \
               -query $QUERIES          \
               -evalue 1e-30 > blastout 2> err > $blast_dir/out
    ) 2> $blast_dir/runtime.txt

    echo "Running class on $size chromosome copies"
    class_dir=$OUT/class_$size
    mkdir -p $class_dir || exit
    time -p (                             \
        $CLASS -db $IN/refdb/chr22        \
            -vdb $IN/$size/classdb/       \
            -chr 22                       \
            -ref $IN/ref/chr22.fa         \
            -in  $QUERIES                 \
            -full $IN/$size/fulldb/all    \
            -v -p 2> $class_dir/err > $class_dir/out \
    ) 2> $class_dir/runtime.txt

    echo "Running CaBLAST on $size chromosome copies"
    cablast_dir=$OUT/cablast_$size
    mkdir -p $cablast_dir || exit
    time -p (                                  \
        $CABLAST -fulldb  $IN/$size/fulldb/all \
                 -db      $IN/$size/cablastdb/uniques \
                 -links   $IN/$size/cablastdb/links.dat \
                 -uniques $IN/$size/cablastdb/uniques.fa \
                 -edits   $IN/$size/cablastdb/edits.dat \
                 -in $QUERIES 2> $cablast_dir/err > $cablast_dir/out \
    ) 2> $cablast_dir/runtime.txt

done