#!/bin/bash

CABLAST='/home/apagan/code/cast_v0.9/cablast'
CLASS='/home/apagan/code/class/src/cpp/search/search'
NCBI='/home/apagan/misc/ncbi_cxx--12_0_0/GCC447-Release64'
BLAST=$NCBI/bin/blastn
MAKEBLASTDB=$NCBI/bin/makeblastdb
SEED=1
COMBINED=all.fa
CHR=chr22

# holds a database for each method
OUT='/datadisk1/apagan/outbench'

# raw fasta, bitvectors, reference genome
IN='/datadisk1/apagan/benchin'

if ! [[ -d $IN ]];
    then
        echo "Input dir '$IN' not found"
        exit 1
fi

mkdir -p $OUT

if ! [[ -e $OUT/ref ]]
    then
        echo "Copying reference"
        mkdir $OUT/ref            || exit 1
        cp $IN/ref/$CHR* $OUT/ref || exit 1
fi

if ! [[ -e $OUT/refdb ]]
    then
        echo "Making reference DB"
        mkdir $OUT/refdb || exit
        $MAKEBLASTDB -in $OUT/ref/$CHR.fa \
            -out $OUT/refdb/$CHR          \
            -dbtype nucl                  \
            -title $CHR || exit 1
fi

for size in 1 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150
do
    echo "Preparing database of $size chromosomes"

    # set up directory to hold databases
    benchdir=$OUT/$size
    mkdir -p $benchdir

    # clean out
    rm -rf $benchdir/*

    if [[ -e $IN/raw ]];
        then
            mkdir $benchdir/raw || exit
            # select some random fasta files
            fasta_files=`ls $IN/raw/*.fa | sort -R --random-source=/dev/zero | head -n $size`
            for f in $fasta_files
            do
                echo "Copying $f"
                cp $f $benchdir/raw

                echo "Replacing defline with filename"
                sed -i '1s/.*/>'$(basename $f)'/' $f

                # echo "Adding newline at end of each file"
                # sed -i -e '$a\' $f

                echo #done...
            done
        else
            echo "No raw fasta files found".
            exit 1
    fi

    # concatenate raw fasta files
    cat $benchdir/raw/*.fa > $benchdir/$COMBINED

    # sanity check
    [[ -e $benchdir/$COMBINED ]] || exit 1

    echo "Making full blast DB"
    mkdir $benchdir/fulldb || exit
    $MAKEBLASTDB -in $benchdir/$COMBINED \
        -out $benchdir/fulldb/all        \
        -dbtype nucl                     \
        -title $CHR-$size || exit 1

    echo "Generating CaBLAST DB"
    mkdir $benchdir/cablastdb || exit
    $CABLAST/gen_compressed_db          \
        $benchdir/cablastdb/uniques.fa  \
        $benchdir/cablastdb/links.dat   \
        $benchdir/cablastdb/edits.dat   \
        smushed < $benchdir/$COMBINED | \
            tee $benchdir/cablastdb/gencomp_stdout.log

    echo "Preparing uniques blast DB"
    $MAKEBLASTDB -in $benchdir/cablastdb/uniques.fa \
        -out $benchdir/cablastdb/uniques            \
        -dbtype nucl                                \
        -title $CHR-cablast-uniques-$size || exit 1

    echo "Copying over CLASS DB"
    mkdir $benchdir/classdb
    for f in $fasta_files
        do
            bvfile=`echo $f | grep -o -e 'NA.*.[12]' -e 'HG.*.[12]' | sed 's/\./-/'`
            echo $bvfile
            cp $IN/bv/$bvfile $benchdir/classdb || exit 1
    done
    for f in `ls $IN/bv/vt*`
        do
            cp $f $benchdir/classdb
    done

done
