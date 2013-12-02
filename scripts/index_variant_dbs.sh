#a!/bin/bash

DBDIR=/data/apagan/vdb
TABIX_ROOT=/home/apagan/code/class/vcflib/tabixpp
DBFILE="vt.db"

for in_path in `ls $DBDIR`
do
    chrnum=`echo $(basename $in_path) | egrep -o "[[:digit:]]*"`
    echo "indexing chromosome $chrnum"
    db_uncompressed=$DBDIR/$in_path/$DBFILE
    db_extended=$db_uncompressed".chradd"
    db_compressed=$db_uncompressed".gz"
    awk '{print "'$chrnum'\t"$0}' $db_uncompressed > $db_extended
    $TABIX_ROOT/bgzip -c $db_extended > $db_compressed \
        && rm $db_extended \
        && $TABIX_ROOT/tabix -s 1 -b 4 -e 4 $db_compressed &
done
wait