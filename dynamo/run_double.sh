#!/bin/bash
startfile=$1
for gzfilelong in /net/annapurna/data3/dynamo_toga_v2/cruise2/PUF/PPI/to1111*2.3.puf.gz
do
    cp $gzfilelong .
    gzfile=`basename $gzfilelong`
    echo $gzfile $startfile
    python double_unfolding_ufgz.py $gzfile $startfile
    rm $gzfile
    file=`echo $gzfile | sed -e 's/.gz$//'`
    startfile=$file.unfold.nc
done
