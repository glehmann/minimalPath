#!/bin/sh

export SEG=$HOME/Projects/ITK/minimalPath/artery.sh
export TARGDIR=/tmp/carotidL

cat $1 | (
    while read raw mark ; 
    do 
    #pref=$(dirname $mark)
    pref=${mark%%.*}
    pref=$(basename $pref)
    output=$TARGDIR/$pref
    echo $output
    $SEG $raw $mark $output
    done
)



