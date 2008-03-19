#!/bin/sh


TMP=`mktemp -td jugular.XXXX`

#trap  'rm -rf $TMP; exit 0' 0 1 2 3 14 15
echo $TMP

#BUILD32=Build32
#SEGPROG=/home/richardb/Projects/ITK/minimalPath/${BUILD32}/segJugular
SEGPROG=/tmp/Build/segJugular
SNAP=/usr/local/bin/InsightSNAP

RM=/bin/rm

#CONV=/usr/local/bin/mnc2nii
#CONV=/usr/local/minc2/bin/mnc2nii
CONV=mnc2nii
RESJUGMARK=$TMP/jugmarker.mnc
JUGMARK=$TMP/jugmark.nii

CONTROL=$TMP/control.nii
RERAW=$TMP/reorient.nii
# $1 will be the original image, $2 the marker


if [ -e $2 ]; then
   # resample the marker to be like the input
   mincresample -nearest_neighbour -keep_real_range -like $1 $2 $RESJUGMARK
   $CONV -byte $RESJUGMARK $JUGMARK
else
  echo "No such file $2"
  exit
fi

if [ -e $1 ]; then
    $CONV  $1 $CONTROL
else
  echo "No such file $1"
  exit
fi

echo "------------------------------------"
echo "Finished converting"

SEGOUT=$TMP/result.nii.gz

$SEGPROG -d $TMP -i $CONTROL  -j $JUGMARK -o $SEGOUT -r $RERAW

dd=$(dirname $1)
bb=$(basename $dd)

if [ -e $SEGOUT ]; then 
    echo "Done segmenting"
    $SNAP --grey $RERAW --segmentation $SEGOUT  
else
    echo "Some sort of failure - check error messages"
fi

#$RM $LAB $MARK $CONTROL $SEGOUT $CONTOUT

#$RM $LAB $MARK $CONTROL

#OUTPREF=$3

#$STLPROG $SEGOUT 1 $OUTPREF".stl"


#mv $SEGOUT $OUTPREF"_seg.nii.gz"
#mv $CONTOUT $OUTPREF"_raw.nii.gz"
#mv $TMP/marker.nii.gz $OUTPREF"_mark.nii.gz"
#mv $TMP/grad.nii.gz $OUTPREF"_grad.nii.gz"
#mv $TMP/cost.nii.gz $OUTPREF"_cost.nii.gz"
#mv $TMP/points.nii.gz $OUTPREF"_points.nii.gz"
