#!/bin/sh


TMP=`mktemp -td carotid.XXXX`

trap  'rm -rf $TMP; exit 0' 0 1 2 3 14 15
echo $TMP

LABPROG=$HOME/Projects/ITK/connectedComponent/testLabelling
SEGPROG=$HOME/Projects/ITK/minimalPath/segArtery2
STLPROG=$HOME/Projects/ITK/meshing/SurfaceExtraction2
SNAP=/usr/local/bin/InsightSNAP

RM=/bin/rm

#CONV=/usr/local/bin/mnc2nii
#CONV=/usr/local/minc2/bin/mnc2nii
CONV=mnc2nii
RESMARK=$TMP/marker.mnc
LAB=$TMP/lab.nii
MARK=$TMP/mark.nii
CONTROL=$TMP/control.nii

SEGOUT=$TMP/result.nii.gz
CONTOUT=$TMP/sub_control.nii.gz
# $1 will be the original image, $2 the marker


if [ -e $2 ]; then
   # resample the marker to be like the input
   mincresample -nearest_neighbour -keep_real_range -like $1 $2 $RESMARK
   $CONV -byte $RESMARK $MARK
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

#$LABPROG $MARK $LAB 0 3 4 1 1
$SEGPROG -i $CONTROL -m $MARK -r 10 -l 1 -o $SEGOUT -s $CONTOUT 1 2 3 4 0 2 5 6

#if [ -e $CONTOUT ]; then 
#    echo "Done segmenting"
#    $SNAP --grey $CONTOUT --segmentation $SEGOUT & 
#    $SNAP --grey $CONTOUT --segmentation startmark.nii.gz
#else
#    echo "Some sort of failure - check error messages"
#fi

#$RM $LAB $MARK $CONTROL $SEGOUT $CONTOUT

#$RM $LAB $MARK $CONTROL

OUTPREF=$3

#$STLPROG $SEGOUT 1 $OUTPREF".stl"

mv $SEGOUT $OUTPREF"_seg.nii.gz"
mv $CONTOUT $OUTPREF"_raw.nii.gz"
mv /tmp/marker.nii.gz $OUTPREF"_mark.nii.gz"
mv /tmp/grad.nii.gz $OUTPREF"_grad.nii.gz"
mv /tmp/cost.nii.gz $OUTPREF"_cost.nii.gz"