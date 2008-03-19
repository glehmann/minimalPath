#!/bin/sh


TMP=`mktemp -td carotid.XXXX`

#trap  'rm -rf $TMP; exit 0' 0 1 2 3 14 15
echo $TMP

LABPROG=$HOME/Projects/ITK/connectedComponent/testLabelling
SEGPROG=$HOME/Projects/ITK/minimalPath/segArtery2
#SEGPROG=/tmp/Build/segArtery2
STLPROG=$HOME/Projects/ITK/meshing/SurfaceExtraction2
SNAP=/usr/local/bin/InsightSNAP

RM=/bin/rm

#CONV=/usr/local/bin/mnc2nii
#CONV=/usr/local/minc2/bin/mnc2nii
CONV=mnc2nii
RESMARK=$TMP/marker.mnc
RESJUGMARK=$TMP/jugmarker.mnc
LAB=$TMP/lab.nii
MARK=$TMP/mark.nii
JUGMARK=$TMP/jugmark.nii

CONTROL=$TMP/control.nii

SEGOUT=$TMP/result.nii.gz
CONTOUT=$TMP/sub_control.nii.gz
# $1 will be the original image, $2 the marker


if [ -e $2 ]; then
   # resample the marker to be like the input
   mincresample -nearest_neighbour -keep_real_range -like $1 $2 $RESMARK
   mincresample -nearest_neighbour -keep_real_range -like $1 $3 $RESJUGMARK
   $CONV -byte $RESMARK $MARK
   $CONV -byte $RESJUGMARK $JUGMARK
else
  echo "No such file $2"
  exit
fi

if [ -e $1 ]; then
    $CONV -signed -short $1 $CONTROL
else
  echo "No such file $1"
  exit
fi

echo "------------------------------------"
echo "Finished converting"

#$LABPROG $MARK $LAB 0 3 4 1 1
$SEGPROG -d $TMP -i $CONTROL -m $MARK -j $JUGMARK -r 10 -l 1 -o $SEGOUT -s $CONTOUT 1 2 3 4 0 2 5 6

dd=$(dirname $1)
bb=$(basename $dd)

BNAME=/tmp/${bb}.nii.gz
MNAME=/tmp/${bb}_jug.nii.gz

mv /tmp/raw.nii.gz $BNAME
mv /tmp/mask.nii.gz $MNAME

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


#mv $SEGOUT $OUTPREF"_seg.nii.gz"
#mv $CONTOUT $OUTPREF"_raw.nii.gz"
#mv $TMP/marker.nii.gz $OUTPREF"_mark.nii.gz"
#mv $TMP/grad.nii.gz $OUTPREF"_grad.nii.gz"
#mv $TMP/cost.nii.gz $OUTPREF"_cost.nii.gz"
#mv $TMP/points.nii.gz $OUTPREF"_points.nii.gz"
