#!/bin/sh
. /etc/profile
. $HOME/.bash_profile


export SEG=$HOME/Projects/ITK/minimalPath/artery.sh
export TARGDIR=/$HOME/Projects/carotidL
#export TARGDIR=/$HOME/Projects/carotidR
TARGETDIR=$HOME/ScriptOutputs

TMPFILE=`mktemp` || exit 1
trap 'rm -f $TMPFILE; exit 0' 0 1 2 3 14 15

echo $TMPFILE

function submitSeg()
{
  # raw image is $1
  # marker image is $2
  # $3 is the output image

jname=`basename $2`

cat >$TMPFILE <<-EOF
#!/bin/sh
#\$ -S /bin/sh
#\$ -o $TARGETDIR/$jname".stdout"
#\$ -e $TARGETDIR/$jname".stderr"
. /etc/profile
source $HOME/.bash_profile
$SEG $1 $2 $3
EOF

qsub -N $jname -l h_rt=00:15:00 $TMPFILE
#qsub -N $jname -q bque $TMPFILE

}

cat $1 | (
    while read raw mark ; 
    do 
    #pref=$(dirname $mark)
    pref=${mark%%.*}
    pref=$(basename $pref)
    output=$TARGDIR/$pref
    #echo $output
    submitSeg $raw $mark $output
    done
)



