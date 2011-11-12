#!/bin/sh

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

my_list_first()
{
if [ $# -eq 0 ]; then
    echo "";exit 0
fi

echo $1
return 0
}

my_list_last()
{
if [ $# -eq 0 ]; then
    echo "";exit 0
fi

while [ $# -ne 1 ]; do
  shift
done
echo $1
return 0
}

my_list_remove_head()
{
shift
echo $*
return 0
}

my_list_count()
{
count=0
while [ $# -ne 0 ]; do
  shift
  count=`echo "$count + 1"|bc`
done
echo $count
return 0
}


if (test $# -ne "3"); then
        echo "usage: "$self" WDIR/INDEX X Y"
	echo "       creates time series WDIR/INDEX-X-Y.xytneu for coordinates (X,Y)"
	echo "       from files ???-relax-{east,north,up}.grd in directory WDIR."
	echo "       requires file WDIR/time.txt to add time information."
	exit 1
fi


WDIR=$(dirname $1)
TSEFILE=$WDIR/$(basename $1 .xytneu)-$2-$3.xytneu
COORD=$2" "$3

echo $self": Using directory "$WDIR
echo $self": Creating file "$TSEFILE
echo $self": Sampling at coordinates "$COORD

LIST=`ls $WDIR/???-relax-up.grd | awk -F "[/-]" '{print $(NF-2)}'`

TMPFILE=$WDIR/temp
rm -f $TSEFILE
rm -f $TMPFILE"C"
rm -f $TMPFILE"E"
rm -f $TMPFILE"N"
rm -f $TMPFILE"U"

for ITEM in $LIST
do

GRDFILEE=$WDIR/$ITEM-relax-east.grd
GRDFILEN=$WDIR/$ITEM-relax-north.grd
GRDFILEU=$WDIR/$ITEM-relax-up.grd

echo $self": Processing item "$ITEM" at coordinates "$COORD

echo "$COORD" >> $TMPFILE"C"

grdtrack -G$GRDFILEE -Z <<EOF >> $TMPFILE"E"
$COORD
EOF
grdtrack -G$GRDFILEN -Z <<EOF >> $TMPFILE"N"
$COORD
EOF
grdtrack -G$GRDFILEU -Z <<EOF >> $TMPFILE"U"
$COORD
EOF

done
TIMES=`cat $WDIR/time.txt`
N1=`my_list_count $TIMES`
N2=`my_list_count $LIST`

echo $self": Done extracting from "$N2" files."

if (test $N1 -eq $N2); then
	paste $TMPFILE"C" $WDIR/time.txt $TMPFILE"N" $TMPFILE"E" $TMPFILE"U" >$TSEFILE
	echo $self": File "$TSEFILE" successfully processed with time info"
else
	paste $TMPFILE"C" $TMPFILE"N" $TMPFILE"E" $TMPFILE"U" >$TSEFILE
	echo $self": Time mismatch; file "$TSEFILE" has no time info"
fi
rm -f $TMPFILE"C" $TMPFILE"N" $TMPFILE"E" $TMPFILE"U"

