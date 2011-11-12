#!/bin/sh

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

while getopts "i:h:" flag
do
  case "$flag" in
    i) iset=1;inc=`echo "" | awk -v a=$OPTARG '{print a*3.14159265/180.0}'`;;
    h) hset=1;hdg=`echo "" | awk -v a=$OPTARG '{pring a*3.14159265/180.0}'`;;
  esac
done
for item in $iset $hset;do
	shift;shift
done
if (test $# -lt "1"); then
        echo "$self projects a set of grd files in the InSAR"
	echo "line-of-sight direction (all angles in degrees)."
	echo ""
	echo "usage: $self -i incidence -h heading index"
	echo "computes out.grd =  cos(incidence) * in-up.grd"
	echo "                   +sin(incidence) * sin(heading) * in-north.grd"
	echo "                   -sin(incidence) * cos(heading) * in-east.grd"
	echo ""
	echo "for ENVISAT, assume a constant incidence angle of"
	echo "23 deg. For descending, use a heading of 166 deg.;"
	echo "for ascending, use a heading of -13 deg."
	exit
	exit
fi
if [ "$iset" != "1" ]; then
	echo "$self: using default incidence angle of 23 deg."
	inc=`echo "" | awk '{print 23*3.14159265/180.0}'`
fi
if [ "$hset" != "1" ]; then
	echo "$self: using default heading of 166 deg. (ENVISAT descending)"
	hdg=`echo "" | awk '{print -166*3.14159265/180.0}'`
fi

while [ "$#" != "0" ];do

INDEX=$(basename $1)

WDIR=$(dirname $1)
echo $self": projecting index $INDEX in $WDIR in the radar line of sight"

GRDE=$WDIR/$INDEX-east.grd
GRDN=$WDIR/$INDEX-north.grd
GRDU=$WDIR/$INDEX-up.grd
GRDL=$WDIR/$INDEX-los.grd

ch=`echo "" | awk -v a=$hdg '{print cos(a) }'`
sh=`echo "" | awk -v a=$hdg '{print sin(a) }'`
ci=`echo "" | awk -v a=$inc '{print cos(a) }'`
si=`echo "" | awk -v a=$inc '{print sin(a) }'`

grdmath $ci $GRDU MUL $si $sh MUL $GRDN MUL ADD -$si $ch MUL $GRDE MUL ADD = $GRDL

shift
if [ "$#" != "0" ];then
	echo ""
fi
done


