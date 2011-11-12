#!/bin/sh

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

while getopts "d:s:" flag
do
  case "$flag" in
    d) dset=1;dip=`echo "" | awk -v a=$OPTARG '{print a*3.1415926535897932385/180.0}'`;;
    s) sset=1;str=`echo "" | awk -v a=$OPTARG '{print a*3.1415926535897932385/180.0}'`;;
  esac
done
for item in $dset $sset;do
	shift;shift
done

if (test $# -lt "1"); then
        echo "$self projects a stress tensor along a hypothetic"
	echo "fault defined by strike and dip orientation"
	echo ""
	echo "usage: $self -s strike -d dip index"
	echo ""
	echo "a fault normal, dip and strike vectors are "
	echo "n = [+sin(str) * sin(dip), -cos(str) * sin(dip), -cos(dip)]"
	echo "d = [+sin(str) * cos(dip), -cos(str) * cos(dip), +sin(dip)]"
	echo "s = [cos(str), sin(str), 0]"
	echo
	echo "the local traction components are"
	echo "t1 = Sxx * n(1) + Sxy * n(2) + Sxz * n(3)"
	echo "t2 = Sxy * n(1) + Syy * n(2) + Syz * n(3)"
	echo "t3 = Sxz * n(1) + Syz * n(2) + Szz * n(3)"
	echo
	echo "the strike, normal and dip components are"
	echo "tn = t1 * n(1) + t2 * n(2) + t3 * n(3)"
	echo "ts = (t1 - tn * n(1)) * s(1) + (t2 - tn * n(2)) * s(2)"
	echo "td = (t1 - tn * n(1)) * d(1) + (t2 - tn * n(2)) * d(2) + (t3 - tn * n(3)) * d(3)"
	echo ""
	echo "by convention, ts is positive for left lateral strike shear."
	echo "the 1, 2 & 3 directions are north, east and down, respectively."
	echo ""
	exit
fi
if [ "$dset" != "1" ]; then
	echo "$self: using default dip angle of 90 degrees."
	dip=1.57079632679489661923132169164
else
	echo "$self: using dip angle "$dip "radian."
fi
if [ "$sset" != "1" ]; then
	echo "$self: using default strike angle of 0 degrees."
	str=0
else
	echo "$self: using strike angle "$str "radian."
fi

while [ "$#" != "0" ];do

INDEX=$(basename $1)

WDIR=$(dirname $1)
echo $self": projecting index $INDEX in $WDIR in traction"

# input stress components
GRD11=$WDIR/$INDEX-s11.grd
GRD12=$WDIR/$INDEX-s12.grd
GRD13=$WDIR/$INDEX-s13.grd
GRD22=$WDIR/$INDEX-s22.grd
GRD23=$WDIR/$INDEX-s23.grd
GRD33=$WDIR/$INDEX-s33.grd

# output traction components
GRDTS=$WDIR/$INDEX-ts.grd
GRDTD=$WDIR/$INDEX-td.grd
GRDTN=$WDIR/$INDEX-tn.grd

# temporary traction components
T1=$WDIR/$INDEX-t1.grd
T2=$WDIR/$INDEX-t2.grd
T3=$WDIR/$INDEX-t3.grd

# cosine and sine of angles
s1=`echo "" | awk -v d=$dip -v s=$str '{print cos(s) }'`
s2=`echo "" | awk -v d=$dip -v s=$str '{print sin(s) }'`
n1=`echo "" | awk -v d=$dip -v s=$str '{print +sin(s)*sin(d) }'`
n2=`echo "" | awk -v d=$dip -v s=$str '{print -cos(s)*sin(d) }'`
n3=`echo "" | awk -v d=$dip -v s=$str '{print -cos(d) }'`
d1=`echo "" | awk -v d=$dip -v s=$str '{print +sin(s)*cos(d) }'`
d2=`echo "" | awk -v d=$dip -v s=$str '{print -cos(s)*cos(d) }'`
d3=`echo "" | awk -v d=$dip -v s=$str '{print +sin(d) }'`

grdmath $GRD11 $n1 MUL $GRD12 $n2 MUL ADD $GRD13 $n3 MUL ADD = $T1
grdmath $GRD12 $n1 MUL $GRD22 $n2 MUL ADD $GRD23 $n3 MUL ADD = $T2
grdmath $GRD13 $n1 MUL $GRD23 $n2 MUL ADD $GRD33 $n3 MUL ADD = $T3

# projection in fault-aligned reference system

# normal component of traction vectors 
# Tn=t1 * n(1) + t2 * n(2) + t3* n(3)
grdmath $T1 $n1 MUL $T2 $n2 MUL ADD $T3 $n3 MUL ADD = $GRDTN

# strike component of traction
# ts = (t1 - tn * n(1)) * s(1) + (t2 - tn * n(2)) * s(2)
grdmath $T1 $GRDTN $n1 MUL SUB $s1 MUL $T2 $GRDTN $n2 MUL SUB $s2 MUL ADD = $GRDTS

# dip component of traction
# td = (t1 - tn * n(1)) * d(1) + (t2 - tn * n(2)) * d(2) + (t3 - tn * n(3)) * d(3)
grdmath $T1 $GRDTN $n1 MUL SUB $d1 MUL $T2 $GRDTN $n2 MUL SUB $d2 MUL ADD $T3 $GRDTN $n3 MUL SUB $d3 MUL ADD = $GRDTD

shift
if [ "$#" != "0" ];then
	echo ""
fi
done


