#!/bin/sh

set -e

self=$(basename $0)


if [ $# -lt "2" ]; then
	echo usage: $self profile.xy index
	echo
	echo creates the .xy files index-prof-{north,east,up}.xy
	echo containing the data across profile.xy
	echo
	exit
fi

PROF=$1
shift

while [ "$#" != "0" ]; do

WDIR=$(dirname $1)
INDEX=$(basename $1)

echo $self: using directory $WDIR

U1=$WDIR/$INDEX-north.grd
U2=$WDIR/$INDEX-east.grd
U3=$WDIR/$INDEX-up.grd

P1=$WDIR/$INDEX-prof-north.xy
P2=$WDIR/$INDEX-prof-east.xy
P3=$WDIR/$INDEX-prof-up.xy

grdtrack $PROF -G$U1 | awk '{print $1,$3}' > $P1
grdtrack $PROF -G$U2 | awk '{print $1,$3}' > $P2
grdtrack $PROF -G$U3 | awk '{print $1,$3}' > $P3

echo $self: Created $P1 $P2 $P3

shift
done
