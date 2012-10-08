#!/bin/bash

# this script is called with grdmap.sh using 
#
#   grdmap.sh -e ./econtour.sh file.grd
#

set -e
self=$(basename $0)
selfdir=$(dirname $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

while getopts "b:c:gp:v:H:" flag
do
	case "$flag" in
	b) bset=1;bds=$OPTARG;;
	c) cset=1;carg=$OPTARG;;
	g) gset=1;;
	p) pset=1;U3=$OPTARG;;
	v) vset=1;SIZE=$OPTARG;VECTOR=$OPTARG"c";;
	H) Hset=1;HEIGHT=$OPTARG;;
	esac
done
for item in $bset $cset $pset $vset $Hset; do
	shift;shift
done
for item in $gset; do
	shift
done

if [ "$#" -lt "1" ]; then
	echo ${self} overlays contours of .grd file on a GMT map.
	echo
	echo usage: $self -b xmin/xmax/ymin/ymax file.ps
	exit 1
fi

echo $self: $cmdline
PSFILE=$1
WDIR=$(dirname $PSFILE)
INDEX=`echo $(basename $PSFILE) | awk '{print substr($0,0,3)}'`
GRD=$WDIR/$(basename $PSFILE -plot.ps)

grdcontour $GRD -O -K -C0.05 -JX -R$bds -P \
         -W0.3p/100/100/100 >> $PSFILE

