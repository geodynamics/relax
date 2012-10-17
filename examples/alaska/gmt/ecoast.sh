#!/bin/bash

# extra script for mapping tool grdmap.sh
# overlays coastlines.
#
#   grdmap.sh -e ecoast.sh WDIR/file.grd
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
	echo ${self} overlays coastlines on a GMT map.
	echo
	echo usage: $self -b xmin/xmax/ymin/ymax file.ps
	exit 1
fi

echo $self: $cmdline
PSFILE=$1
WDIR=./gmt

#psxy -O -K -JX -R$bds -P -M \
#         -W1.5p/050/050/050  \
#       $WDIR/twcoasts_km.xy >> $PSFILE

if [ "$gset" == "" ]; then
	psxy -O -K -JX -R$bds -P -M \
		-W0.8p/050/050/050  \
		$WDIR/coasts_km.xyz>> $PSFILE
fi


