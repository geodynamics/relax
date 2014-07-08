#!/bin/bash

# extra script for mapping tool grdmap.sh 
# overlays the contours of fault patches projected in map view
#
#   grdmap.sh -e erpatch.sh WDIR/file.grd
#
# assumes files WDIR/rfaults-???.dat contain slip model, as
# created by program relax and alike.

set -e
self=$(basename $0)
selfdir=$(dirname $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

while getopts "b:c:gp:v:H:J" flag
do
	case "$flag" in
	b) bset=1;bds=$OPTARG;;
	c) cset=1;carg=$OPTARG;;
	g) gset=1;;
	p) pset=1;U3=$OPTARG;;
	v) vset=1;SIZE=$OPTARG;VECTOR=$OPTARG"c";;
	J) Jset="-J";;
	H) Hset=1;HEIGHT=$OPTARG;;
	esac
done
for item in $bset $cset $pset $vset $Hset; do
	shift;shift
done
for item in $gset $Jset; do
	shift
done

if [ "$#" -lt "1" ]; then
	echo ${self} overlays the contours of fault patches on a GMT map.
	echo
	echo usage: $self -b xmin/xmax/ymin/ymax file.ps
	exit 1
fi

echo $self: $cmdline
PSFILE=$1
WDIR=$(dirname $PSFILE)

# plot a coseismic model (do nothing in geographic coordinates)
if [ "$gset" != "1" ]; then
	if [ -e "$WDIR/rfaults-001.xy" ]; then
		psxy -O -K -JX -R$bds -L -M \
			-W1.0p/20/20/20 \
			<<EOF >> $PSFILE
`awk '{if (">"==$1){print $1}else{print $1,$2}}' $WDIR/rfaults-???.xy`
EOF
	fi
	if [ -e "$WDIR/rdykes-001.xy" ]; then
		psxy -O -K -JX -R$bds -L -M \
			-W1.0p/20/20/20 \
			<<EOF >> $PSFILE
`awk '{if (">"==$1){print $1}else{print $1,$2}}' $WDIR/rdykes-???.xy`
EOF
	fi
fi


