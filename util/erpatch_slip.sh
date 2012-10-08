#!/bin/bash

# extra script for mapping tool grdmap.sh 
# overlays the slip model projected in map view
#
#   grdmap.sh -e erpatch_slip.sh WDIR/file.grd
#
# assumes files WDIR/rfaults-???.dat contain slip model, as
# created by program relax and alike.

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
	echo ${self} overlays the slip model on a GMT map.
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
		psxy -O -K -JX -R$bds -C$WDIR/palette.cpt -W1p0/0/0 -L -M \
	        	 -W0.5p/20/20/20 \
	        	<<EOF >> $PSFILE
`awk '{if (">"==$1){print $0}else{print $2,$1}}' $WDIR/rfaults-???.xy`
EOF
	fi
	if [ -e "$WDIR/rdykes-001.xy" ]; then
		psxy -O -K -JX -R$bds -C$WDIR/palette.cpt -W1p0/0/0 -L -M \
	        	 -W0.5p/20/20/20 \
	        	<<EOF >> $PSFILE
`awk '{if (">"==$1){print $0}else{print $2,$1}}' $WDIR/rdykes-???.xy`
	fi
fi

