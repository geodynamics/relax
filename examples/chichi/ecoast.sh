#!/bin/sh

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if [ "$#" -lt "2" ]; then
	echo $self": expected $self file.ps xmin/xmax/ymin/ymax. exiting."
	exit 1
fi

echo $self: $*
PSFILE=$1
bds=$2

WDIR=./gmt

#psxy -O -K -JX -R$bds -P -M \
#         -W1.5p/050/050/050  \
#       $WDIR/twcoasts_km.xy >> $PSFILE


psxy -O -K -JX -R$bds -P -M \
         -W0.8p/050/050/050  \
       $WDIR/coasts_km.xyz>> $PSFILE
