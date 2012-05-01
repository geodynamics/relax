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

WDIR=$(dirname $PSFILE)
INDEX=`echo $(basename $PSFILE) | awk '{print substr($0,0,3)}'`
GRD=$WDIR/$(basename $PSFILE -plot.ps)

echo $GRD
grdcontour $GRD -O -K -C0.05 -JX -R$bds -P \
         -W0.3p/100/100/100 >> $PSFILE

