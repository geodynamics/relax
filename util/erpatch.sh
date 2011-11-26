#!/bin/sh

# extra script for mapping tool map.sh plotting the contours
# of fault patches projected in map view
#
#   map.sh -e rpatch.sh WDIR/file.grd
#
# assumes files WDIR/rfaults-???.dat contain slip model, as
# created by program relax and alike.

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if [ "$#" -lt "2" ]; then
	echo $self": extra script for map.sh."
	echo $self": expected $self file.ps xmin/xmax/ymin/ymax. exiting."
	exit 1
fi

echo $self: $*
PSFILE=$1
bds=$2

WDIR=$(dirname $PSFILE)

# plot a coseismic model
psxy -O -K -JX -R$bds -L -M \
         -W1.0p/20/20/20 \
        <<EOF >> $PSFILE
`awk '{if (">"==$1){print $1}else{print $2,$1}}' $WDIR/rfaults-???.xy`
EOF

