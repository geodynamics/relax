#!/bin/sh

# extra script for mapping tool map.sh plotting the slip model
# projected in map view
#
#   map.sh -e rpatch_slip.sh WDIR/file.grd
#
# assumes files WDIR/rfaults-???.dat contain slip model, as
# created by program relax and alike.

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

psxy -O -K -JX -R$bds -C$WDIR/palette.cpt -W1p0/0/0 -L -M \
         -W0.5p/20/20/20 \
        <<EOF >> $PSFILE
`awk '{if (">"==$1){print $0}else{print $2,$1}}' $WDIR/rfaults-???.xy`
EOF

