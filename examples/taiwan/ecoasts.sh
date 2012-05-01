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

WDIR=gmt

psxy -O -K -JX -R$bds -P -m \
         -W0.5p/10/10/10 \
        <<EOF >> $PSFILE
`cat $WDIR/twcoasts_km.xyz`
EOF

#psxy -O -K -JX -R$bds -P -m \
#         -W0.5p/10/10/10 -St0.2c \
#        <<EOF >> $PSFILE
#`awk '{print $4,$3}' $WDIR/cont_gps_km.dat`
#EOF

#pstext -O -K  -R$bds -JX -P \
#        -G0/0/0 -D0/-0.1 \
#        <<EOF >> $PSFILE
#`awk '{print $4,$3," 12 0 4 CT ",$2}' $WDIR/cont_gps_km.dat`
#EOF

#pstext -O -K  -R$bds -JX -P \
#        -G0/0/0 -D0/-0.1 \
#        <<EOF >> $PSFILE
#`cat landers_km.flt | awk 'BEGIN{c=0;}{if ($5 == 0){c=c+1;print $4,$3," 12 0 4 CT ",c}}'`
#EOF


