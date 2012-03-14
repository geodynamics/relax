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
iscale=$3
oscale=`echo $iscale | awk -F "c" '{print 1/$1}'`

WDIR=$(dirname $PSFILE)
INDEX=`echo $(basename $PSFILE) | awk '{print substr($0,0,3)}'`
DDIR=./gps

psvelo -O -K -JX -R$bds -P \
        -A0.05/0.3/0.15 -W0.5px,0/0/0 \
        -G158/33/106 -Se${oscale}/0.91/06 \
        <<EOF >> $PSFILE
`awk '{print $2,$1,$3/1e3,$4/1e3," 0 0 0"}' $DDIR/postseismic_neu_today_km.dat`
EOF

psxy -O -K -JX -R$bds -P -D0.4c/0 \
        -Sc0.4c -W0.5px,0/0/0 -C$WDIR/palette.cpt \
        <<EOF >> $PSFILE
`awk '{print $2,$1,$5/1e3}' $DDIR/postseismic_neu_today_km.dat`
EOF

if [ "1" == "0" ]; then
pstext -O -K -JX -R$bds -P \
        -G0/0/0 -D0/0.1i \
        <<EOF >> $PSFILE
`awk '{print $2,$1," 08 0 4 CM ",$4}' $WDIR/opts.dat`
EOF
fi

psvelo -O -K -JX -R$bds -P -A0.05/0.3/0.15 -W0.5px,0/0/0 -G158/33/106 -Se${oscale}/0.91/06 \
<<EOF >> $PSFILE
-90 145 0.2 0 0 0 0
EOF

pstext  -K -O -R$bds -JX -P -G0/0/0 \
  <<EOF >>$PSFILE
-80 150 10 0 0.5 CM data
EOF
