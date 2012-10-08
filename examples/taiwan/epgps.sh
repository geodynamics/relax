#!/bin/bash

# this script is called with grdmap.sh using 
#
#   grdmap.sh -e ./ecoast.sh file.grd
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
	echo ${self} overlays coast lines on a GMT map.
	echo
	echo usage: $self -b xmin/xmax/ymin/ymax file.ps
	exit 1
fi

echo $self: $cmdline
PSFILE=$1
iscale=$VECTOR
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
