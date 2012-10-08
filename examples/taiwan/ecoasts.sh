#!/bin/bash

# this script is called with grdmap.sh using 
#
#   grdmap.sh -e ./ecoasts.sh file.grd
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


