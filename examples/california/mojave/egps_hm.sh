#!/bin/bash

# this script is called with grdmap.sh using 
#
#   grdmap.sh -e ./egps_hm.sh file.grd
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
	echo ${self} overlays Hector Mine GPS offsets on a GMT map.
	echo
	echo usage: $self -b xmin/xmax/ymin/ymax file.ps
	exit 1
fi

echo $self: $cmdline
PSFILE=$1
iscale=$VECTOR
oscale=`echo $iscale | awk -F "c" '{print 1/$1}'`
WDIR=$(dirname $PSFILE)
DDIR=gps

# gps coseismic offsets for Hector Mine
psvelo -O -K -JX -R$bds -P \
        -A0.1/0.18/0.15 \
        -G20/20/20 -Se${oscale}/0.91/06 \
        <<EOF >> $PSFILE
`awk '{if (1<NR){print $2,$3,$4,$5}}' $DDIR/HM_coseis.dat | \
	proj +proj=utm +zone=11 | \
	awk '{print ($1-566940.91)/1e3,($2-3828373.73)/1e3,$3/1e2,$4/1e2," 0 0 0"}'`
EOF
# -A arrow width, length, width
#        -A0.03/0.18/0.09 \
# legend: -116.2 33.55 0.05 0 0 0 0 5 cm

