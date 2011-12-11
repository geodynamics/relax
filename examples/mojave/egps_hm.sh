#!/bin/sh

# this script is called with grdmap.sh using 
#
#   grdmap.sh -e ./egps_hm.sh file.grd
#

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if [ "$#" -lt "5" ]; then
	echo $self": expected $self file.ps xmin/xmax/ymin/ymax scale u3 height. exiting."
	exit 1
fi

echo $self: $*
PSFILE=$1
bds=$2
iscale=$3
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

