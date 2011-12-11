#!/bin/sh

# this script is called with grdmap.sh using 
#
#   grdmap.sh -e ./eopts.sh file.grd
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

psxy -O -K -JX -R$bds -P -m \
         -W0.5p/10/10/10 -St0.2c -G10/10/10 \
        <<EOF >> $PSFILE
`awk '{print $2,$1}' $WDIR/opts.dat`
EOF

if [ "1" != "0" ]; then
pstext -O -K -JX -R$bds -P \
        -G0/0/0 -D0/0.1i \
        <<EOF >> $PSFILE
`awk '{print $2,$1," 08 0 4 CM ",$4}' $WDIR/opts.dat`
EOF
fi

for i in `awk '{print $4}' $WDIR/opts.dat | xargs `; do
	DISP=`grep -v "#" $WDIR/$i.txt | awk '{print $3,$2}'`
	COOR=`grep "$i" $WDIR/opts.dat | awk '{print $2,$1}'`
	echo $i $COOR $DISP
psvelo -O -K -JX -R$bds -P \
        -A0.1/0.18/0.15 -W0.5px,0/0/0 \
        -G220/220/220 -Se${oscale}/0.91/06 \
        <<EOF >> $PSFILE
`echo $COOR $DISP 0 0 0`
EOF
done

