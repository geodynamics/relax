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

psxy -O -K -JX -R$bds -P -m \
         -W0.5p/10/10/10 -St0.2c -G10/10/10 \
        <<EOF >> $PSFILE
`awk '{print $2,$1}' $WDIR/opts.dat`
EOF

# name of observation points
pstext -O -K -JX -R$bds -P \
        -G0/0/0 -D0/0.1i \
        <<EOF >> $PSFILE
`awk '{print $2,$1," 08 0 4 CM ",$4}' $WDIR/opts.dat`
EOF

# cumulative displacement vectors (coseismic+relaxation)
for i in `awk '{print $4}' $WDIR/opts.dat | xargs `; do
	DISP=`grep -v "#" $WDIR/${i}.txt | awk -v i=$INDEX '{if ((NR-1)==i){print $3,$2}}'`
	COOR=`grep "$i" $WDIR/opts.dat | awk '{print $2,$1}'`
	# echo $i $COOR $DISP
psvelo -O -K -JX -R$bds -P \
        -A0.1/0.18/0.15 -W0.5px,0/0/0 \
        -G200/100/0 -Se${oscale}/0.91/06 \
        <<EOF >> $PSFILE
`echo $COOR $DISP 0 0 0`
EOF
done

# legend
psvelo -O -K -JX -R$bds -A0.1/0.18/0.15 \
	-W0.5px,0/0/0 -G200/100/0 \
	-Se${oscale}/0.91/06 \
<<EOF >> $PSFILE
-100 60 0.03 0 0 0 0
EOF

pstext  -K -O -R$bds -JX -P -G0/0/0 \
  <<EOF >>$PSFILE
-100 60 10 0 0.5 RM model (3 cm)
EOF
