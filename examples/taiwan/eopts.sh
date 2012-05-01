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

psxy -O -K -JX -R$bds -P \
         -W0.5p/10/10/10 -St0.1c -G10/10/10 \
        <<EOF >> $PSFILE
`awk '{print $2,$1}' $WDIR/opts.dat`
EOF

#if [ "0" == "0" ]; then
#pstext -O -K -JX -R$bds -P \
#        -G0/0/0 -D0/0.1i \
#        <<EOF >> $PSFILE
#`awk '{print $2,$1," 08 0 4 CM ",$4}' $WDIR/opts.dat`
#EOF
#fi

for i in `awk '{print $4}' $WDIR/opts.dat | xargs `; do
	DISP=`grep -v "#" $WDIR/${i}-relax.txt | awk -v i=$INDEX '{if (NR==i){print $3,$2}}'`
	COOR=`grep "$i" $WDIR/opts.dat | awk '{print $2,$1}'`
	#echo $i $COOR $DISP
psvelo -O -K -JX -R$bds  \
        -A0.05/0.3/0.15 -W0.5px,0/0/0 \
        -G200/100/0 -Se${oscale}/0.91/06 \
        <<EOF >> $PSFILE
`echo $COOR $DISP 0 0 0`
EOF
done

psvelo -O -K -JX -R$bds -A0.05/0.3/0.15 -W0.5px,0/0/0 -G200/100/0 -Se${oscale}/0.91/06 \
<<EOF >> $PSFILE
-90 160 0.2 0 0 0 0
EOF

pstext  -K -O -R$bds -JX -P -G0/0/0 \
  <<EOF >>$PSFILE
-80 165 10 0 0.5 CM model
EOF
