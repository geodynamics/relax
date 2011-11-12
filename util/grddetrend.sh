#!/bin/sh

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit 1' ERR

if [ "$#" -lt "1" ]; then
	echo usage: $self file.grd 
	exit 1
fi
INDEX=$(basename $1 .grd)

GRD=${INDEX}.grd
XYZ=${INDEX}.xyz
XYZW=${INDEX}.xyzw
XYR=${INDEX}.xyr
OUT=${INDEX}_detrended.grd

echo $self: convert $GRD to $XYZ
grd2xyz $GRD -S > $XYZ

echo $self: detrends $(basename $GRD) to $(basename $XYR)
# with weights
#awk '{if ($2>35.75){w=0}else{w=1};print $0,w}' $XYZ > $XYZW
#trend2d $XYZW -Fxyr -N4r -W > $XYR
# no weights
trend2d $XYZ -Fxyr -N4 > $XYR

# resample detrended file to GRD format
echo $self: converts $(basename $XYR) to $(basename $OUT)
bds=`grdinfo -C $GRD | awk '{print $2"/"$3"/"$4"/"$5}'`
inc=`grdinfo -C $GRD | awk '{print $8"/"$9}'`
xyz2grd $XYR -G$OUT -I$inc -R$bds
echo $self: $(basename $OUT) created successfully

