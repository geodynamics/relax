#!/bin/bash

# PURPOSE:
# 
# rotate the spatial coordinates of the three
# grd files index-north.grd index-east.grd and
# index-up.grd by a given angle.

set -e
self=$(basename $0)
trap 'echo ${self}: Some errors occurred. Exiting.; exit 1' ERR

rotate_horizontal(){

	# check if files exists
	if [ ! -e $IGRDN ]; then
		echo $self": could not find "$IGRDN
		exit 1
	fi
	if [ ! -e $IGRDE ]; then
		echo $self": could not find "$IGRDE
		exit 1
	fi
	
	XYZN=$WDIR/$(basename $IGRDN .grd).xyz
	XYZE=$WDIR/$(basename $IGRDE .grd).xyz
	bds=`grdinfo -C $IGRDN | awk '{print $2"/"$3"/"$4"/"$5}'`
	inc=`grdinfo -C $IGRDN | awk '{print $8"/"$9}'`
	rad=`grdinfo -C $IGRDN | awk '{print sqrt($8^2+$9^2)*1.5}'`
	
	grd2xyz $IGRDN > $XYZN
	grd2xyz $IGRDE > $XYZE
	
	echo $self": creating $(basename $OGRDE)"
	
	# nearest neighbor gridding
	paste $XYZN $XYZE | \
		awk -v v=$PHI '{x=cos(v)*$1-sin(v)*$2; \
				y=sin(v)*$1+cos(v)*$2; \
				n=cos(v)*$3+sin(v)*$6; \
				print x,y,n}' | \
			nearneighbor -G$OGRDN -R$bds -I$inc -S$rad

	echo $self": creating $(basename $OGRDN)"

	paste $XYZN $XYZE | \
		awk -v v=$PHI '{x=cos(v)*$1-sin(v)*$2; \
				y=sin(v)*$1+cos(v)*$2; \
				e=-sin(v)*$3+cos(v)*$6; \
				print x,y,e}' | \
			nearneighbor -G$OGRDE -R$bds -I$inc -S$rad
	
}

rotate(){

	# check if file exists
	if [ ! -e $IGRD ]; then
		echo $self": could not find "$IGRD
		exit 1
	fi

	XYZ=$WDIR/$(basename $IGRD .grd).xyz
	bds=`grdinfo -C $IGRD | awk '{print $2"/"$3"/"$4"/"$5}'`
	inc=`grdinfo -C $IGRD | awk '{print $8"/"$9}'`
	rad=`grdinfo -C $IGRD | awk '{print sqrt($8^2+$9^2)*1.5}'`
	grd2xyz $IGRD > $XYZ

	echo $self": creating $(basename $OGRD) "

	# nearest neighbor gridding
	awk -v v=$PHI '{x=cos(v)*$1-sin(v)*$2;y=sin(v)*$1+cos(v)*$2;print x,y,$3}' $XYZ | \
		nearneighbor -G$OGRD -R$bds -I$inc -S$rad

}

if [ $# -lt "2" ]; then
	echo ${self} angle index1 ... indexN
        echo rotates the grd files index-{north,east,up}.grd
	echo by a given angle.
        exit 1
fi

ANG=$1
# radian angle
PHI=`echo $ANG | awk '{print $1/180*3.14159265}'`
shift

# loop over the input files
while [ "$#" != "0" ];do

INDEX=$(basename $1)
WDIR=$(dirname $1)


IGRD=$WDIR/$INDEX
if [ -e $IGRD ]; then
	# single scalar-valued file of the form file.grd
	# transformed to file-rot.grd

	echo $self": rotating $(basename $IGRD) by $ANG degrees counter clockwise"

	OGRD=$WDIR/$(basename $INDEX .grd)-rot.grd
	rotate
else
	# vector-valued grid with three files of the form 
	# index-east.grd index-north.grd and index-up.grd 
	# converted to index-rot-up.grd, index-rot-east.grd
	# and index-rot-up.grd, with rotation of position
	# and horizontal vectors.

	echo $self": rotating $(basename $INDEX) by $ANG degrees counter clockwise"

	IGRD=$WDIR/$INDEX-up.grd
	OGRD=$WDIR/$INDEX-rot-up.grd
	if [ -e $IGRD ]; then
		rotate
	else
		echo $self": could not find "$(basename $IGRD)" in "$WDIR
		exit 1
	fi
	IGRDE=$WDIR/$INDEX-east.grd
	OGRDE=$WDIR/$INDEX-rot-east.grd
	IGRDN=$WDIR/$INDEX-north.grd
	OGRDN=$WDIR/$INDEX-rot-north.grd
	rotate_horizontal
fi

shift
if [ "$#" != "0" ];then
	echo ""
fi
done
