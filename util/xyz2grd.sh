#!/bin/sh

# author: sylvain barbot
#
# last modified: 06/17/2010
#
# xyz2grd.sh converts one or a set of three files with extensions
# -north.xyz -east.xyz and -up.xyz to GMT's binary file .grd.

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if (test $# -lt "1"); then
	echo "usage: "$self" [-b -10/10/-10/10] [-i 0.1] [-s 1] WDIR/INDEX1 WDIR/INDEX2 ... WDIR/INDEXN"
	echo ""
	echo "options:"
	echo "  -b sets bounds of the GRD files. default is span of .xyz file"
	echo "  -i sets the sampling interval of the GRD files. default is y interval in .xyz file"
	echo "  -s sets the scaling factor between XYZ and GRD files (GRD = s * XYZ). default is 1"
        exit
fi


while getopts "b:i:s:" flag
do
  case "$flag" in
    b) bset=1;BDS=$OPTARG;;
    i) iset=1;DIST=$OPTARG;;
    s) sset=1;SCALE=$OPTARG;;
  esac
done
for item in $iset $bset $sset; do
	shift;shift;
done

while [ "$#" != "0" ];do

WDIR=`dirname $1`
ITEM=`basename $1`

XYZFILE1=$WDIR/$ITEM-east.xyz
XYZFILE2=$WDIR/$ITEM-north.xyz
XYZFILE3=$WDIR/$ITEM-up.xyz

if [ -e $XYZFILE1 ]; then
	XYZFILE=$XYZFILE1
else
	XYZFILE=$WDIR/$ITEM
fi

# defaults
if [ "$bset" != "1" ]; then
        BDS=`awk 'BEGIN{xmin=1e20;xmax=-1e20;ymin=1e20;ymax=-1e20}{if ($1<xmin){xmin=$1};if ($1>xmax){xmax=$1};if ($2>ymax){ymax=$2};if ($2<ymin){ymin=$2}}END{printf("%f/%f/%f/%f\n",xmin,xmax,ymin,ymax)}' $XYZFILE`
	xset=1
fi
if [ "$iset" != "1" ]; then
        DIST=`grep -v 'NaN' $XYZFILE | head -n 2 | awk 'function max(a,b){if (a<b){return b}else{return a}};BEGIN{c1=0;c2=0;}{if (2==NR){print max(sqrt(($2-c2)^2),sqrt(($1-c1)^2))};c1=$1;c2=$2}'`
fi

if [ "$xset" != "1" ]; then
	echo ${self}: Processing $WDIR/$ITEM with surface
else
	echo ${self}: Processing $WDIR/$ITEM with xyz2grd
fi
echo $self": Using bounds: "$BDS" and increment "$DIST

	
if [ -e $XYZFILE3 ]; then
	# three files with east, north and up

	GRDFILE1=$WDIR/$ITEM-east.grd
	GRDFILE2=$WDIR/$ITEM-north.grd
	GRDFILE3=$WDIR/$ITEM-up.grd

	surface $XYZFILE1 -G$GRDFILE1 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
	surface $XYZFILE2 -G$GRDFILE2 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
	surface $XYZFILE3 -G$GRDFILE3 -I$DIST/$DIST -R$BDS -T0.25 -C0.001

	if [ "$sset" == "1" ]; then
		if [ "$xset" != "1" ]; then
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE1 | \
				surface -G$GRDFILE1 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE2 | \
				surface -G$GRDFILE2 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE3 | \
				surface -G$GRDFILE3 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
		else
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE1 | \
				xyz2grd -G$GRDFILE1 -I$DIST/$DIST -R$BDS
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE2 | \
				xyz2grd -G$GRDFILE2 -I$DIST/$DIST -R$BDS
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE3 | \
				xyz2grd -G$GRDFILE3 -I$DIST/$DIST -R$BDS
		fi
	else
		if [ "$xset" != "1" ]; then
			surface $XYZFILE1 -G$GRDFILE1 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
			surface $XYZFILE2 -G$GRDFILE2 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
			surface $XYZFILE3 -G$GRDFILE3 -I$DIST/$DIST -R$BDS -T0.25 -C0.001
		else
			xyz2grd $XYZFILE1 -G$GRDFILE1 -I$DIST/$DIST -R$BDS
			xyz2grd $XYZFILE2 -G$GRDFILE2 -I$DIST/$DIST -R$BDS
			xyz2grd $XYZFILE3 -G$GRDFILE3 -I$DIST/$DIST -R$BDS
		fi
	fi

	echo $self": Created files "$GRDFILE1", "$GRDFILE2", "$GRDFILE3
else
	# single file
	XYZFILE=$WDIR/$(basename $ITEM .xyz).xyz
	GRDFILE=$WDIR/$(basename $ITEM .xyz).grd
	if [ "$sset" == "1" ]; then
		if [ "$xset" != "1" ]; then
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE | \
				surface -G$GRDFILE -I$DIST/$DIST -R$BDS -T0.25 -C0.001
		else
			awk -v s=$SCALE '{$3=$3*s;print $0}' $XYZFILE | \
				xyz2grd -G$GRDFILE -I$DIST/$DIST -R$BDS
		fi
	else
		if [ "$xset" != "1" ]; then
			surface $XYZFILE -G$GRDFILE -I$DIST/$DIST -R$BDS -T0.25 -C0.001
		else
			xyz2grd $XYZFILE -G$GRDFILE -I$DIST/$DIST -R$BDS
		fi
	fi
	echo $self": Created file "$GRDFILE
fi

shift
if [ "$#" != "0" ];then
	echo ""
fi
done

