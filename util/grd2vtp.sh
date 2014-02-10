#!/bin/bash

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

grd2vtp_header(){
	SX1=`grdinfo -C $GRDFILE | awk '{print $10}'`
	SX2=`grdinfo -C $GRDFILE | awk '{print $11}'`
	N=`echo $SX1 $SX2 | awk '{print $1*$2}'`
	TITLE=`grdinfo -C $GRDFILE | awk '{print $1}'`
	NPOLY=`echo $SX1 $SX2 | awk '{print ($1-1)*($2-1)*2}'`
	NSIDE=`echo $SX1 $SX2 | awk '{print ($1-1)*($2-1)*2*4}'`

	echo "# vtk DataFile Version 3.0"
	echo $TITLE
	echo "ASCII"
	echo "DATASET POLYDATA"
	echo "POINTS $N float"
	if [ "$tset" != "1" ]; then
		grd2xyz $GRDFILE | awk -v d=$depth '{print $2,$1,d}'
	else
		grd2xyz $GRDFILE | awk -v d=$depth '{print $2,$1,d-$3}'
	fi
	echo ""
	echo "POLYGONS $NPOLY $NSIDE"
	echo "" | awk -v SX1=$SX1 -v SX2=$SX2 '{ \
		for (x2=1;x2<SX2;x2++){ \
			for (x1=1;x1<SX1;x1++){ \
				print 3,(x2-1)*SX1+x1-1,(x2-1)*SX1+x1,x2*SX1+x1; \
				print 3,(x2-1)*SX1+x1-1,x2*SX1+x1-1,x2*SX1+x1 \
			} \
		} \
	}'
	echo ""
	echo "POINT_DATA $N"
}

grd2vtp_footer(){
	TITLE=`grdinfo -C $GRDFILE | awk '{print $1}'`
	echo "SCALARS $(basename $TITLE .grd) float"
	echo "LOOKUP_TABLE default"
	grd2xyz $GRDFILE | awk '{print $3}'
}



usage(){
	echo "$self converts a list of grd files to a vtp polygon polydata file"
	echo "in ASCII legacy format for visualization in Paraview."
	echo ""
	echo "usage: $self file.grd > file.vtp"
	echo ""
	echo "or"
	echo ""
	echo "       $self file1.grd file2.grd file3.grd > file.vtp"
	echo ""
	echo "options:"
	echo "        -d depth provide the uniform depth of the VTK file (0)."
	echo "        -t indicates that z-value is topography."
	echo "           For multiple files, the first file contains the topography."
	echo "           In combination with -d option, topography is shifted by depth."
	echo ""
	echo "note: in Paraview, the X axis is north, the Y axis is east,"
	echo "      and the Z axis is depth, forming a right-handed system."
	echo ""
}

# print usage is no arguments are supplied or if redirection pipe is detected
if [ ! -t 0 ] || [ $# -eq 0 ]; then
	usage
	exit
fi

# optional command-line parameters
while getopts "d:t" flag
do
  case "$flag" in
    r) dset=1;depth=$OPTARG;;
    t) tset=1;;
  esac
done
for item in $tset $Kset $Oset ;do
	shift
done
for item in $dset ;do
	shift;shift
done

if [ "$dset" != "1" ]; then
	depth=0
fi

GRDFILE=$1
grd2vtp_header
grd2vtp_footer
shift

# loop over list of files to convert
while [ $# -ne 0 ];do
	GRDFILE=$1
	if [ ! -e $GRDFILE ]; then
		echo $self: could not find $GRDFILE. exiting. >&2
		exit 2
	fi
	grd2vtp_footer
	shift
done

