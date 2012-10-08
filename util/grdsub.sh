#!/bin/bash

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if (test $# -lt "3"); then
        echo "usage: $self item1 item2 item3"
	echo "       computes item3-{north,east,up} = item1-{north,east,up} - item2-{north,east,up}."
	echo "       use $self item1 0 item3 to copy three grd files to item3-{north,east,up}"
	exit 1
fi

ITEMA=$1
ITEMB=$2
ITEMC=$3
echo $self": Using directory "$(dirname $ITEMA)

GRDFILEA1=$ITEMA-east.grd
GRDFILEA2=$ITEMA-north.grd
GRDFILEA3=$ITEMA-up.grd

GRDFILEB1=$ITEMB-east.grd
GRDFILEB2=$ITEMB-north.grd
GRDFILEB3=$ITEMB-up.grd

GRDFILEC1=$ITEMC-east.grd
GRDFILEC2=$ITEMC-north.grd
GRDFILEC3=$ITEMC-up.grd

if [ -e $ITEMA ]; then
	if [ "$ITEMB" == "0" ]; then
		echo ${self}: copying to $(dirname $ITEMC), file $(basename $ITEMC)
		cp $ITEMA $ITEMC
	else
		echo ${self}: subtracting to file $(basename $ITEMC)
		grdmath $ITEMA $ITEMB SUB = $ITEMC
	fi
else
	if [ "$ITEMB" == "0" ]; then
		echo $self": copying to "$(dirname $ITEMC)", files "$(basename $GRDFILEC1), $(basename $GRDFILEC2), $(basename $GRDFILEC3)
		cp $GRDFILEA1 $GRDFILEC1
		cp $GRDFILEA2 $GRDFILEC2
		cp $GRDFILEA3 $GRDFILEC3
	else
		echo $self": subtracting to files $(basename $GRDFILEC1), $(basename $GRDFILEC2), $(basename $GRDFILEC3)"

		grdmath $GRDFILEA1 $GRDFILEB1 SUB = $GRDFILEC1
		grdmath $GRDFILEA2 $GRDFILEB2 SUB = $GRDFILEC2
		grdmath $GRDFILEA3 $GRDFILEB3 SUB = $GRDFILEC3
	fi
fi




