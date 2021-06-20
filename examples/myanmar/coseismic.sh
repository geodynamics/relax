#!/bin/bash

# Mw 6.9 2011 Tarlay (Myanmar) earthquake
# coseismic displacements and stress
#
# origin position
#
#   99.949, 20.705, UTM zone 47
#
# echo 99.949 20.705 | proj +proj=utm +zone=47
# 598824.34       2289789.12
#

FLT=faults/wang+13.flt

# output directory based on script name.
WDIR=$(basename $0 .sh)

if [ ! -e $WDIR ]; then
	echo $(basename $0): creating directory $WDIR
	mkdir $WDIR
fi

SX=256
DX=0.5

LEN=`echo $SX $DX | awk '{print ($1-1)*$2/2}'`

OMP_NUM_THREADS=2 time relax $* <<EOF | tee $WDIR/in.param
# grid dimension (sx1,sx2,sx3)
$SX $SX $SX
# sampling (dx1,dx2,dx3), smoothing (beta, nyquist)
$DX $DX $DX 0.2 0
# origin position (x0,y0) and rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
#73.584008666 34.273921999 6 1e3
# observation depth (displacement and stress)
0 5
# output directory
$WDIR
# lambda, mu, gamma (gamma = (1 - nu) rho g / mu)
30e3 30e3 8.33e-4
# time interval, (positive time step) or (negative skip, scaling)
0 -1 1 
# number of observation planes
0 
# number of observation points
0
# number of Coulomb planes
0
# number of prestress interfaces
0
# number of linear viscous interfaces
0
# number of nonlinear viscous interfaces
0
# number of fault creep interfaces
0
# number of inter-seismic strike-slip segments
0
# number of inter-seismic opening segments
0
# number of events
1
# number of coseismic shear-slip segments (ns)
`grep -v "#" $FLT | wc`
# no.     slip       xs       ys       zs  length   width strike   dip   rake
`grep -v "#" $FLT`
# number of tensile cracks 
0
# number of mogi source
0
# number of surface traction
0
EOF
