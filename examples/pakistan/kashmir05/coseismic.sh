#!/bin/bash

# Mw 7.6 2005 Kashmir (Pakistan) earthquake
# coseismic displacements and stress
#
# coseismic slip distribution from:
#
#   Avouac, J. P., F. Ayoub, S. Leprince, A. O. Konca, and D. V. Helmberger, 
#   "The 2005, Mw 7.6 Kashmir earthquake: sub-pixel correlation of ASTER images and seismic waveforms analysis", 
#   Earth And Planetary Science Letters, 2006. 
#
# origin position
#
#   73.584008666 34.273921999, UTM zone 43
#
# echo 73.584008666 34.273921999 | proj +proj=utm +zone=43
# 369652.60       3793435.75
#

FLT=faults/avouac+06.flt

# output directory based on script name.
WDIR=$(basename $0 .sh)

if [ ! -e $WDIR ]; then
	echo $(basename $0): creating directory $WDIR
	mkdir $WDIR
fi

SX=512
DX=0.8

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
`cat gps/gps_km.dat | awk -v l=$LEN 'function abs(x){return (x>0)?x:-x}{if (abs($3)<l && abs($4)<l){print $0}}' | wc`
# number NAME x1 x2 x3
`cat gps/gps_km.dat | awk -v l=$LEN 'function abs(x){return (x>0)?x:-x}BEGIN{c=1}{if (abs($3)<l && abs($4)<l){$1=c;print $0;c=c+1}}'`
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
