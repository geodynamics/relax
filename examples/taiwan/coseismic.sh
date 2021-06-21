#!/bin/bash

# 1999 Mw 7.6 Chi Chi coseismic displacements
# geographic origin (longitude, latitude, UTM zone)
# 120.982 23.772 51

FLT=faults/chichi.flt
# output directory based on script name.
WDIR=$(basename $0 .sh)

if [ ! -e $WDIR ]; then
	echo $(basename $0): creating directory $WDIR
	mkdir $WDIR
fi

time relax $* <<EOF | tee $WDIR/in.param
# grid dimension (sx1,sx2,sx3)
512 512 512
# sampling (dx1,dx2,dx3), smoothing (beta, nyquist)
0.7 0.7 0.7 0.2 1
# origin position (x0,y0) and rotation
0 0 0
# observation depth (displacement and stress)
0 5
# output directory
$WDIR
# lambda, mu, gamma (gamma = (1 - nu) rho g / mu)
1 1 1e-5
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
$(grep -cv "#" "$FLT")
# no.     slip       xs       ys       zs  length   width strike   dip   rake
$(grep -v "#" "$FLT")
# number of tensile cracks 
0
# number of mogi source
0
# number of surface traction
0
EOF
