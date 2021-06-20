#/bin/bash

# Mw 8.6 03/28/2005 Nias earthquake
#
# reference point: 100E, 0N, utm zone 47, 611280.65 E 0 N
#

FLT=faults/shao+ji_km.dat

# output directory based on script name.
WDIR=$(basename $0 .sh)

if [ ! -e $WDIR ]; then
	echo $(basename $0): creating directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=4 time relax $* <<EOF | tee $WDIR/in.param
# grid dimension (sx1,sx2,sx3)
512 512 512
# sampling (dx1,dx2,dx3), smoothing (beta, nyquist)
4 4 4 0.2 1
# origin position (x0,y0) and rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
#100 0 47 1e3
# observation depth (displacement and stress)
0 5
# output directory
$WDIR
# lambda, mu, gamma (gamma = (1 - nu) rho g / mu)
30 30 1e-5
# time interval, (positive time step) or (negative skip, scaling)
10 -1 1 
# number of observation planes
0 
# number of observation points
0 
# number of Coulomb planes
0
# number of prestress interfaces
0
# number of linear viscous interfaces
1
# nb depth gammadot0 cohesion
   1    40         1        0
# number of linear ductile zones
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
