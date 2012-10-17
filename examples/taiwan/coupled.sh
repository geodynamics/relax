#!/bin/bash

# this is Rousset et al. (in prep. 2012) preferred model for the 1999 Mw 7.1 
# Chi-Chi earthquake including viscoelastic relaxation in a three-dimensional 
# domain and afterslip. The original model is referred to as
#
#   mixed_H30_T40_W15_EC10_FW20_51017-51018_V030
#
#

FLT=faults/chichi.flt
WDIR=coupled
OBS=gps/gps.dat

if [ ! -e $WDIR ]; then
	echo $(basename $0): creating directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=4 time relax --no-proj-output --no-vtk-output --no-stress-output $* <<EOF | tee $WDIR/in.param
# grid dimension (sx1,sx2,sx3)
512 512 512
# sampling (dx1,dx2,dx3), smoothing (beta, nyquist)
1 1 1 0.2  1
# origin position (x0,y0) and rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
#120.982 23.772 51 1e3
# observation depth (displacement and stress)
0 5
# output directory
$WDIR
# lambda, mu, gamma (gamma = (1 - nu) rho g / mu)
1 1 8.33e-4
# time interval, (positive time step) or (negative skip, scaling)
20 -1 1 
# number of observation planes
1
# no. x1  x2 x3 length width strike dip
    1 30 -60  0    100    50    120  90
# number of observation points
`grep -v "#" $OBS | wc`
# no NAME x1 x2 x3
`grep -v "#" $OBS`
# number of Coulomb planes
0
# number of prestress interfaces
0
# number of linear viscous interfaces
1
# no. x3 gammadot0 (1/tm) cohesion
    1 30       0.5               0 
# number of viscous zone
1
# no. dgammadot0      x1    x2 x3 length width thickness strike dip
   1           5 -101.48 -0.90 15    200    15        40     20  90
# number of nonlinear viscous interfaces
0
# number of fault creep interfaces
1
# no. x3  Vo (a-b)sigma friction cohesion
    1  0  30          1      0.6        0
# number of creeping faults
2
# no.     x1     x2 x3 length width strike dip rake
    1 -25.00 -32.00  0     80    20      5  30   90
    2 -26.51 -14.75 10     80    20      5   5   90
# number of inter-seismic strike-slip segments
0
# number of inter-seismic opening segments
0
# number of events
1
# number of coseismic shear-slip segments (ns)
`grep -v "#" $FLT | wc `
# no.     slip       xs       ys       zs  length   width strike   dip   rake
`grep -v "#" $FLT`
# number of tensile cracks 
0
# number of mogi source
0
# number of surface traction
0
EOF
