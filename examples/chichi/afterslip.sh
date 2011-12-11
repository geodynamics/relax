#/bin/sh

# model of afterslip on down-dip extension of the Chelungpu Fault.
# notice the reduced sampling size to better resolve deformation
# localization on faults.

FLT=faults/chichi.flt
WDIR=$(basename $0 .sh)

if [ ! -e $WDIR ]; then
	echo $(basename $0): creating directory $WDIR
	mkdir $WDIR
fi

time ../../relax --no-proj-output --no-stress-output $* <<EOF | tee $WDIR/in.param
# grid dimension (sx1,sx2,sx3)
512 512 512
# sampling (dx1,dx2,dx3), smoothing (beta, nyquist)
0.35 0.35 0.35 0.2 1
# origin position (x0,y0) and rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
#120.982 23.772 51 1e3
# observation depth (displacement and stress)
0 5
# output directory
$WDIR
# lambda, mu, gamma (gamma = (1 - nu) rho g / mu)
1 1 1e-5
# time interval, (positive time step) or (negative skip, scaling)
1 -1 1 
# number of observation planes
0 
# number of observation points
0 
# number of prestress interfaces
0
# number of linear viscous interfaces
0
# number of nonlinear viscous interfaces
0
# number of fault creep interfaces
1
# no.    depth   gamma0 (a-b)sig friction cohesion
    1        1        5      1e1      0.6        0
# number of afterslip planes
1
# no.        x1       x2       x3   length    width strike    dip   rake
    1     -6.85   -14.62     5.24       60       20     10      5     90
# number of inter-seismic strike-slip segments
0
# number of inter-seismic opening segments
0
# number of events
1
# number of coseismic shear-slip segments (ns)
`grep -v "#" $FLT | wc | awk '{print $1}'`
# no.     slip       xs       ys       zs  length   width strike   dip   rake
`grep -v "#" $FLT`
# number of tensile cracks 
0
# number of mogi source
0
# number of surface traction
0
EOF
