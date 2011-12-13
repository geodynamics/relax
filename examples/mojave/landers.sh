#!/bin/sh

# Landers earthquake 
# Newtonian viscosity
# weak lower crust   H = 25-30 km
# weak asthenosphere H = 60+ km
# Origin of coordinates system: -116.27E, 34.595N and 566940.91E, 3828373.73N in UTM zone 11.

WDIR=landers
FLT=faults/landers_km.flt

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time ../../relax <<EOF $* --no-vtk-output --no-proj-output \
	--no-stress-output | tee $WDIR/in.param
# grid size (sx1,sx2,sx3)
512 512 512
# dx1  dx2  dx3 beta nyquist
  0.6  0.6  0.6 0.30     2.0
# origin position and rotation
0 0 0
# geographic origin zone unit
#-116.27 34.595 11 1e3
# observation depth
0 5
# output directory
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (t1) and time steps
20 -1 1
# number of observation planes
0
# number of observation points
0
# number of stress observation planes
0
# number of prestress interfaces with depth
0
# number of linear viscous interfaces
0
# number of powerlaw viscous interfaces
5
# no depth gammadot0 power cohesion
   1    25       0.3   1.0      0.0
   2    30       0.3   1.0      0.0
   3    30         0   1.0      0.0
   4    60         0   1.0      0.0
   5    60         3   1.0      0.0
# number of nonlinear ductile zones
0
# number of fault creep interface
0
# number of interseismic dislocation
0
# number of interseismic dykes
0
# number of coseismic events
1
# number of shear dislocations
`awk 'BEGIN{c=0}{if ($5 > 2 && $5 < 20){c=c+1}}END{print c}' $FLT`
# index slip x1 x2 x3 length width strike dip rake
`awk 'BEGIN{c=1}{if ($5 > 2 && $5 < 20){$1=c;print $0;c=c+1}}' $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
EOF

