#!/bin/sh

# nonlinear viscoelastic relaxation (power exponent n=3)
# following slip on a strike-slip fault
#
# output every computational step (Dt = -2). 
# The time step is computed automatically and the value unchanged (scale = 1)
#
# to visualize coseismic deformation (requires GRD output, or manual conversion to GRD format):
# map.sh -b -3/3/-3/3 output2/000
#
# to visualize postseismic deformation (requires GRD output): 
# map.sh -b -5/5/-5/5 output2/0{02,04,06,08}-relax
#
# type map.sh for a description of command-line options.
#
# the output projected in geographical coordinates is cancelled (--no-proj-output)

WDIR=./output2

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time ../relax --no-proj-output $* <<EOF | tee output2/in.param
# grid size (sx1,sx2,sx3)
256 256 256
# sampling size, smoothing & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.2 2
# origin position & rotation
0 0 0
# observation depth (displacements and stress) (stress is only exported in GRD)
0 0.5
# output directory
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (t1)
50 -2 1
# number of observation planes
0
# number of observation points
0
# number of prestress interfaces
0
# number of linear viscous interfaces
0
# number of powerlaw viscous interfaces
2
# no depth gammadot0 power cohesion
   1   3.0       5e3   3.0      0.0
   2   8.0       5e3   3.0      0.0
# number of power-law ductile zones
0
# number of friction faults
0
# number of interseismic loading strike-slip and opening
0
0
# number of coseismic events
1
# number of shear dislocations
1
# index slip x1 x2 x3 length width strike dip rake
      1    1  1  0  0      2     1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface loads
0
EOF
