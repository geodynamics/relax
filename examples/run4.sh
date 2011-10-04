#!/bin/sh

#
# comments:
#
# linear viscoelastic relaxation in the lower crust and 
# in a confined damage zone around the deep extension of a vertical fault
# following right-lateral strike-slip on the upper segment of the fault.
#
# output is exported in geographic coordintes in .xyz format
#
# output at each computation step (Dt = -1)
#
# to visualize coseismic deformation (requires GRD output, or manual conversion to GRD format):
#
#   map.sh -b -3/3/-3/3 output1/000
#
# to convert .xyz output in .grd output in post processing, type (requires GMT installed on your machine)
#
#   xyz2grd.sh output1/000
#
# to visualize a time series of postseismic deformation (requires GRD output): 
#
#   map.sh -b -5/5/-5/5 -p -0.002/0.002/0.0001 -v 0.005 output1/00{1,2,3,4}-relax
#
# type map.sh for a description of command-line options.
# type ../relax --help for more info about command-line options.
#
# to visualize in 3-D with Paraview (www.paraview.org),
# load the following files:
# 
#   output4/cgrid.vtp           for the computational grid spatial extent
#   output4/rfaults-001.vtp     for the coseismic slip distribution 
#   output4/disp-000.vtr        for the coseismic displacement field
#   output4/linearlayer-001.vtp for the depth of brittle-ductile transition
#   output4/weakzone-001.vtp    for the spatial extent of the anomalous viscous zone 
#   output4/vel-001.vtr         for the instantaneous velocity field

WDIR=./output4

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time ../relax --no-stress-output $* <<EOF
# grid size (sx1,sx2,sx3)
256 256 256 
# sampling size, smoothing & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.2 1
# origin position & rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
-120 34 11 1e3
# observation depth (displacement and stress) (stress in only exported in GRD)
0 0
# output directory
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (t1), time step and scaling
2 -1 1
# number of observation planes
0
# number of observation points
0
# number of Coulomb patches
0
# number of prestress interfaces
0
# number of linear viscous interfaces
3
# gammadot0 is assumed to be 0 between depths of 0 and 2.
# the last value is extended down to infinity.
# no depth gammadot0 cohesion
   1     2         1        0
   2     3         1        0
   3     3         0        0
# number of ductile zones (this defines volumes where viscosity changes)
1
# a confined area below the ruptured fault is allowed to flow viscously.
# no gammadot0 x1 x2  x3 length width thickness strike dip
   1         1 -1  0 1.1      2   0.9         1      0  90
# number of powerlaw viscous interfaces
0
# number of friction faults
0
# number of interseismic loading strike-slip
0
# number of interseismic loading opening cracks
0
# number of coseismic events
1
# number of shear dislocations
1
# no slip  x1 x2 x3 length width strike dip rake
   1   -1  -1  0  0      2     1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
EOF
