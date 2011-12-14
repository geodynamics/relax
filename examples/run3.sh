#!/bin/sh

# comments:
#
# afterslip following left-lateral strike-slip on a vertical fault.
#
# output at each computation step (Dt = -1).
#
# to visualize coseismic deformation (requires GRD output, or manual conversion to GRD format):
#
#   map.sh -b -3/3/-3/3 -e rpatch.sh output3/000
#
# the option -e rpatch.sh adds the trace of the fault to the map. 
# to visualize the displacement due to afterslip only:
#
#   map.sh -b -5/5/-5/5 -s 0.75 -e rpatch.sh output3/001-relax
#
# the -s 0.75 option optimizes the spacing between horizontal vectors.
# to convert .xyz output in .grd output in post processing, type (requires GMT installed on your machine)
# xyz2grd.sh output3/000
#
# to visualize a time series of postseismic deformation (requires GRD output): 
# map.sh -b -5/5/-5/5 -p -0.002/0.002/0.0001 -v 0.005 output3/0{01,02,03,04,05,06,07,08,09,10}-relax
#
# type map.sh for a description of command-line options.
#
# the conversion to geographic is aborted (--no-proj-output) and 
# the stress component are not exported (--no-stress-output)

WDIR=./output3

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time ../relax --no-proj-output --no-stress-output $* <<EOF | tee $WDIR/in.param
# grid size (sx1,sx2,sx3)
256 256 256 
# sampling size, smoothing & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.2 1
# origin position & rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
# displacements and stress are converted to lon/lat geographic coordinates
# unit corresponds to a scaling from dx1,dx2,dx3 to real unit
# use unit = 1   if dimensions are described in units of m
# use unit = 1e3 if dimensions are described in units of km
#-115 37.5 11 1e3
# observation depth (displacement and stress) (stress in only exported in GRD)
0 0
# output directory
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (t1)
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
0
# number of powerlaw viscous interfaces
0
# number of friction interfaces
1
# no    depth   gamma0 (a-b)sig friction cohesion
   1       1      1e3      1e3      0.6        0
# number of creeping faults
1
# define the fault planes where afterslip occurs. it is possible to
# constrain the rake slip to be within 180 degrees of given direction by
# providing a rake between -180 and 180. for values outside this range,
# the constraint of rake is not applied. this afterslip plane starts
# immediately below the coseismic rupture.
# no  x1 x2 x3 length width strike dip rake (slip rake is +-180 degrees from this value)
   1  -1  0  1     2      1      0  90  180
# number of interseismic loading strike-slip
0
# number of interseismic loading opening cracks
0
# number of coseismic events
1
# number of shear dislocations
1
# no slip  x1 x2 x3 length width strike dip rake
   1   +1  -1  0  0      2     1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface loads
0
EOF
