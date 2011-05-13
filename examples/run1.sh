#!/bin/sh

#
# comments:
#
# linear viscoelastic relaxation following a strike-slip fault
# followed by a second strike-slip event and a relaxation of
# both stress perturbations.
#
# run this example with:
#
#   ./run1.sh
#
# you can also abort computation to produce only geometry information with:
#
#   ./run1.sh --dry-run
#
# viscous substrate starts at 3 times the seismogenic depth.
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
#   map.sh -b -5/5/-5/5 -p -0.002/0.002/0.0001 -v 0.005 output1/0{01,02,03,04,05,06,07,08,09,10}-relax
#
# type map.sh for a description of command-line options. 
# the command used to generate a map can be retrieve from the .ps file with
#
#   tail -n 1 output1/000-plot.ps
#
# to visualize in 3-D with Paraview (www.paraview.org),
# load the following files:
# 
#   output1/cgrid.vtp           for the computational grid spatial extent
#   output1/rfaults-001.vtp     for the slip distribution of the first event
#   output1/rfaults-002.vtp     for the slip distribution of the second event
#   output1/disp-000.vtr        for the coseismic displacement field
#   output1/linearlayer-001.vtp for the depth of brittle-ductile transition
#   output1/vel-001.vtr         for the instantaneous velocity field
#
# note that we save the entire output to output1/in.param. this is convenient
# to document how the simulation was computed.
#
# type relax --help to display options.
#
# modify the number of threads used with
#
#   export OMP_NUM_THREADS=N with sh and ksh shells
#
# or
#
#  setenv OMP_NUM_THREADS N with csh
#

WDIR=./output1

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time ../relax $* <<EOF | tee output1/in.param
# use '#' character to include comments in your input file
# grid size (sx1,sx2,sx3)
256 256 256
# sampling size (in unit of length), smoothing (0-0.5) & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.2 2
# origin position & rotation
0 0 0
# geographic origin (longitude and latitude), UTM zone and real length unit (usually m or km)
# displacements and stress are converted to lon/lat geographic coordinates
# unit corresponds to a scaling from dx1,dx2,dx3 to real unit
# use unit = 1   if dimensions are described in units of m
# use unit = 1e3 if dimensions are described in units of km
120 22 51 1e3
# observation depth for displacement and for stress (stress in only exported in GRD)
0 0.5
# output directory (all output written here)
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (in unit of time), step (negative for automatic) and scaling of computed value
20 -1 1
# number of observation planes
0
# number of observation points
0
# number of prestress interfaces
0
# number of linear viscous interfaces (where viscosity changes)
2
# no depth gammadot0 cohesion (gammadot0 is shear modulus divided by viscosity)
   1   3.0       1.0      0.0
   2   9.0       1.0      0.0
# number of linear ductile zones
0
# number of powerlaw viscous interfaces
0
# number of friction faults
0
# number of interseismic loading strike-slip and opening
0
0
# number of coseismic events (when slip distribution is prescribed)
2
# number of shear dislocations (strike-slip and dip-slip faulting)
1
# no slip x1 x2 x3 length width strike dip rake
   1    1 -1  0  0      1     1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
# time of second event
1
# number of shear dislocations
1
# no slip x1 x2 x3 length width strike dip rake
   1    1  0  0  0      1     1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
EOF
