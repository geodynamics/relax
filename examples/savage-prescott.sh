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

count(){
	echo $#
}

WDIR=./savage-prescott

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

event="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19"

export OMP_NUM_THREADS=8

time ../relax --no-proj-output --no-stress-output $* <<EOF | tee $WDIR/in.param
# grid dimension (sx1,sx2,sx3)
4 1024 512
# sampling (dx1,dx2,dx3), smoothing (beta, nyquist)
0.05 0.05 0.05 0.2 2
# origin position (x0,y0) and rotation
0 0 0
# observation depth (displacement and stress)
0 2
# output directory
$WDIR
# lambda, mu, gamma (gamma = (1 - nu) rho g / mu)
1 1 8.33e-4
# time interval, (positive time step) or (negative skip, scaling)
20 0.1
# number of observation planes
0
# number of observation points
0
# number of stress observation segments
0
# number of prestress interfaces
0
# number of linear viscous interfaces
1
# n. depth gammadot0 cohesion
   1     4         1        0
# number of linear weak zones
0
# number of nonlinear viscous interfaces
0
# number of fault creep interfaces
0
# number of inter-seismic strike-slip segments
2
# n. slip/time xs   ys zs length width strike dip rake
   1        -1 -2 -7.0  0      4     3      0  90    0
   2        -1 -2  7.0  0      4     3      0  90    0
# number of inter-seismic tensile segments
0
# number of events
$(count $event)
# number of coseismic strike-slip segments
1
# n. slip   xs ys zs length width strike dip rake
   1    0 -6.4  0  0   12.8     1      0  90    0
# number of coseismic tensile segments
0
# number of coseismic dilatation point sources
0
# number of surface loads
0
# generate a sequence of earthquakes
`for e in $event; do
	echo \# time of next coseismic event
	echo $e
	echo \# number of coseismic strike-slip segments
	echo 1
	echo \# n. slip   xs ys zs length width strike dip rake
	echo     1    1 -3.2  0  0    6.4     1      0  90    0
	echo \# number of coseismic tensile segments
	echo 0
	echo \# number of coseismic dilatation point sources
	echo 0
	echo \# number of surface loads
	echo 0
done`
EOF
