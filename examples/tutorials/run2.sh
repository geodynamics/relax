#!/bin/sh

# nonlinear viscoelastic relaxation (power exponent n=3)
# following slip on a strike-slip fault
#
# output every two computational steps (Dt = -2). 
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

OMP_NUM_THREADS=2 time ../../build/relax --no-proj-output $* < run2.input | tee $WDIR/in.param
