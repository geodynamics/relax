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
# modify the number of threads used with all parallel programs in the session with
#
#   export OMP_NUM_THREADS=N with sh and ksh shells
#
# or
#
#   setenv OMP_NUM_THREADS N with csh
#
# alternatively, alter the number of threads used for the current command line with
#
#   OMP_NUM_THREADS=N command
#
# where N is the desired number of threads.

WDIR=./output1

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=4 ../../build/relax $* < run1.input | tee $WDIR/in.param
