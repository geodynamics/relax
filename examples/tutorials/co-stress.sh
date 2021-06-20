#!/bin/bash

#
# comments:
#
# coseismic stress perturbation (map view & cross section) due to
# 1 m slip on a 20x20 km thrust fault.
#
# run this example with:
#
#   ./co-stress.sh
#
# you can also abort computation to produce only geometry information with:
#
#   ./co-stress.sh --dry-run
#
# The output units are in MPa
#
# to visualize coseismic deformation (requires GRD output, or manual conversion to GRD format):
#
#   grdmap.sh -b -3/3/-3/3 co-stress/000
#
# to convert .xyz output in .grd output in post processing, type (requires GMT installed on your machine)
#
#   xyz2grd.sh co-stress/000
#
# to visualize a time series of postseismic deformation (requires GRD output): 
#
#   grdmap.sh -b -5/5/-5/5 -p -0.002/0.002/0.0001 -v 0.005 co-stress/0{01,02,03,04,05,06,07,08,09,10}-relax
#
# type map.sh for a description of command-line options. 
# the command used to generate a map can be retrieve from the .ps file with
#
#   tail -n 1 co-stress/000-plot.ps
#
# to visualize in 3-D with Paraview (www.paraview.org),
# load the following files:
# 
#   co-stress/cgrid.vtp           for the computational grid spatial extent
#   co-stress/rfaults-001.vtp     for the slip distribution of the first event
#   co-stress/rfaults-002.vtp     for the slip distribution of the second event
#   co-stress/disp-000.vtr        for the coseismic displacement field
#   co-stress/linearlayer-001.vtp for the depth of brittle-ductile transition
#   co-stress/vel-001.vtr         for the instantaneous velocity field
#
# note that we save the entire output to co-stress/in.param. this is convenient
# to document how the simulation was computed.
#
# type relax --help to display options.
#
# modify the number of threads used with
#
#   export OMP_NUM_THREADS=N with bash, sh and ksh shells
#
# or
#
#  setenv OMP_NUM_THREADS N with csh
#

WDIR=./co-stress

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=4 relax $* < co-stress.input | tee $WDIR/in.param
