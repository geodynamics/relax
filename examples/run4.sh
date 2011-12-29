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

time ../relax --no-stress-output --no-proj-output $* < run4.input
