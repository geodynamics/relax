#!/bin/sh

# comments:
#
# afterslip following left-lateral strike-slip on a vertical fault.
#
# output at each computation step (Dt = -1) with a reduced time step (1/100 factor) for
# better numerical accuracy.
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
#
#   xyz2grd.sh output3/000
#
# to visualize a time series of postseismic deformation (requires GRD output): 
#
#   map.sh -b -5/5/-5/5 -p -0.001/0.001/0.00001 -v 0.001 -s 1 -e rpatch.sh output3/0{00,01,02,03,04,05,06,07,08,09,10}-relax
#
# to visualize a time series of fault creep:
#
#   map.sh -t 1 -s 0.1 -p -0.2/0.2/0.01 output3/0??.s00001.creep-east.grd
#
# type map.sh for a description of command-line options.
#
# the conversion to geographic is cancelled (--no-proj-output) and 
# the stress component are not exported (--no-stress-output)
#
# TIP:
# the following line in the input file
#
#   # integration time (t1)
#   0.5 -1 0.01
#
# reduces the time step suggested by the code by a factor of 1/100. This choice
# of time step leads to slowly decreasing power density
#                                                                |
#                                                               \./
#    I  |   Dt   | tm(ve) | tm(pl) | tm(as) |     t/tmax     | power  |  C:E^i |
#   000* 0.00E+00 0.00E+00 0.00E+00 0.00E+00 0.00E+00/5.00E-1 0.00E+00 1.03E+00
#   001* 2.00E-02 1.00E+07 1.00E+07 2.00E+01 2.00E-02/5.00E-1 6.47E-01 1.04E+00
#   002* 2.00E-02 1.00E+07 1.00E+07 2.00E+01 4.00E-02/5.00E-1 6.11E-01 1.05E+00
#   003* 2.00E-02 1.00E+07 1.00E+07 2.00E+01 6.00E-02/5.00E-1 5.77E-01 1.06E+00
#   004* 2.00E-02 1.00E+07 1.00E+07 2.00E+01 8.00E-02/5.00E-1 5.47E-01 1.07E+00
#   005* 2.00E-02 1.00E+07 1.00E+07 2.00E+01 1.00E-01/5.00E-1 5.18E-01 1.08E+00
#                                                               /^\
#                                                                |
#
# It is informative to run the same code with a larger time step, for example
# with the following parameters:
#
#   # integration time (t1)
#   5 -1 0.1
#
# which does not lead to numerical convergence.
#                                                                |
#                                                               \./
#    I  |   Dt   | tm(ve) | tm(pl) | tm(as) |     t/tmax     | power  |  C:E^i |
#   000* 0.00E+00 0.00E+00 0.00E+00 0.00E+00 0.00E+00/5.00E+0 0.00E+00 1.03E+00
#   001* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 2.00E-01/5.00E+0 4.69E-01 1.12E+00
#   002* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 4.00E-01/5.00E+0 3.13E-01 1.18E+00
#   003* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 6.00E-01/5.00E+0 2.18E-01 1.22E+00
#   004* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 8.00E-01/5.00E+0 1.55E-01 1.26E+00
#   005* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 1.00E+00/5.00E+0 1.18E-01 1.28E+00
#   006* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 1.20E+00/5.00E+0 9.96E-02 1.30E+00
#   007* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 1.40E+00/5.00E+0 9.47E-02 1.31E+00
#   008* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 1.60E+00/5.00E+0 1.02E-01 1.33E+00
#   009* 2.00E-01 1.00E+07 1.00E+07 2.00E+01 1.80E+00/5.00E+0 1.35E-01 1.35E+00
#                                                               /^\
#                                                                |
# in general, a slowly decreasing power density (about 10-20% per step) is
# a sign of an accurate computation.

WDIR=./output3

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time ../relax --no-proj-output --no-stress-output $* < run3.input | tee $WDIR/in.param
