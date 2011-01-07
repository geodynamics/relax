#!/bin/sh

# linear viscoelastic relaxation following a strike-slip fault
# followed by a second strike-slip event and a relaxation of
# both stress perturbations.

# output is exported in geographic coordintes in .xyz format

# viscous substrate starts at 3 times the seismogenic depth.
# output at each computation step (Dt = -1)

# to visualize coseismic deformation (requires GRD output, or manual conversion to GRD format):
# map.sh -b -3/3/-3/3 output1/000

# to convert .xyz output in .grd output in post processing, type (requires GMT installed on your machine)
# xyz2grd.sh output1/000

# to visualize a time series of postseismic deformation (requires GRD output): 
# map.sh -b -5/5/-5/5 -p -0.002/0.002/0.0001 -v 0.005 output1/0{01,02,03,04,05,06,07,08,09,10}-relax
# or simply
# map.sh -b -5/5/-5/5 -p -0.002/0.002/0.0001 -v 0.005 output1/0??-relax

# type map.sh for a description of command-line options.

time ../relax <<EOF
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
./output3
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (t1)
2 -1
# number of observation planes
0
# number of observation points
0
# number of prestress interfaces
0
# number of linear viscous interfaces
0
# number of powerlaw viscous interfaces
0
# number of friction faults
1
# no    depth   gamma0 (a-b)sig friction cohesion
   1       1      1e3      1e3      0.6        0
# number of creeping faults
1
# no  x1 x2 x3 length width strike dip
   1 -1  0 1     2    1      0  90
# number of interseismic loading strike-slip
0
# number of interseismic loading opening cracks
0
# number of coseismic events
1
# number of shear dislocations
1
# no slip  x1 x2 x3 length width strike dip rake
   1   -1 -1  0  0     2    1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
EOF
