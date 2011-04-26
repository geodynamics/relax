#!/bin/sh

# linear viscoelastic relaxation following a strike-slip fault
# followed by a second strike-slip event and a relaxation of
# both stress perturbations.

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

time ../relax <<EOF | tee output1/in.param
# grid size (sx1,sx2,sx3)
256 256 256
# sampling size, smoothing & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.2 2
# origin position & rotation
0 0 0
# geographic origin, UTM zone and length unit
120 22 51 1e3
# observation depth (displacement and stress) (stress in only exported in GRD)
0 0.5
# output directory
./output1
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (t1), step and scaling
20 -1 1
# number of observation planes
0
# number of observation points
0
# number of prestress interfaces
0
# number of linear viscous interfaces
2
# no depth gammadot0 cohesion
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
# number of coseismic events
2
# number of shear dislocations
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
