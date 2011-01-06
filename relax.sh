#!/bin/sh

time ./relax <<EOF
# grid size (sx1,sx2,sx3)
256 256 256
# sampling size & smoothing (dx1,dx2,dx3,beta)
0.05 0.05 0.05 0.2
# origin position
0 0
# observation depth
0
# output directory
./output
# elastic parameters (lambda,mu)
1 1
# integration time (t1)
20 0.5
# number of observation planes
0
# number of observation points
0
# number of prestress interfaces
0
# number of linear viscous interfaces
2
1 1.0 0 0.0
2 9.0 0 0.0
# number of powerlaw viscous interfaces
2
1 1.0 1e1 3.0 0.0
2 9.0 1e1 3.0 0.0
# number of friction faults
0
# number of interseismic loading stuff
0
0
# number of coseismic events
2
# number of shear dislocations
1
# index slip   x1 x2 x3 length width strike dip rake
      1    1 -1.0  0  0      1   0.8      0  90    0
# number of tensile cracks
0
# time of second event
10
# number of shear dislocations
1
# index slip  x1 x2 x3 length width strike dip rake
      1 0.02 0.0  0  0      1   0.8      0  90    0
# number of tensile cracks
0
EOF
