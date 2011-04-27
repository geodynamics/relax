#!/bin/sh

time ./relax $* --no-proj-output <<EOF
# grid size (sx1,sx2,sx3)
256 256 256
# sampling size & smoothing (dx1,dx2,dx3,beta,nyquist)
5 5 5 0.25 2
# origin position, rotation, lon lat
0 0 0 
# geographic origin (longitude, latitude and zone)
-120 34 11 1000
# observation depth
0 0
# output directory
./output1
# elastic parameters and gamma 
1 1 0
# integration time (t1) and time steps
0 -1 1.000
# number of observation planes
0
# number of observation points
0
# number of prestress interfaces with depth
0
# number of linear viscous interfaces
2
# no  x3 gammadot0 cohesion
   1 200         1        0
   2 300         1        0
# number of ductile shear zone
2
# no dgammadot0 x1 x2 x3 length width thickness strike dip
   1          1  0  0 100   100   100        50      0  90
   2          1  0  0 100   100   200        50     40  70
# number of nonlinear viscous interfaces
0
# number of fault creep interfaces
0
# no depth gamma0 (a-b)sigma friction cohesion
#   1   15      1       5e-1      0.6        0
# number of afterslip planes
#1
# no  x1 x2 x3 length width strike dip
#   1 -50 -5 15     90    15  -25.4  90
# number interseismic shear disloc
0
# number interseismic tensile cracks
0
# number of coseismic events
1
# number of shear dislocations
1
# no slip   x1 x2 x3 length width strike dip rake
   1    1 -100  0  0    2e2   1e2      0  90   90
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface tractions
0
EOF

