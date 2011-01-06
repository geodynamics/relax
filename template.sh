#!/bin/sh

time ./relax <<EOF
# grid size (sx1,sx2,sx3)
512 512 512
# sampling size & smoothing (dx1,dx2,dx3,beta)
1.0 1.0 1.0 0.2
# origin position and rotation
0 0 0
# observation depth
0
# output directory
./output
# elastic parameters (lambda,mu)
3e1 3e1
# integration time (t1)
1000 10
# number of observation points
12
# index name x1 x2 x3
      1 GPS1 25 10  0
      2 GPS2 25 20  0
      3 GPS3 25 30  0
      4 GPS4 25 40  0
      5 GPS5 25 50  0
      6 GPS6 25 60  0
      7 GPS7 50 10  0
      8 GPS8 50 20  0
      9 GPS9 50 30  0
     10 GP10 50 40  0
     11 GP11 50 50  0
     12 GP12 50 60  0
# number of layers
2
# index depth lambda mu gammadot0 
      1     0      1  1       0.0
      1    20      1  1       1.0
# number of shear dislocations
4
# index  slip x1  x2  x3 length width strike dip rake
      1     2  0 -40   0     80     5     90  90    0
      2     1  0 -40   5     80     5     90  90    0
      3   0.5  0 -40  10     80     5     90  90    0
      4  0.01  0 -40  15     80     5     90  90    0
# number of tensile cracks
1
# index opening x1  x2 x3 length width strike dip
      1      -1  0 -40  5     80    40     90  40
EOF


