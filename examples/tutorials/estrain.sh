#!/bin/bash

# displacement and stress kernels for distributed strain

WDIR=./estrain-e33-dip-60

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=12 relax --no-proj-output --with-eigenstrain $* <<EOF | tee $WDIR/in.param
# use '#' character to include comments in your input file
# grid size (sx1,sx2,sx3)
512 512 512
# sampling size (in unit of length), smoothing (0-0.5) & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.1 2
# origin position & rotation
0 0 0
# observation depth for displacement and for stress (stress in only exported in GRD)
0 2
# output directory (all output written here)
$WDIR
# lambda, mu, gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 0
# integration time (in unit of time), step (negative for automatic) and scaling of computed value
0 -1 1
# number of observation planes
2
# n x1 x2 x3 length width strike dip
  1  0 -6  0     12     6     90  90
  2 -6  0  0     12     6      0  90
# number of observation points
2
# n name x1 x2 x3
  1 GPS1  0  1  0
  2 KER1  0  0  2
# number of Coulomb patches
0
# number of prestress interfaces
0
# number of linear viscous interfaces (where viscosity changes)
0
# number of powerlaw viscous interfaces
0
# number of friction faults
0
# number of interseismic loading strike-slip and opening
0
0
# number of coseismic events (when slip distribution is prescribed)
1
# number of shear dislocations (strike-slip and dip-slip faulting)
0
# number of tensile cracks
0
# number of dilatation sources
0
# number of distributed eigenstrain
1
# n     e11  e12  e13     e22  e23     e33 x1 x2 x3 length width thickness strike dip
  1 0.00e-3 0e-3 0e-3 1.00e-3 0e-3 0.00e-3 -1  1  1      2     2         2      0  90
# number of surface traction
0
EOF

