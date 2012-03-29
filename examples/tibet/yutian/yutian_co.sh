#!/bin/sh

# reference:
#
#   Elliott et al., "Extension on the Tibetan plateau: recent normal faulting measured by InSAR and body wave seismology
#   Geophys. J. Int., vol. 183, pp. 503-535, 2010, doi:10.1111/j.1365-246X.2010.04754.x
# 
# the slip model is in local Cartesian coordinates centered on
#
#   E81.392,N35.445, UTM zone=44 north
#

WDIR=./yutian_co

FLT=faults/yutian.dat

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=4 ../../../build/relax --no-proj-output $* <<EOF | tee $WDIR/in.param
# use '#' character to include comments in your input file
# grid size (sx1,sx2,sx3)
512 512 512
# sampling size (in unit of length), smoothing (0-0.5) & nyquist (dx1,dx2,dx3,beta,nq)
0.3 0.3 0.3 0.2 2
# origin position & rotation
0 0 0
# geographic origin (longitude and latitude), UTM zone and real length unit (usually m or km)
# displacements and stress are converted to lon/lat geographic coordinates
# unit corresponds to a scaling from dx1,dx2,dx3 to real unit
# use unit = 1   if dimensions are described in units of m
# use unit = 1e3 if dimensions are described in units of km
#96.629 33.271 47 1e3
# observation depth for displacement and for stress (stress in only exported in GRD)
0 5
# output directory (all output written here)
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (in unit of time), step (negative for automatic) and scaling of computed value
0 -1 1
# number of observation planes
1
# no x1 x2 x3 length width strike dip
   1  0 -3  0      6     6     90  90
# number of observation points
0
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
`wc $FLT`
# no slip x1 x2 x3 length width strike dip rake
`grep -v "#" $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
EOF
