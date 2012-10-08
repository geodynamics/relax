#!/bin/bash

# reference for the slip model:
#
#   Xiaopeng Tong, David T. Sandwell, and Yuri Fialko,
#   "Coseismic slip model of the 2008 Wenchuan earthquake derived from joint inversion of interferometric synthetic aperture radar, GPS, and field data"
#   JGR, Vol 115, B04314, 19 PP., 2010, doi:10.1029/2009JB006625
#
# the slip model is in local Cartesian coordinates centered on
#
#   E104.2,N31.4, UTM zone=48 north
#

WDIR=./wenchuan

#FLT=faults/tong+10.flt
FLT=faults/shen+09.flt

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

SX=512
DX=2

R=`echo $SX $DX | awk '{print ($1-1)*$2/2}'`

OMP_NUM_THREADS=2 relax --no-proj-output $* <<EOF | tee $WDIR/in.param
# use '#' character to include comments in your input file
# grid size (sx1,sx2,sx3)
$SX $SX $SX
# sampling size (in unit of length), smoothing (0-0.5) & nyquist (dx1,dx2,dx3,beta,nq)
$DX $DX $DX 0.2 0
# origin position & rotation
0 0 0
# geographic origin (longitude and latitude), UTM zone and real length unit (usually m or km)
# displacements and stress are converted to lon/lat geographic coordinates
# unit corresponds to a scaling from dx1,dx2,dx3 to real unit
# use unit = 1   if dimensions are described in units of m
# use unit = 1e3 if dimensions are described in units of km
#96.500 36.0 47 1e3
# observation depth for displacement and for stress (stress in only exported in GRD)
0 10
# output directory (all output written here)
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
30 30 8.33e-4
# integration time (in unit of time), step (negative for automatic) and scaling of computed value
1 -1 1
# number of observation planes
0
# no   x1 x2 x3 length width strike dip
# number of observation points
1
# i NAME x1 x2 x3
  1 GPS1  0  0  0
# number of Coulomb patches
0
# number of prestress interfaces
0
# number of linear viscous interfaces (where viscosity changes)
1
# nb depth gammadot0 cohesion
   1    20         1        0
# number of ductile zones
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
`grep -v "#" $FLT | wc `
# no slip x1 x2 x3 length width strike dip rake
`grep -v "#" $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
EOF
