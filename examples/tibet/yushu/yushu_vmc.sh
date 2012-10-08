#!/bin/bash

# reference:
#
#   Li et al., "The 2010 MW 6.8 Yushu (Qinghai, China) earthquake: 
#               constraints provided by InSAR and body wave seismology"
#   J. Geophys. Res., vol. 116, B10302, 16 pp., 2011, doi:10.1029/2011JB008358
# 
# the slip model is in local Cartesian coordinates centered on
#
#   E96.629,N33.271, UTM zone=47 north
#
# forward model of viscoelastic flow in the middle crust 
# 

WDIR=./yushu_vmc

FLT=faults/yushu.dat

DDIR=gps

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

SX=512
DX=0.3
L=`echo $SX $DX | awk '{print $1*$2/2}'`

OMP_NUM_THREADS=4 relax --no-proj-output $* <<EOF | tee $WDIR/in.param
# use '#' character to include comments in your input file
# grid size (sx1,sx2,sx3)
$SX $SX $SX
# sampling size (in unit of length), smoothing (0-0.5) & nyquist (dx1,dx2,dx3,beta,nq)
$DX $DX $DX 0.2 2
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
1 -1 1
# number of observation planes
1
# no x1 x2 x3 length width strike dip
   1  0 -3  0      6     6     90  90
# number of observation points
`grep -v "#" $DDIR/gps_km.dat | awk -v l=$L 'function abs(x){return (x<0)?-x:x}{if (abs($3)<l && abs($4)<l){print $0}}' | wc`
# no name x1 x2 x3
`grep -v "#" $DDIR/gps_km.dat | awk -v l=$L 'function abs(x){return (x<0)?-x:x}BEGIN{i=0}{if (abs($3)<l && abs($4)<l){i=i+1;$1=i;print $0}}'`
# number of Coulomb patches
0
# number of prestress interfaces
0
# number of linear viscous interfaces (where viscosity changes)
1
# nb depth gammadot0 cohesion
   1    13        10        0
# number of linear ductile zone
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
