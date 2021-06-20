#!/bin/bash

# reference:
#
#   Mohamed Chlieh, Jean-Philippe Avouac, Vala Hjorleifsdottir, Teh-Ru Alex Song, 
#   Chen Ji, Kerry Sieh, Anthony Sladen, Helene Hebert, Linette Prawirodirdjo, 
#   Yehuda Bock, and John Galetzka, 
#   
#   "Coseismic Slip and Afterslip of the Great Mw 9.15 Sumatra-Andaman Earthquake of 2004", 
#
#   J. Geophys. Res., vol. 117, B03401, 36 pp., 2012, doi:10.1029/2011JB008868
# 
# the slip model is in local Cartesian coordinates centered on
#
#   E95.854, N3.338298 zone 46
#
# the reference location of the fault patches is not specified in the publication
# therefore the slip model location is approximated.

WDIR=./coseismic

FLT=faults/chlieh+07_km.dat

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=4 relax $* <<EOF | tee $WDIR/in.param
# use '#' character to include comments in your input file
# grid size (sx1,sx2,sx3)
512 512 512
# sampling size (in unit of length), smoothing (0-0.5) & nyquist (dx1,dx2,dx3,beta,nq)
6.0 6.0 6.0 0.2 2
# origin position & rotation
0 0 0
# geographic origin (longitude and latitude), UTM zone and real length unit (usually m or km)
# displacements and stress are converted to lon/lat geographic coordinates
# unit corresponds to a scaling from dx1,dx2,dx3 to real unit
# use unit = 1   if dimensions are described in units of m
# use unit = 1e3 if dimensions are described in units of km
#95.854 3.338298 46 1e3
# observation depth for displacement and for stress (stress in only exported in GRD)
0 5
# output directory (all output written here)
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
30 30 8.33e-4
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
`grep -v "#" $FLT | wc`
# nb slip  x1 x2 x3 length width strike dip rake
`grep -v "#" $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
EOF
