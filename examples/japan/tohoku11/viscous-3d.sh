#!/bin/bash

# postseismic relaxation from the April 2012 Offshore Sumatra earthquake
#
# reference point: 142.369E, 38.322N, UTM zone 54 (619669.75,4242429.17)
#

FLT=faults/wei+11.flt

# output directory based on script name.
WDIR=./viscous-3d

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

SX=256
DX=10

LEN=`echo $SX $DX | awk '{print ($1-1)*$2/2}'`

OMP_NUM_THREADS=2 time relax --no-proj-output --no-stress-output $* <<EOF | tee $WDIR/in.param
# grid size (sx1,sx2,sx3)
$SX $SX $SX
# sampling size, smoothing & nyquist (dx1,dx2,dx3,beta,nq)
$DX $DX $DX 0.2 0
# origin position & rotation
0 0 0
# observation depth (displacements and stress) (stress is only exported in GRD)
0 10
# output directory (all output files are stored in the following directory)
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (t1)
10 -1 0.5
# number of observation planes
0
# number of observation points
0
# number of stress observation segments
0
# number of prestress interfaces
0
# number of linear viscous interfaces
1
# no depth gammadot0 cohesion
   1  80.0       1e0      0.0
# number of linear ductile zones
`grep -h -v "#" faults/{izu,kurils,ryukyu}_30km.flt | wc`
# nb dgammadot0   x1   x2   x3   length   width   thickness   strike   dip 
`grep -h -v "#" faults/{izu,kurils,ryukyu}_30km.flt | awk '{print NR,-10,$2,$3,$4,$5,80,$6,$7,$8+90}'`
# number of powerlaw viscous interfaces
0
# number of friction interfaces
0
# number of interseismic loading strike-slip and opening
0
0
# number of coseismic events
1
# number of shear dislocations (strike-slip and dip-slip faults)
`grep -v "#" $FLT | wc`
# index slip x1 x2 x3 length width strike dip rake
`grep -v "#" $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface loads
0
EOF
