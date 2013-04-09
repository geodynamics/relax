#!/bin/bash

# linear viscoelastic relaxation following the 2006 Mw 8.3 Kuril earthquake
# including the effect of an elastic slab for the USGS slab geometry.

# the reference point is at 153.0 47.0, UTM zone 56.

WDIR=viscous3d

FLT=faults/steblov+08.flt
GPS=gps/gps_Kuril.dat

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

OMP_NUM_THREADS=2 time relax --no-proj-output --no-txt-output --no-stress-output $* <<EOF | tee $WDIR/in.param
# grid size (sx1,sx2,sx3)
512 512 512
# sampling size, smoothing & nyquist (dx1,dx2,dx3,beta,nq)
4 4 2 0.2 0
# origin position & rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
#153.0 47.0 56 1e3
# observation depth (displacement and stress) (stress in only exported in GRD)
0 0.5
# output directory
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km (9.147e-4)
1 1 8.33e-4 
# integration time (t1)
1.7 -1 1.0
# number of observation planes
0
# number of observation points
`grep -v "#" $GPS | wc -l`
# index NAME x1 x2 x3
`cat $GPS`
# number of Coulomb patches
0
# number of prestress interfaces
0
# number of linear viscous interfaces
3
# no depth gammadot0 cohesion
    1  63.0      10.5718      0.0
    2 220.0      10.5718      0.0
    3 220.0       0.02114     0.0
# number of linear ductile zones
`grep -v "#" faults/kuril.flt | wc`
# nb dgammadot0   x1   x2   x3   length   width   thickness   strike   dip 
`grep -v "#" faults/kuril.flt | awk '{print NR,-10,$2,$3,$4,$5,42,$6,$7,$8+90}'`
# number of powerlaw viscous interfaces
0
# number of friction faults
0
# number of interseismic loading strike-slip and opening
0
0
# number of coseismic events
1
# number of shear dislocations
`grep -v "#" $FLT | wc`
# no slip     x1         x2        x3  length width strike dip rake
`grep -v "#" $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface loads
0
EOF

