#!/bin/bash

# Mw 7.6 2005 Kashmir (Pakistan) earthquake
# afterslip model
#
# plot the afterslip-induced displacement with
#
#   map.sh -b -40/60/-30/70 -p -0.001/0.001/0.00001 -v 0.0005 -e erpatch.sh afterslip/0{00,01,02,03,04,05,06,07}-relax
#
# in Paraview, load the creep history on the two afterslip planes
#
#   creep-0001s-0*
#   creep-0002s-0*
#
# the geometry of the static afterslip planes can be checked before calculation (with --dry-run option) from the files
#
#   aplane-0001.vtp
#   aplane-0002.vtp
#

FLT=faults/avouac+06.flt

# output directory based on script name.
WDIR=$(basename $0 .sh)

if [ ! -e $WDIR ]; then
	echo $(basename $0): creating directory $WDIR
	mkdir $WDIR
fi

SX=256
DX=0.9

LEN=`echo $SX $DX | awk '{print ($1-1)*$2/2}'`

OMP_NUM_THREADS=4 time relax $* <<EOF | tee $WDIR/in.param
# grid dimension (sx1,sx2,sx3)
$SX $SX $SX
# sampling (dx1,dx2,dx3), smoothing (beta, nyquist)
$DX $DX $DX 0.2 0
# origin position (x0,y0) and rotation
0 0 0
# geographic origin (longitude, latitude, UTM zone, unit)
#-147.73 61.04 6 1e3
# observation depth (displacement and stress)
0 5
# output directory
$WDIR
# lambda, mu, gamma (gamma = (1 - nu) rho g / mu)
30e3 30e3 8.33e-4
# time interval, (positive time step) or (negative skip, scaling)
10 -1 0.5
# number of observation planes
0 
# number of observation points
`cat gps/gps_km.dat | awk -v l=$LEN 'function abs(x){return (x>0)?x:-x}{if (abs($3)<l && abs($4)<l){print $0}}' | wc`
# number NAME x1 x2 x3
`cat gps/gps_km.dat | awk -v l=$LEN 'function abs(x){return (x>0)?x:-x}BEGIN{c=1}{if (abs($3)<l && abs($4)<l){$1=c;print $0;c=c+1}}'`
# number of Coulomb planes
0
# number of prestress interfaces
0
# number of linear viscous interfaces
0
# number of nonlinear viscous interfaces
0
# number of fault creep interfaces
1
# nb depth gamma0 (a-b)sigma friction cohesion
   1     0   1e-3          5      0.6        0
# number of afterslip planes
2
# nb        x1        x2       x3 length width strike dip rake
   1 -21.46101  18.09357        0     60    45    -40  29   90
   2  3.837752  48.24346 21.81643     60    45    -40  15   90
# number of inter-seismic strike-slip segments
0
# number of inter-seismic opening segments
0
# number of events
1
# number of coseismic shear-slip segments (ns)
`grep -v "#" $FLT | wc`
# no.     slip       xs       ys       zs  length   width strike   dip   rake
`grep -v "#" $FLT`
# number of tensile cracks 
0
# number of mogi source
0
# number of surface traction
0
EOF
