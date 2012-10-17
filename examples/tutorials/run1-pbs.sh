#!/bin/bash

# example input file for a queueing system (for a single node of a cluster)

#PBS -l nodes=1:ppn=8
#PBS -l walltime=0:02:60
#PBS -q debug
#PBS -N output1
#PBS -o ./output1/out

# change the working directory (default is home directory)
# echo working directory: $PBS_O_WORKDIR
cd $PBS_O_WORKDIR

export MX_RCACHE=0

WDIR=./output1
if [ ! -e $WDIR ]; then
        echo adding directory $WDIR
	mkdir $WDIR
fi

# run this example on pangu, or another PBS queued cluster, with
#
#   qsub ./run1-pbs.sh
#

OMP_NUM_THREADS=12 relax <<EOF
# use '#' character to include comments in your input file
# grid size (sx1,sx2,sx3)
256 256 256
# sampling size (in unit of length), smoothing (0-0.5) & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.2 2
# origin position & rotation
0 0 0
# geographic origin (longitude and latitude), UTM zone and real length unit (usually m or km)
# displacements and stress are converted to lon/lat geographic coordinates
# unit corresponds to a scaling from dx1,dx2,dx3 to real unit
# use unit = 1   if dimensions are described in units of m
# use unit = 1e3 if dimensions are described in units of km
#120 22 51 1e3
# observation depth for displacement and for stress (stress in only exported in GRD)
0 0.5
# output directory (all output written here)
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time (in unit of time), step (negative for automatic) and scaling of computed value
20 -1 1
# number of observation planes
0
# number of observation points
0
# number of Coulomb patches
0
# number of prestress interfaces
0
# number of linear viscous interfaces (where viscosity changes)
2
# no depth gammadot0 cohesion (gammadot0 is shear modulus divided by viscosity)
   1   3.0       1.0      0.0
   2   9.0       1.0      0.0
# number of linear ductile zones
0
# number of powerlaw viscous interfaces
0
# number of friction faults
0
# number of interseismic loading strike-slip and opening
0
0
# number of coseismic events (when slip distribution is prescribed)
2
# number of shear dislocations (strike-slip and dip-slip faulting)
1
# no slip x1 x2 x3 length width strike dip rake
   1    1 -1  0  0      1     1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
# time of second event
1
# number of shear dislocations
1
# no slip x1 x2 x3 length width strike dip rake
   1    1  0  0  0      1     1      0  90    0
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface traction
0
EOF
