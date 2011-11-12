#!/bin/sh

# nonlinear viscoelastic relaxation (power exponent n=3)
# following slip on a strike-slip fault
#
# output every two computational steps (Dt = -2). 
# The time step is computed automatically and the value unchanged (scale = 1)
#
# to visualize coseismic deformation (requires GRD output, or manual conversion to GRD format):
# map.sh -b -3/3/-3/3 output2/000
#
# to visualize postseismic deformation (requires GRD output): 
# map.sh -b -5/5/-5/5 output2/0{02,04,06,08}-relax
#
# type map.sh for a description of command-line options.
#
# the output projected in geographical coordinates is cancelled (--no-proj-output)

WDIR=./output2

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time ../relax --no-proj-output $* <<EOF | tee $WDIR/in.param
# grid size (sx1,sx2,sx3)
256 256 256
# sampling size defines the grid spacing in units of distance
# smoothing is a roll-off parameter 0 <= beta <= 0.5 that tapers 
# the slip distribution on faults beta=0.2 and beta=0.3 are good values
# nyquist is the minimum size of fault patches, relative to the sampling size
# sampling size, smoothing & nyquist (dx1,dx2,dx3,beta,nq)
0.05 0.05 0.05 0.2 2
# origin position & rotation
0 0 0
# observation depth (displacements and stress) (stress is only exported in GRD)
0 0.5
# output directory (all output files are stored in the following directory)
$WDIR
# the elastic parameters (Lambda and mu) are the Lame parameter and the shear modulus
# gamma=(1-nu)*rho*g/mu is a wavelength relevant to buoyancy.
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
1 1 8.33e-4
# integration time refers to the length in units of time of the simulation.
# integration time (t1)
100 -1 1
# observation planes are arbitrary planes where the inelastic strain rate
# is sampled for output.
# number of observation planes
0
# observation points are GPS points, for example, where time series of
# displacement are output.
# number of observation points
0
# stress observation segments are patches where Coulomb and other
# stress components are evaluated.
# number of stress observation segments
0
# prestress interfaces corresponds to depths where the prestress
# tensor changes (default is no prestress).
# number of prestress interfaces
0
# viscous interfaces are depths where the viscous properties change.
# number of linear viscous interfaces
0
# number of powerlaw viscous interfaces
2
# no depth gammadot0 power cohesion
   1   2.0       5e3   3.0      0.0
   2   8.0       5e3   3.0      0.0
# ductile zones corresponds to volumes where viscous properties change.
# number of power-law ductile zones
0
# depths where friction properties change.
# number of friction interfaces
0
# number of interseismic loading strike-slip and opening
0
0
# number of coseismic events
1
# number of shear dislocations (strike-slip and dip-slip faults)
1
# index slip x1 x2 x3 length width strike dip rake
      1    1 -1  0  0      2     1      0  90    0
# number of tensile cracks
0
# volumetric deformation (Mogi source)
# number of dilatation sources
0
# surface loads can be used to model loading of lakes or dams
# number of surface loads
0
EOF
