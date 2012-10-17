#!/bin/bash

# viscoelastic relaxation following a dilatation (Mogi) source
# viscous flow governed by a power-law viscosity with
# stress exponent of n = 3.
#
# run this example with
#
#   ./mogi.sh
#
# or with
#
#   ./mogi.sh --no-grd-output
#

WDIR=./mogi

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

time relax --no-stress-output --no-proj-output $* < mogi.input
