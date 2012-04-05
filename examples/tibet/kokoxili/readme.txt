
#
# co- and post-seismic deformation associated with the 2001, 
# Mw 7.8 Kokoxili earthquake
#
# the slip distributions of Wang et al. (2003) and Lasserre 
# et al. (2005) are in the folder ./faults
#
# location of some local GPS stations is in the folder ./gps
#


# The coseismic deformation can be modeled with the script
./kokoxili_co.sh

# then, to plot a cross-section of the shear stress across the fault, use
map.sh -b -50/50/-100/0 -p -1/1/0.01 -t 25 -u MPa kokoxili_co/001.op001-s12.grd

# to map the shear stress at 10 km depth, use
map.sh -b -300/300/-300/300 -p -1/1/0.01 -t 50 -u MPa kokoxili_co/000-s12.grd

# and to map the 3-d surface displacement, use
map.sh -b -250/250/-300/300 -p -0.2/0.2/0.01 -v 2 -s 20 -u m -e rpatch.sh kokoxili_co/000



# the viscoelastic relaxation using the layered model of 
#
#   Ryder et al. "Lower crustal relaxation beneath the Tibetan 
#   Plateau and Qaidam Basin following the 2001 Kokoxili earthquake"
#   Geophys. J. Int., 187, pp. 613-630, 2011
#
# can be computed with
./kokoxili_vs.sh

# the displacement is computed in map view and at the GPS 
# location of Ren and Wang, Quat. Sc., 25(1) 34-44, 2005.
# To compute the relaxation component (cumulative - coseismic),
# use
obsrelax.sh kokoxili_vs/????.txt


# the viscoelastic relaxation using a 3d visco-elastic
# model with a step in the depth of the brittle-ductile
# transition depth across the Kunlun fault can be modeled
# with
./kokoxili_vs3d.sh

# followed by
obsrelax.sh kokoxili_vs3d/????.txt

# a comparative plot of the layered and the 3-d viscous models
# can be obtained with
map.sh -b -260/340/-240/260 -p -0.01/0.01/0.0001 -v 0.03 -s 25 -e rpatch.sh -e ./eopts-relax.sh kokoxili_{vs,vs3d}/002-relax-up.grd
