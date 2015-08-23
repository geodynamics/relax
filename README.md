![Relax: It's not an acronym, it's a motto!](/graphics/Icon-Relax.png) 

**Relax: It's not an acronym, it's a motto!**

"RELAX - time-dependent postseismic deformation with afterslip and viscoelastic flow."

## INSTALLATION:

The code is written in Fortran90 and is optimized for the the gfortran and the 
INTEL ifort compiler. The gmt 4.5+ library is required to export results to GMT 
for post-processing.

## DOCUMENTATION:

Generate a .pdf file of the documentation with the command:

pdflatex latex/documentation.tex

Generate a browsable version of the code with:

doxygen .doxygen

## RUN:

Some examples are available in the examples directory. Look up the *.sh files for
comments and explanations.

## VISUALIZATION:

Many outputs are exported in the General Mapping Tools (GMT) format, deformation
maps can be obtained with typical GMT post-processing. Check the post-processing
and visualization tools available in the 'util' directory.
Make sure to have the .ps file viewer 'gv' or the .pdf file viewer 'xpdf' installed.
Simulations can be visualized in 3D with the free software Paraview (paraview.org).

## PURPOSE:

RELAX computes nonlinear time-dependent viscoelastic deformation with 
powerlaw rheology and rate-strengthening friction in a cubic grid due to coseismic 
stress changes, initial stress, surface loads, and/or moving faults.

## DESCRIPTION:

Computation is done semi-analytically inside a cartesian grid. The grid is defined
by its size sx1*sx2*sx3 and the sampling intervals dx1, dx2 and dx3. Rule of thumb
is to allow for at least five samples per fault length or width, and to have the 
tip of any fault at least 10 fault widths away from any edge of the computational
grid.

Coseismic stress changes and initial coseismic deformation results from the 
presence of dislocations in the brittle layer. Fault geometry is prescribed 
following Okada or Wang's convention, with the usual slip, strike, dip and rake and
is converted to a double-couple equivalent body-force analytically. Current 
implementation allows shear fault (strike slip and dip slip), dykes, Mogi source, 
and surface traction. Faults and dykes can be of arbitrary orientation in the half
space.


## INPUT:

Static dislocation sources are discretized into a series of planar segments. Slip
patches are defined in terms of position, orientation, and slip, as illustrated in
the following figure:

                     N (x1)
                    /
                   /| Strike
       x1,x2,x3 ->@------------------------      (x2)
                  |\        p .            \ W
                  :-\      i .              \ i
                  |  \    l .                \ d
                  :90 \  S .                  \ t
                  |-Dip\  .                    \ h
                  :     \. | Rake               \
                  |      -------------------------
                  :             L e n g t h
                  Z (x3)

Dislocations are converted to double-couple equivalent body-force analytically.
Solution displacement is obtained by application of the Green's functions in the 
Fourier domain.

For friction faults where slip rates are evaluated from stress and a constitutive 
law, the rake corresponds to the orientation of slip. That is, if r_i is the rake
vector and v_i is the instantaneous velocity vector, then r_j v_j >= 0. 

## REFERENCES:

More information about parameters and constitutive laws can be found in

S. Barbot and Fialko Y., "Fourier-Domain Green's Function for an Elastic Semi-
Infinite Solid under Gravity, with Applications to Earthquake and Volcano 
Deformation", Geophysical Journal International, v. 182, no. 2, pp. 568-582, 2010,
doi:10.1111/j.1365-246X.2010.04655.x

and

S. Barbot and Fialko Y., "A Unified Continuum Representation of Postseismic 
Relaxation Mechanisms: Semi-Analytic Models of Afterslip, Poroelastic Rebound and
Viscoelastic Flow", Geophysical Journal International, v. 182, 3, p. 1124-1140, 
2010, doi:10.1111/j.1365-246X.2010.04678.x

Please cite these papers in publications or public presentations when referring
to this method.


