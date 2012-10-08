Deformation associated with the Mw 7.6 September 20th, 1999 Chi Chi earthquake. The coseismic slip distribution is from:

Y. Hsu, J.-P. Avouac, S.-B. Yu, C.-H. Chang, Y.-M. Wu, and J. Woessner, "Spatio-temporal slip, and stress level on the faults within the Western Foothills of Taiwan: Implications for fault frictional properties", Pure Appl. Geophys., 166 (2009)
http://mh-gps-p1.caltech.edu/~avouac/publications/Hus_GJI2009.pdf

To compute the coseismic displacement, run 

  ./coseismic.sh

Create a map of the coseismic slip displacements with

  grdmap.sh -b -100/100/-100/100 -p -2/2/0.01 -v 15 -e erpatch.sh -e ./ecoasts.sh -e ./efaults.sh coseismic/000

To compute the viscoelastic deformation, assuming a stratified viscoelastic earth with a viscoelastic substrate below 25km, run 

  ./visco1d.sh

To visualize in paraview, load 

  rfaults-001.vtp (color by slip)
  linear-layer-001.vtp
  vel-..vtk (time series of instantaneous velocity)
  power-..vtk (time series of power density)

You can compute the norm of the power density tensor with a Calculator Filter with the formula 

  sqrt(power_0^2+2*power_1^2+2*power_2^2+power_4^2+2*power_5^2+power_8^2)

then plot a contour of the norm to see where the source of the deformation occurs (Contour Filter)

To create a map of the viscoelastic relaxation, use, for example

  grdmap.sh -b -100/100/-100/100 -p -0.05/0.05/0.001 -v 0.05 -s 5 -e erpatch.sh -e ./ecoasts.sh -e ./efaults.sh visco1d/001-relax

To compute a stress-driven afterslip model, run

  ./afterslip.sh

Create a map of the afterslip solution with

  grdmap.sh -b -60/60/-40/80 -p -0.05/0.05/0.001 -v 0.20 -s 5 -e erpatch.sh -e ./ecoasts.sh -e ./efaults.sh -i gmt/relief_as_grad_div3.grd afterslip/009-relax

To localize the 3-d deformation in paraview, load

  rfaults-001.vtp
  aplane-001.vtp (afterslip plane)
  vel-..vtk

Our preferred model can be simulated with

./coupled.sh

a figure of the observed and modeled deformation can be generated with:

grdmap.sh -t 50 -b -110/110/-150/180 -p -0.2/0.2/0.01 -v 0.15 -s 10 -e ./efaults.sh -e ./ecoast.sh -e rpatch.sh -e ./econtour.sh -e ./epgps.sh -e ./eopts.sh -Y -0.6 -u m coupled/200-relax-up.grd
