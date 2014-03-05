Viscoelastic deformation in the Mojave Desert due to the 1992 Landers and the 1999 Hector Mine earthquakes.

The script landers.sh compute the postseismic viscoelastic relaxation following the Landers earthquake, based on the model of Fialko (2004). 

The script hectorm.sh computes the relaxation with the same Earth model due to the Hector Mine earthquake, using the model of Simons, Fialko and Rivera (2002). The output is exported every 0.033 year in files prefixed by 000,001,002... corresponding to t=0,0.033,0.066 yr, etc... 

For a sanity check, plot the coseismic displacement. For example - assuming erpatch.sh is in a directory listed in the environment variable $PATH - you can use the command:

  grdmap.sh -b -60/60/-60/60 -e ./efaults.sh hectorm/000

You can also plot the cumulative postseismic displacements after 3 years:

  grdmap.sh -b -100/100/-100/100 -s 5 -e ./efaults.sh hectorm/090-relax

To compare with insar, you can convert the three-component displacements to line-of-sight with the command:

  grdinsar.sh hectorm/090-relax

and plot the line-of-sight model with the command:

  grdmap.sh -b -150/150/-150/150 -e ./efaults.sh hectorm/090-relax-los.grd

Origin of coordinates system: -116.27E, 34.595N and 566940.91E, 3828373.73N in UTM zone 11.

