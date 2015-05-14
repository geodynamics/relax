
# reference point
echo 94 2 | proj +proj=utm +zone=47
-56748.29       221910.44

# compute deviatoric stress tau from individual stress components
wdir=coseismic;index=000.op001;grdmath 5 9 DIV $wdir/$index-s11.grd 2 POW MUL -8 9 DIV $wdir/$index-s11.grd $wdir/$index-s22.grd MUL MUL ADD -2 9 DIV $wdir/$index-s11.grd $wdir/$index-s33.grd MUL MUL ADD 5 9 DIV $wdir/$index-s22.grd 2 POW MUL ADD -2 9 DIV $wdir/$index-s22.grd $wdir/$index-s33.grd MUL MUL ADD 11 9 DIV $wdir/$index-s33.grd 2 POW MUL ADD 2 $wdir/$index-s12.grd 2 POW MUL ADD 2 $wdir/$index-s13.grd 2 POW MUL ADD 2 $wdir/$index-s23.grd 2 POW MUL ADD = $wdir/$index-tau.grd

# get geographic bounds of Relax models
grd2xyz ./coseismic/000-east.grd | awk '{printf "%f %f\n", $1*1e3-56748.29,$2*1e3+221910.44}' | invproj +proj=utm +zone=47 -f '%f' | minmax
N = 65536      <82.28486/105.575262>   <-9.571978/13.495301>

# project the gradient of topography to the local coordinate system
grd2xyz -R80/105/-10/15 /Volumes/data/shared/Topography/GLOBAL/ETOPO1/ETOPO1_Bed_g.grad | proj +proj=utm +zone=47 | awk '{print ($1+56748.29)/1e3, ($2-221910.44)/1e3, $3}' | xyz2grd -R-1280/1270/-1280/1270 -I10/10 -Ggmt/topo_hs.grd
