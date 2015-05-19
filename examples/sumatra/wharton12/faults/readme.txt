
# convert slip model to relax format
awk 'NR>3' chlieh+07.dat | proj +proj=utm +zone=47 | awk 'BEGIN{print "# n      slip         x1         x2        x3 length width strike dip    rake"}{printf "%3.3d %6.3e %10.3e %10.3e %6.3e     %d    %d    %d  %d %7.3f\n",NR,$3/1e2,($2-221910.44)/1e3,($1+56748.29)/1e3,$7,20,16,$5,$6,$4}' > chlieh+07_km.dat

# convert to paraview format
flt2vtk.sh < chlieh+07_km.flt > ../paraview/chlieh+07_km.vtp
