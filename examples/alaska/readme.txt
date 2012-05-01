

# reference coordinates: N61.04 W147.73
echo -147.73 61.04 | proj +proj=utm +zone=6
460567.56       6767462.54

# create coastline
pscoast -R-162/-138/56/64 -JM -Di -W1 -m | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | proj +proj=utm +zone=6 | awk '{if ("#"==substr($0,0,1)){print ">"}else{x=($1-460567.56)/1e3;y=($2-6767462.54)/1e3;print x,y,0}}' > coasts_km.xyz
xyz2vtk.sh coasts_km.xyz

# compute the coseismic displacement
./coseismic.sh

# plot the coseismic surface displacements
map.sh -b -700/500/-700/400 -s 40 -c anatolia.cpt -e ./gmt/ecoast.sh -e rpatch.sh coseismic/000
