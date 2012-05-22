
# reference point
echo 142.369 38.322 | proj +proj=utm +zone=54
619669.75       4242429.17

# create Japan coast lines
pscoast -R138/148/32/46 -JM -Di -W1 -m | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | proj +proj=utm +zone=54 | awk '{if ("#"==substr($0,0,1)){print ">"}else{x=($1-619669.75)/1e3;y=($2-4242429.17)/1e3;print x,y,0}}' > coasts_km.xyz

# create slab geometry
ofile=kur_slab1.0_clip_100km.csv;echo \"East\",\"North\",\"Depth\",\"Easting\",\"Northing\",\"Elevation \(m\)\" > $ofile; grep -v "NaN" ~/Documents/work/taiwan/slab/kurils/kur_slab1.0_clip.xyz | proj +proj=utm +zone=54 | awk '{x=($1-619669.75)/1e3;y=($2-4242429.17)/1e3;if ($3>-100 && x<600){printf("%f,%f,%f,%f,%f,%f\n",x,y,-$3,x,y,$3*1e3)}}' >> $ofile

# run the coseismic model
./coseismic.sh

# map view of the surface displacements
map.sh -b -350/250/-400/400 -e rpatch.sh -e ./gmt/ecoast.sh coseismic/000
