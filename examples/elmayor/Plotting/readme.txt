
# extract coast lines from GMT and convert to vtk
pscoast -R-119/-112/30/35 -JM -Di -W1 -m | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | proj +proj=utm +zone=11 | awk '{if ("#"==substr($0,0,1)){print ">"}else{x=($1-640915.61)/1e3;y=($2-3596850.34)/1e3;print x,y,0}}' > coasts_km.xyz
xyz2vtk.sh coasts_km.xyz

# extract political boundaries from GMT and convert to vtk
pscoast -R-119/-112/30/35 -JM -Di -Na -m | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | proj +proj=utm +zone=11 | awk '{if ("#"==substr($0,0,1)){print ">"}else{x=($1-640915.61)/1e3;y=($2-3596850.34)/1e3;print x,y,0}}' > boundary_km.xyz
xyz2vtk.sh boundary_km.xyz

# convert trace of El-Mayor Cucapah to vtk
cat cucapamayorfaults.lonlat | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | proj +proj=utm +zone=11 | awk '{if ("#"==substr($0,0,1)){print ">"}else{x=($1-640915.61)/1e3;y=($2-3596850.34)/1e3;print x,y,0}}' > cucapamayorfaults_km.xyz
xyz2vtk.sh cucapamayorfaults_km.xyz

# join the fit from matlab to database of station location
join -o 1.1,1.2,1.3,2.3,2.2,2.4 gps_ll.dat gps_post_fit.dat > gps_post_fit_ll.dat

# convert gps_ll.dat to _km.dat
awk '{print $2,$3,$1,$4,$5,$6}' gps_post_fit_200_aftershocks_ll.dat | proj +proj=utm +zone=11 | awk '{print $3,($1-640915.61)/1e3,($2-3596850.34)/1e3,$4,$5,$6}' > gps_post_fit_200_aftershocks_km.dat

# export gps to paraview
ofile=gps_post_fit_km.csv; echo \"east\(km\)\",\"north\(km\)\",\"height\(km\)\",\"ve\(km\)\",\"vn\(m\)\",\"vu\(m\)\" > $ofile; awk '{if ($3<400){print $2","$3",0,"$4","$5","$6}}' gps_post_fit_km.dat >> $ofile
