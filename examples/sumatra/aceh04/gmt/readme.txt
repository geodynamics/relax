
# generate coast line data
pscoast -R90/100/0/20 -JM -Df -W1 -m | proj +proj=utm +zone=46 +datum=WGS84 +geoc | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | awk '{if (substr($0,1,1)=="*"){print ">"} else{x=($1-817174.63)/1e3;y=($2-369446.69)/1e3;print x,y,0}}' > coasts_km.xyz

# convert coast lines to Paraview
xyz2vtk.sh < coasts_km.xyz > ../paraview/coasts_km.vtp
