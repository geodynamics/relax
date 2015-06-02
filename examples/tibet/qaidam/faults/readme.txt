
# reference coordinates: 37.589°N, 95.836°E
echo 95.836 37.589 | proj +proj=utm +zone=46
750402.82       4163997.29

# change convention for slip model:
grep -v "#" fault_corners_qaidam08.txt | awk 'BEGIN{pi=atan2(1,0)*2}{x2=($1-750402.82)/1e3;x1=($2-4163997.29)/1e3;strike=$3;dip=$4;r=$5;s=$6;len=$7;x3=$8;width=($9-$8)/sin(dip*pi/180);print NR,s,x1,x2,x3,len,width,strike,dip,r}' > qaidam08_km.flt

# plot the coseismic displacement
map.sh -b -40/40/-50/50 -s 2 -e rpatch.sh qaidam08_co/000
