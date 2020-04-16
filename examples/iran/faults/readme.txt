
# reference point
echo 58.4 29 | proj +proj=utm +zone=40
636367.15       3208793.43

# convert to .flt format
grep -v "#" fialko+05.dat | awk 'NR>13{pi=2*atan2(1,0);slip=sqrt($9**2+$10**2);rake=atan2($10,$9)*180/pi;$NF="";$1=NR-13; w=$4-$3; l=$2*2; x3=$3; dip=$5+15; str=$6; sv1=cos(str/180*pi);sv2=sin(str/180*pi); x1=$7-sv1*l/2; x2=$8-sv2*l/2; print $1,slip,x1,x2,x3,l,w,str,dip,rake}' | awk 'NR<=286'> fialko+05_km.flt
