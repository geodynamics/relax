
# reference point
echo 84.708 28.147 | proj +proj=utm +zone=45
274917.78       3115610.96

# convert the model of Wei (2015) to the .flt format
grep -v "#" static_out.dat | proj +proj=utm +zone=45 -r | awk 'BEGIN{Pi=3.14159265;print "# n    slip         x1         x2        x3   length  width strike   dip  rake"}{slip=$4/1e2;rake=$5;str0=$6;dip0=$7;L=$10;W=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2-3115610.96)/1e3-s[1]*L/2;x2=($1-274917.78)/1e3-s[2]*L/2;x3=$3;printf "%03.3d %7.3f %10.5f %10.5f %9.5f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,slip,x1,x2,x3,L,W,str0,dip0,rake}' > wei+15_km_1.flt
