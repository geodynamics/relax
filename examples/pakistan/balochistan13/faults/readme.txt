# M7.7 - 61km NNE of Awaran, Pakistan
# http://comcat.cr.usgs.gov/earthquakes/eventpage/usb000jyiv#summary

# reference point
echo 65.501 26.951 | proj +proj=utm +zone=41
748291.43       2983465.15

# convert the .dat format to .flt format
 grep -v "#" static_out.txt | proj +proj=utm +zone=41 -r | awk 'BEGIN{Pi=3.14159265;print "#  n  slip       x1       x2       x3   length  width strike  dip  rake"}{slip=$4/1e2;rake=$5;str0=$6;dip0=$7;L=$10;W=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2-2983465.15)/1e3-s[1]*L/2;x2=($1-748291.43)/1e3-s[2]*L/2;x3=$3;printf "%3d %6.2f %8.3f %8.3f %8.3f   %6.2f %6.2f    %d   %d %5.1f\n", NR,slip,x1,x2,x3,L,W,str0,dip0,rake}' > avouac+14_km.flt
grep -v "#" static_out.txt | proj +proj=utm +zone=41 -r | awk 'BEGIN{Pi=3.14159265;print "# n   slip         x1         x2        x3   length  width strike   dip  rake"}{slip=$4/1e2;rake=$5;str0=$6;dip0=$7;L=$10;W=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2-2983465.15)/1e3-s[1]*L/2;x2=($1-748291.43)/1e3-s[2]*L/2;x3=$3;printf "%3d %6.2f %10.2f %10.2f %9.1f %8.1f %6.1f %6.2f %5.1f %5.1f\n", NR,slip,x1*1e3,x2*1e3,x3*1e3,L*1e3,W*1e3,str0,dip0,rake}' > avouac+14.flt
