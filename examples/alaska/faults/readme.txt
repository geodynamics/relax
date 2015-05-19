
# reference coordinates: N61.04 W147.73
echo -147.73 61.04 | proj +proj=utm +zone=6
460567.56       6767462.54

# change convention for slip model:
grep -v "#" holdahlsauber94.slip | join -o 2.5 -o 2.4 -o 1.2 -o 1.4 -o 2.2 -o 2.3 -o 2.6 -o 2.7 -o 2.8 - holdahlsauber94.geom | proj +proj=utm +zone=6 | awk '{$1=($1-474079.78);$2=($2-6768451.61);print $0}' | awk 'BEGIN{print "# nb slip(m) x1(km)    x2     x3 length(km) width strike dip   rake";pi=atan2(1,0)*2}{str=$5*pi/180;dip=$7/180*pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;n[1]=sin(str)*cos(dip);n[2]=-cos(str)*cos(dip);d[3]=sin(dip);ss=$3;ds=-$4;slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss);len=$6*2;wid=50000;x1=$1-s[1]*len/2-d[1]*wid/2;x2=$2-s[2]*len/2-d[2]*wid/2;x3=0;printf("%3.3d  %6.3f %6.2f %6.2f %5.2f      %5.2f %5.2f %6.1f %3.1f %6.3f\n",NR,slip,x1/1e3,x2/1e3,x3/1e3,len/1e3,wid/1e3,str*180/pi,dip*180/pi,rake*180/pi)}' > test.flt

# convert the Johnson et al (1996) model
grep -v "#" johnsonetal96.geom | awk '{print $4,$3,$2,$5,$6,$7,$8,$9,$10}' | proj +proj=utm +zone=6 | awk 'BEGIN{print "#n slip      x1       x2 x3 len wid str dip rake"}{printf("%2.2d %4.1f %7.2f %8.2f %2.0f %3.0f %3.0f %d  %2.0f  %3.0f\n",NR,$3,($2-6767462.54)/1e3,($1-460567.56)/1e3,$4,$5,$6,$7,$8,$9)}' > johnsonetal96.flt

# create coastline
pscoast -R-162/-138/56/64 -JM -Di -W1 -m | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | proj +proj=utm +zone=6 | awk '{if ("#"==substr($0,0,1)){print ">"}else{x=($1-460567.56)/1e3;y=($2-6767462.54)/1e3;print x,y,0}}' > coasts_km.xyz
xyz2vtk.sh coasts_km.xyz
