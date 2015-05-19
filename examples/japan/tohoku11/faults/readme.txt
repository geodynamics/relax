
# reference point
echo 142.369 38.322 | proj +proj=utm +zone=54
619669.75       4242429.17

# convert the Caltech model to Relax format
grep -v "#" subfault.dat | proj +proj=utm +zone=54 | awk '{$1=($1-619669.75)/1e3;$2=($2-4242429.17)/1e3;print $0}' | awk 'BEGIN{print "# nb slip(m) x1(km)    x2     x3 length(km) width strike dip   rake";pi=atan2(1,0)*2}{str=$6*pi/180;dip=$7/180*pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;n[1]=sin(str)*cos(dip);n[2]=-cos(str)*cos(dip);d[3]=sin(dip);slip=$4;rake=$5;len=25;wid=20;x1=$2-s[1]*len/2-d[1]*wid/2;x2=$1-s[2]*len/2-d[2]*wid/2;x3=$3-s[3]*len-d[3]*wid/2;printf("%3.3d  %6.3f %6.2f %6.2f %5.2f      %5.2f %5.2f %6.1f %3.1f %6.2f\n",NR,slip,x1,x2,x3,len,wid,str*180/pi,dip*180/pi,rake)}' > wei+11.flt

# create Japan coast lines
pscoast -R138/148/32/46 -JM -Di -W1 -m | awk '{if (substr($0,0,1)==">"){print "#"}else{print $0}}' | proj +proj=utm +zone=54 | awk '{if ("#"==substr($0,0,1)){print ">"}else{x=($1-619669.75)/1e3;y=($2-4242429.17)/1e3;print x,y,0}}' > coasts_km.xyz

# create slab geometry
ofile=kur_slab1.0_clip_100km.csv;echo \"East\",\"North\",\"Depth\",\"Easting\",\"Northing\",\"Elevation \(m\)\" > $ofile; grep -v "NaN" ~/Documents/work/taiwan/slab/kurils/kur_slab1.0_clip.xyz | proj +proj=utm +zone=54 | awk '{x=($1-619669.75)/1e3;y=($2-4242429.17)/1e3;if ($3>-100 && x<600){printf("%f,%f,%f,%f,%f,%f\n",x,y,-$3,x,y,$3*1e3)}}' >> $ofile

ofile=ryu_slab1.0_clip_100km.csv;echo \"East\",\"North\",\"Depth\",\"Easting\",\"Northing\",\"Elevation \(m\)\" > $ofile; grep -v "NaN" ~/Documents/work/taiwan/slab/ryukyu/ryu_slab1.0_clip.xyz | proj +proj=utm +zone=54 | awk '{x=($1-619669.75)/1e3;y=($2-4242429.17)/1e3;if ($3>-100 && x<100){printf("%f,%f,%f,%f,%f,%f\n",x,y,-$3,x,y,$3*1e3)}}' >> $ofile
