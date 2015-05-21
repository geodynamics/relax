
# reference coordinates: W120.370 N35.8150, UTM zone 11
echo -120.370 35.8150 | proj +proj=utm +zone=11
195506.33       3968673.17

# convert Ziv (2012) model to relax format
echo -120.365997 35.8181992 | proj +proj=utm +zone=11 | awk '{print ($1-195506.33)/1e3+6.5,($2-3968673.17)/1e3-7.5}'
6.87404 -7.15744
paste quad4okada best_slip | awk 'BEGIN{pi=atan2(1,0)*2;strike=319.5;cs=cos(strike/180*pi);ss=sin(strike/180*pi); print "# i  slip    x1     x2    x3 length width strike dip rake"}{x1i=-7.15744;x2i=6.87404;x1=x1i+($1-$3/2)*cs;x2=x2i+($1-$3/2)*ss;x3=$2-$4/2;len=$3;wid=$4;rake=180;slip=$7;dip=90;printf("%3.3i %5.3f %5.2f %6.2f %5.2f %6.2f %5.2f %6.1f %3.2i %4.3i\n",NR,slip,x1,x2,x3,len,wid,strike,dip,rake)}' > ziv12.dat

paste quad4okada best_slip | awk 'BEGIN{pi=atan2(1,0)*2;strike=319.5;cs=cos(strike/180*pi);ss=sin(strike/180*pi); print "# i    x1     x2    x3 length width strike dip rake"}{x1i=-8;x2i=6.5;x1=x1i+($1-$3/2)*cs;x2=x2i+($1-$3/2)*ss;x3=$2-$4/2;len=$3;wid=$4;rake=180;slip=$7;dip=90;printf("%3.3i %5.2f %6.2f %5.2f %6.2f %5.2f %6.1f %3.2i %4.3i\n",NR,x1,x2,x3,len,wid,strike,dip,rake)}' > ziv12.flt


# rotate a coseismic slip distribution
grep -v "#" $FLT | awk 'BEGIN{x1c=5.77;x2c=-3.69 ;pi=atan2(1,0)*2;r=5;rot=r*pi/180;}{\
                                                                x1=+($3-x1c)*cos(rot)-($4-x2c)*sin(rot);\
                                                                x2=+($3-x1c)*sin(rot)+($4-x2c)*cos(rot);\
                                                                $3=x1c+x1;$4=x2c+x2;$8=$8+r;print $0}'

# convert Ji (2004) to the .flt format
grep -v "#" s2004PARKFI01JIxx.dat | proj +proj=utm +zone=11 -r | awk 'BEGIN{Pi=3.14159265;print "# n   slip       x1       x2       x3 length width strike  dip  rake"}{str0=$10;dip0=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%03.3d %6.3f %8.3f %8.3f %8.3f %6.2f %5.2f %6.1f %4.1f %5.1f\n",NR,$6,($2-3968673.17)/1e3-$8*s[1]/2,($1-195506.33)/1e3-$8*s[2]/2,$5,$8,$9,$10,$11,$7}' > ji04_km.flt
grep -v "#" s2004PARKFI01JIxx.dat | proj +proj=utm +zone=11 -r | awk 'BEGIN{Pi=3.14159265;print "# n   slip         x1         x2         x3  length   width strike  dip  rake"}{str0=$10;dip0=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%03.3d %6.3f %10.3f %10.3f %10.3f %7.2f %7.2f %6.1f %4.1f %5.1f\n",NR,$6,(($2-3968673.17)/1e3-$8*s[1]/2)*1e3,(($1-195506.33)/1e3-$8*s[2]/2)*1e3,$5*1e3,$8*1e3,$9*1e3,$10,$11,$7}' > ji04.flt

# convert Dregert et al. (2004) to .flt format
grep -v "#" s2004PARKFI01DREG.dat | proj +proj=utm +zone=11 -r | awk 'BEGIN{Pi=3.14159265;print "# n   slip       x1       x2      x3   length    width strike  dip  rake"}{str0=$9;dip0=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%03.3d %6.3f %8.3f %8.3f %7.4f %8.3f %8.3f %6.1f %4.1f %5.1f\n",NR,$6,($2-3968673.17)/1e3-$7*s[1]/2,($1-195506.33)/1e3-$7*s[2]/2,$5,$7,$8,$9,$10,$11}' > dreger+05_km.flt
grep -v "#" s2004PARKFI01DREG.dat | proj +proj=utm +zone=11 -r | awk 'BEGIN{Pi=3.14159265;print "# n   slip         x1         x2      x3    length     width strike  dip  rake"}{str0=$9;dip0=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%03.3d %6.3f %10.3f %10.3f %7.1f %9.3f %9.3f %6.1f %4.1f %5.1f\n",NR,$6,(($2-3968673.17)/1e3-$7*s[1]/2)*1e3,(($1-195506.33)/1e3-$7*s[2]/2)*1e3,$5*1e3,$7*1e3,$8*1e3,$9,$10,$11}' > dreger+05.flt
