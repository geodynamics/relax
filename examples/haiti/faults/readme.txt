
echo -73.4992 18.0 | proj +proj=utm +zone=18
658893.42       1990828.72

# convert slip distribution of Symithe et al. (2013) to Relax .flt format
grep -v "#" slip.sjs2 | awk 'BEGIN{pi=atan2(1,0)*2;print "#  n      slip         x1         x2        x3    length     width    strike     dip    rake"}{printf "%04.4d %9.4f %10.3f %10.3f %9.3f %9.3f %9.3f %9.3f %7.3f %7.3f\n", NR,sqrt(($8)^2+($9)^2),$7,$6,$3,$1,$2,$5,$4,atan2($9,$8)*180/pi}'

# convert the slip distribution of Calais et al., (2010) from SRCMOD to Relax .flt format
grep -v "#" s2010HAITIx01CALA.fsp | proj +proj=utm +zone=18 -r | awk 'BEGIN{Pi=atan2(1,0)*2;print "# n   slip      x1      x2      x3   length  width    strike    dip  rake"}{L=$8;W=$9;slip=$6;rake=$7;str0=$10;dip0=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2-1990828.72)/1e3-s[1]*L/2;x2=($1-658893.42)/1e3-s[2]*L/2;x3=$5;printf "%03.3d %6.3f %7.2f %7.2f %7.2f   %6.2f %6.2f    %6.2f %6.2f %5.1f\n", NR,slip,x1,x2,x3,L,W,str0,dip0,rake}' > calais+10_km.flt

# convert the slip distribution of Sladen et al., (2010) from SRCMOD to Relax .flt format
grep -v "#" s2010HAITIx01SLAD.fsp | proj +proj=utm +zone=18 -r | awk 'BEGIN{Pi=atan2(1,0)*2;print "# n   slip      x1      x2      x3   length  width    strike    dip  rake"}{L=$8;W=$9;slip=$6;rake=$7;str0=$10;dip0=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);x1=($2-1990828.72)/1e3-s[1]*L/2;x2=($1-658893.42)/1e3-s[2]*L/2;x3=$5;printf "%03.3d %6.3f %7.2f %7.2f %7.2f   %6.2f %6.2f    %6.2f %6.2f %5.1f\n", NR,slip,x1,x2,x3,L,W,str0,dip0,rake}' > sladen+10_km.flt
