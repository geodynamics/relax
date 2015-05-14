
# reference point
echo -70.0600 -22.3400 | proj +proj=utm +zone=19
390844.18       -2470844.80

# convert the srcmod .fsp source distribution of Bejar-Pizzaro et al. (2010) to the .flt format
grep -v "#" s2007TOCOPI01BEJA.dat | proj +proj=utm +zone=19 -r | awk 'BEGIN{Pi=3.14159265;print "#  n   slip       x1       x2       x3 length width strike  dip  rake"}{str0=$10;dip0=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%04.4d %6.3f %8.3f %8.3f %8.3f %6.1f %5.1f %6.1f %4.1f %5.1f\n",NR,$6,($2+2470844.80)/1e3-$8*s[1]/2,($1-390844.18)/1e3-$8*s[2]/2,$5,$8,$9,$10,$11,$7}' > bejar-pizzaro+10_km.flt
