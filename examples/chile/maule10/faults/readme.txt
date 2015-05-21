
# reference point
echo -72.733  -35.909 | proj +proj=utm +zone=18
704570.60       -3976229.63

# convert the slip distribution of lorito et al. (2001) to .flt format
grep -v "#" lorito+11.dat | awk '{$1=""; print $0}' | proj +proj=utm +zone=18 | awk 'BEGIN{Pi=3.14159265;print "# n   slip       x1       x2       x3 length width strike  dip  rake"}{str0=$4;dip0=$5;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%03.3d %6.3f %8.3f %8.3f %8.3f %6.1f %5.1f %6.1f %4.1f %5.1f\n",NR,$6,(($2+3976229.63)/1e3-25/2*s[1]+25/2*d[1]),(($1-704570.60)/1e3-25/2*s[2]+25/2*d[2]),$3,25,25,$4,$5,$8;}' > lorito+11_km.flt
grep -v "#" lorito+11.dat | awk '{$1=""; print $0}' | proj +proj=utm +zone=18 | awk 'BEGIN{Pi=3.14159265;print "# n   slip       x1       x2       x3 length width strike  dip  rake"}{str0=$4;dip0=$5;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%03.3d %6.3f %12.3f %12.3f %12.3f %6.1f %5.1f %6.1f %4.1f %5.1f\n",NR,$6,(($2+3976229.63)/1e3-25/2*s[1]+25/2*d[1])*1e3,(($1-704570.60)/1e3-25/2*s[2]+25/2*d[2])*1e3,$3*1e3,25*1e3,25*1e3,$4,$5,$8;}' > lorito+11.flt

# convert the srcmod model of Luttrel et al (2011) to .flt format
grep -v "#" s2010MAULEC01LUTT.dat | proj +proj=utm +zone=18 -r | awk 'BEGIN{Pi=3.14159265;print "#  n   slip       x1       x2       x3 length width strike  dip  rake"}{str0=$10;dip0=$11;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip); printf "%04.4d %6.3f %8.3f %8.3f %8.3f %6.1f %5.1f %6.1f %4.1f %5.1f\n",NR,$6,($2+3976229.63)/1e3-$8*s[1]/2,($1-704570.60)/1e3-$8*s[2]/2,$5,$8,$9,$10,$11,90}' > luttrel+11_km.flt
