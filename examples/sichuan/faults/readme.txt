# Mw 7.9 2008 Sichuan earthquake

# reference point used by Mong-Han Huang
echo 103.501 32.501 | proj +proj=utm +zone=48
359179.90       3596959.89

# reference point used by Xiaopeng
echo 104.2 31.4 | proj +proj=utm +zone=48
423946.66       3474210.00

# difference
64766.76  -122749.89

# convert source?.dat to source?.flt
grep -v "#" source1.dat | awk -v xo=100.37 -v yo=120.91 -v str0=-128 -v dip0=70 -v L=95.73 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y-122.74989;x2=xo+s[2]*(x-L/2)+d[2]*y+64.76676;x3=+d[3]*y;printf "%3d %6.2f %7.2f %7.2f %7.2f   %6.2f %6.2f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source1.flt
grep -v "#" source2.dat | awk -v xo=9.91 -v yo=34.70 -v str0=-137 -v dip0=50 -v L=154.82 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y-122.74989;x2=xo+s[2]*(x-L/2)+d[2]*y+64.76676;x3=+d[3]*y;printf "%3d %6.2f %7.2f %7.2f %7.2f   %6.2f %6.2f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source2.flt
grep -v "#" source3.dat | awk -v xo=-77.35 -v yo=-48.99 -v str0=-128 -v dip0=35 -v L=87.6 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y-122.74989;x2=xo+s[2]*(x-L/2)+d[2]*y+64.76676;x3=+d[3]*y;printf "%3d %6.2f %7.2f %7.2f %7.2f   %6.2f %6.2f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source3.flt
grep -v "#" source4.dat | awk -v xo=-11.36 -v yo=-2.31 -v str0=-137 -v dip0=25 -v L=82.07 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y-122.74989;x2=xo+s[2]*(x-L/2)+d[2]*y+64.76676;x3=+d[3]*y;printf "%3d %6.2f %7.2f %7.2f %7.2f   %6.2f %6.2f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source4.flt

# merge all segments
cat source?.flt | grep -v "#" | awk 'BEGIN{print "# n    slip       x1        x2        x3     length    width strike dip  rake"}{$1=sprintf("%3d",NR);printf "%3d %7.3f %8.4f %9.4f %9.4f   %8.4f %8.4f  %d   %d %5.1f\n", $1,$2/1e2,$3,$4,$5,$6,$7,$8,$9,$10}' > tong+10_km.flt

grep -v "#" tong+10_km.flt | awk '{$3=$3*1e3;$4=$4*1e3;$5=$5*1e3;$6=$6*1e3;$7=$7*1e3;printf "%3d %7.3f %12.1f %12.1f %12.1f   %11.1f %11.1f  %d   %d %5.1f\n", $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > tong+10.flt

# export Shen et al to relax
grep -v "#" Shen.prn | awk '{s=sqrt($9*$9+$10*$10);r=atan2($10,$9)*90/atan2(1,0);print $8,$7,$4,s,$2,$3,$6,$5,r}' | proj +proj=utm +zone=48 | awk 'BEGIN{print "# nb   slip      x1      x2      x3   length  width strike dip  rake";Pi=3.14159265}{str=$7*Pi/180;dip=$8/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);l=$5;width=$6;x1=($2-3596959.89)/1e3+width*d[1];x2=($1-359179.90)/1e3+width*d[2];x3=$3+width*d[3];printf "%3d %6.2f %7.2f %7.2f %7.2f   %6.2f %6.2f  %d   %d %5.1f\n", NR,$4,x1,x2,x3,l,width,$7,$8,$9; print d[1],d[2],d[3]}' > shen+000.flt


grep -v "#" source1.dat | awk -v xo=100.37 -v yo=120.91 -v str0=-128 -v dip0=70 -v L=95.73 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y;x2=xo+s[2]*(x-L/2)+d[2]*y;x3=+d[3]*y;printf "%3d %6.2f %9.4f %9.4f %9.4f   %8.4f %8.4f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source1.flt
grep -v "#" source2.dat | awk -v xo=9.91 -v yo=34.70 -v str0=-137 -v dip0=50 -v L=154.82 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y;x2=xo+s[2]*(x-L/2)+d[2]*y;x3=+d[3]*y;printf "%3d %6.2f %9.4f %9.4f %9.4f   %8.4f %8.4f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source2.flt
grep -v "#" source3.dat | awk -v xo=-77.35 -v yo=-48.99 -v str0=-128 -v dip0=35 -v L=87.6 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y;x2=xo+s[2]*(x-L/2)+d[2]*y;x3=+d[3]*y;printf "%3d %6.2f %9.4f %9.4f %9.4f   %8.4f %8.4f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source3.flt
grep -v "#" source4.dat | awk -v xo=-11.36 -v yo=-2.31 -v str0=-137 -v dip0=25 -v L=82.07 'BEGIN{Pi=3.14159265;print "# nb  slip      x1      x2      x3   length  width strike dip  rake"}{if (1==NR){L=L+2*$1};x=$1;y=$5;len=$3-$1;width=$5-$6;ss=$9;ds=$10;str=str0*Pi/180;dip=dip0/180*Pi;s[1]=cos(str);s[2]=sin(str);s[3]=0;d[1]=sin(str)*cos(dip);d[2]=-cos(str)*cos(dip);d[3]=-sin(dip);slip=sqrt(ss^2+ds^2);rake=atan2(ds,ss)*180/Pi;x1=yo+s[1]*(x-L/2)+d[1]*y;x2=xo+s[2]*(x-L/2)+d[2]*y;x3=+d[3]*y;printf "%3d %6.2f %9.4f %9.4f %9.4f   %8.4f %8.4f  %d   %d %5.1f\n", NR,slip,x1,x2,x3,len,width,str0,dip0,rake}' > source4.flt

