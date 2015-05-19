
echo 99.949 20.705 | proj +proj=utm +zone=47
598824.34       2289789.12


grep -v "#" Fault_slip_yellow.txt | proj +proj=utm +zone=47 | awk 'BEGIN{d2r=atan2(1,0)/90.0;L=1;W=0.6;strike=70*d2r;dip=86*d2r;s[1]=cos(strike);s[2]=sin(strike);s[3]=0;n[1]=sin(strike)*sin(dip);n[2]=-cos(strike)*sin(dip);n[3]=cos(dip);d[1]=sin(strike)*cos(dip);d[2]=-cos(strike)*cos(dip);d[3]=-sin(dip);print "# n      slip         x1         x2        x3 length width strike dip rake"}{printf "%03d %5.3e %6.3e %6.3e %5.3e      %d   %3.1f     %d  %d  %d\n", NR,$3,($2-2289789.12)/1e3-L/2*s[1]+W/2*d[1],($1-598824.34)/1e3-L/2*s[2]+W/2*d[2],$4+W/2*d[3],L,W,70,86,180}' | head
