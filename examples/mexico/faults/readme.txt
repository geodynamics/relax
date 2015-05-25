
# slip distributions from
#
# M. Radiguet, F. Cotton, M. Vergnolle, M. Campillo, A. Walpersdorf, N. Cotte, and V. Kostoglodov. 
# Slow slip events and strain accumulation in the Guerrero gap, Mexico. J. Geophys. Res., 117(B4) :B04305, 2012. 
#

grep -v "#" data_SSE_2006.txt | awk 'function cosd(x){return cos(x*atan2(1,0)/90.0)};function sind(x){return sin(x*atan2(1,0)/90.0)}BEGIN{print "# n       slip         x1         x2        x3 length width strike dip rake"}{L=$6;W=$7;strike=$8;dip=-$9;rake=$10;s[1]=cosd(strike);s[2]=sind(strike);s[3]=0;n[1]=sind(strike)*sind(dip);n[2]=-cosd(strike)*sind(dip);n[3]=cosd(dip);d[1]=sind(strike)*cosd(dip);d[2]=-cosd(strike)*cosd(dip);d[3]=-sind(dip);printf "%03d %10.3e %10.3e %10.3e %9.3e   %3.1f  %3.1f  %3.1f %3.1f %3.1f\n", NR,$2/1e3,$3+L*s[1]-W*d[1],$4+L*s[2]-W*d[2],$5-W*d[3],L,W,strike-180,-dip,rake}' > radiguet+12_06_km.flt

grep -v "#" data_SSE_2006.txt | awk 'function cosd(x){return cos(x*atan2(1,0)/90.0)};function sind(x){return sin(x*atan2(1,0)/90.0)}BEGIN{print "# n       slip            x1            x2           x3       length       width strike dip rake"}{L=$6;W=$7;strike=$8;dip=-$9;rake=$10;s[1]=cosd(strike);s[2]=sind(strike);s[3]=0;n[1]=sind(strike)*sind(dip);n[2]=-cosd(strike)*sind(dip);n[3]=cosd(dip);d[1]=sind(strike)*cosd(dip);d[2]=-cosd(strike)*cosd(dip);d[3]=-sind(dip);printf "%03d %10.3e %13.6e %13.6e %12.6e   %6.4f  %6.4f  %3.1f %3.1f %3.1f\n", NR,$2/1e3,($3+L*s[1]-W*d[1])*1e3,($4+L*s[2]-W*d[2])*1e3,($5-W*d[3])*1e3,L*1e3,W*1e3,strike-180,-dip,rake}' > radiguet+12_06.flt
