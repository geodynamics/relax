
# reference coordinates: 35.445°N, 81.392°E
echo 81.392 35.445 | proj +proj=utm +zone=44
535576.00       3922464.27

# X(utm m) Y(utm m) Strike (deg) Dip (deg) Rake (deg) Slip(m) Length(km) Top(km) Bottom(km)
#
# change convention for slip model:
awk '{print $5,$6}' fault_corners_yutian_ds.txt | paste fault_corners_yutian_ss.txt - | awk 'BEGIN{pi=atan2(1,0)*2}{ss=$6*cos($5/180*pi)+$11*cos($10/180*pi);ds=$6*sin($5/180*pi)+$11*sin($10/180*pi);r=atan2(ds,ss);s=sqrt(ss^2+ds^2);x2=($1-535576)/1e3;x1=($2-3922464.27)/1e3;strike=$3;dip=$4;len=$7;x3=$8;width=($9-$8)/sin(dip*pi/180);print NR,s,x1,x2,x3,len,width,strike,dip,r/pi*180}' > yutian.dat

# plot the coseismic displacements
map.sh -b -20/40/-40/40 -p -2/2/0.02 -v 2 -s 2 -e rpatch.sh yutian_co/000
