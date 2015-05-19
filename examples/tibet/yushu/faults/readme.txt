
# reference coordinates: 33.271°N 96.629°E
echo 96.629 33.271 | proj +proj=utm +zone=47
279165.36       3683838.25

# change convention for slip model:
awk '{print $5,$6}' fault_corners_yushu_ds.txt | paste fault_corners_yushu_ss.txt - | awk 'BEGIN{pi=atan2(1,0)*2}{ss=$6*cos($5/180*pi)+$11*cos($10/180*pi);ds=$6*sin($5/180*pi)+$11*sin($10/180*pi);r=atan2(ds,ss);s=sqrt(ss^2+ds^2);x2=($1-279165.36)/1e3;x1=($2-3683838.25)/1e3;strike=$3;dip=$4;len=$7;x3=$8;width=($9-$8)/sin(dip*pi/180);print NR,s,x1,x2,x3,len,width,strike,dip,r/pi*180}' > yushu.dat
