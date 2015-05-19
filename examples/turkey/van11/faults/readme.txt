
# USGS site
# http://earthquake.usgs.gov/earthquakes/eqinthenews/2011/usb0006bqc/#details

# reference point
echo 43.497 38.691 | proj +proj=utm +zone=38
369285.11       4283559.32

# convert to list for paste
cat van-seg1_43.6195_38.6261_253_40_63_18_25.txt | tr ' ' '\n' | awk '{if (""!=$0){print $0}}' > seg1_slip.dat
cat van-seg2_43.4057_38.6118_253_54_92_18_25.txt | tr ' ' '\n' | awk '{if (""!=$0){print $0}}' > seg2_slip.dat


# convert to Relax format (km)
echo 43.6195 38.6261 | proj +proj=utm +zone=38 | awk '{print 1,($2-4283559.32)/1e3,($1-369285.11)/1e3,0,18,38.9,253,40,0,1,1.080555,1,1}' | seg2flt.py | grep -v "#" | paste seg1_slip.dat - | awk 'BEGIN{print "# n      slip        x1        x2        x3    length     width    strike       dip      rake"}{s=$1;i=$2;$1=i;$2=s;printf "%03.3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > elliott+13_km.flt
echo 43.4057 38.6118 | proj +proj=utm +zone=38 | awk '{print 1,($2-4283559.32)/1e3,($1-369285.11)/1e3,0,18,30.9,253,54,0,1,0.858333,1,1}' | seg2flt.py | grep -v "#" | paste seg1_slip.dat - | awk 'BEGIN{print "# n      slip        x1        x2        x3    length     width    strike       dip      rake"}{s=$1;i=$2;$1=i;$2=s;printf "%03.3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",$1+648,$2,$3,$4,$5,$6,$7,$8,$9,$10}' >> elliott+13_km.flt


#
echo 43.6195 38.6261 | proj +proj=utm +zone=38 | awk '{print 1,($2-4283559.32),($1-369285.11),0,18e3,38.9e3,253,40,0,1e3,1.080555e3,1,1}' | seg2flt.py | grep -v "#" | paste seg1_slip.dat - | awk 'BEGIN{print "# n      slip        x1        x2        x3    length     width    strike       dip      rake"}{s=$1;i=$2;$1=i;$2=s;printf "%03.3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}' > elliott+13.flt
echo 43.4057 38.6118 | proj +proj=utm +zone=38 | awk '{print 1,($2-4283559.32),($1-369285.11),0,18e3,30.9e3,253,54,0,1e3,0.858333e3,1,1}' | seg2flt.py | grep -v "#" | paste seg2_slip.dat - | awk 'BEGIN{print "# n      slip        x1        x2        x3    length     width    strike       dip      rake"}{s=$1;i=$2;$1=i;$2=s;printf "%03.3d %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n",$1+648,$2,$3,$4,$5,$6,$7,$8,$9,$10}' >> elliott+13.flt



# coast line
pscoast -R40/50/35/45 -JM -Di -W1 -m | proj +proj=utm +zone=38 -t">" | awk '{if (substr($0,1,1)==">"){print $0} else{x=($1-369285.11)/1e3;y=($2-4283559.32)/1e3;print x,y,0}}' > coasts_km.xyz
