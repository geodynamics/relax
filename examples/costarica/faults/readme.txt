
# reference point
echo -85.56 9.76 | proj +proj=utm +zone=16
657947.74       1079213.96

# convert raw data to .flt format
grep -v "#" FFM_model.txt | awk '{print $4,$3,$5,$6,$7,$8,$10}' | proj +proj=utm +zone=16 | awk 'BEGIN{print "# n   slip            x1            x2            x3       length        width strike    dip    rake"}{printf "%03d %6.3f %13.6e %13.6e %13.6e %12.5e %12.5e %6.2f %6.2f %7.2f\n", NR,$7,($2-1079213.96)/1e3,($1-657947.74)/1e3,$3,7.5,7.39303,$4,$5,$6}' > yue+13_km.flt
