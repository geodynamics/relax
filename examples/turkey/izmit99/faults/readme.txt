
# reference point
echo 29.987 40.702 | proj +proj=utm +zone=36
245446.49       4510043.99

# create _km relax format 
grep -v "#" 1999IZMITT01REIL.dat | proj +proj=utm +zone=36 -r | awk 'BEGIN{print "# n slip      x1      x2      x3 length width strike dip rake"}{print NR,$6,($2-4510043.99)/1e3,($1-245446.49)/1e3,$5,$9,$10,$7,$8,$11}' > reilinger+00_km.flt
