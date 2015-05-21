
# find center of fault patch
echo -118.52930 34.36385 | proj +proj=utm +zone=11 | awk '{print ($2-3828373)/1e3,($1-566940.91)/1e3}'

# fault left corner of fault patch
extrude.sh
# nb   x1   x2   x3 length width thickness strike dip
 1 -24.8134 -207.57 0 -10 0 0 109.6 0
-2.145888e+01 -2.169906e+02 +0.000000e+00

# subsample the fault
echo 1 -2.145888e+01 -2.169906e+02 1 20 26 109.6 41 90 2 2 1 1 | seg2flt.py > geometry

# line up strike slip
cat <<EOF | xargs | awk '{OFS="\n";for (i=1;i<=NF;i++){print -$i}}' > temp-strike-slip
cat <<EOF | xargs | awk '{OFS="\n";for (i=1;i<=NF;i++){print +$i}}' > temp-dip-slip

# create .flt file and convert to m
paste temp-strike-slip temp-dip-slip | awk '{pi=atan2(1,0)*2;print sqrt($1^2+$2^2),atan2($2,$1)/pi*180}' | paste - geometry | awk 'BEGIN{print "# n x1 x2 x3 length width strike dip rake"}{print $3,$1/100,$4,$5,$6,$7,$8,$9,41,$2}' > hudnut+96.flt
