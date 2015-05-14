
# reference point
echo 13.334 42.334 | proj +proj=utm +zone=33
362747.88       4688204.99

# convert the Atzori et al. (2009) model for the 2009 Mw 6.3 L'Aquila earthquake to the Relax format
grep -v "#" grl26191-sup-0002-ds01.txt | awk 'BEGIN{print "# n    slip         x1         x2         x3   length    width     strike       dip        rake"}{printf "%03.3d %7.4f %10.3e %10.3e %10.3e %f %f %f %f %f\n",NR,$9/1e2,($7-4688204.99)/1e3,($6-362747.88)/1e3,$3/1e3,$1/1e3,$2/1e3,$4,$5,$8}' > atzori+09_km.flt
grep -v "#" grl26191-sup-0002-ds01.txt | awk 'BEGIN{print "# n    slip            x1            x2            x3      length       width     strike       dip        rake"}{printf "%03.3d %7.4f %13.6e %13.6e %13.6e %f %f %f %f %f\n",NR,$9/1e2,($7-4688204.99),($6-362747.88),$3,$1,$2,$4,$5,$8}' > atzori+09.flt
