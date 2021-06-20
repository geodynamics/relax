#!/bin/bash

# best-fitting model for the 2010 Mw7.2 El Mayor earthquake
# from Rollins, Barbot & Avouac "Mechanisms of Postseismic 
# Deformation Following the 2010 El Mayor-Cucapah Earthquake",
# Pageoph, submitted, 2014

# El Mayor earthquake 
# Newtonian viscosity
# weak upper mantle (70 km -infinity)

#PBS -N linear_asthenosphere_below75_1024
#PBS -l nodes=12
#PBS -l walltime=80:00:00
#PBS -V

if [ "$PBS_O_WORKDIR" != "" ]; then
	echo "MPI Used:" `which mpirun`

	# change the working directory (default is home directory)
	echo working directory is $PBS_O_WORKDIR
	cd $PBS_O_WORKDIR

	# write out some information on the job
	echo running on host `hostname`
	echo time is `date`

	# define number of processors
	NPROCS=`wc -l < $PBS_NODEFILE`
	echo this job has allocated $NPROCS cpus

	# tell which nodes it is run on
	echo " "
	echo this jobs runs on the following processors:
	echo `cat $PBS_NODEFILE`
	echo " "
fi

FLT=faults/elmayor_wei_km_highslip3.flt

WDIR=linear_asthenosphere_below75_1024
OPTS=gps/baja_km.coord

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

SX=1024
DX=0.6
LEN=`echo $SX $DX | awk '{print $1*$2/3}'`

time relax --no-vtk-output --no-stress-output $* <<EOF | tee $WDIR/in.param
# grid size (sx1,sx2,sx3)
$SX $SX $SX
# dx1  dx2  dx3 beta nyquist
  $DX  $DX  $DX 0.25     1.0
# origin position and rotation
0 0 0
# geographic origin and UTM zone
#-115.5 32.5 11 1e3
# observation depth (displacement and stress)
0 7
# output directory
$WDIR
# elastic parameters and gamma = (1-nu) rho g / mu = 8.33e-7 /m = 8.33e-4 /km
3e1 3e1 8.33e-4
# integration time (t1), time steps and scale
3 -1 0.5
# number of observation planes
2
# nb  x1   x2 x3 length width strike dip
   1 166.6 -107.8 0     250   100    140  90
   2 -55.5 -23.42  0    180   100     40  90
# number of observation points
`grep -v "#" $OPTS | awk -v l=$LEN \
	'function abs(x){return (0>x)?-x:x}
	 BEGIN{i=0}
	 {if (abs($2)<l && abs($3)<l){i=i+1;}}
	 END{print i}'`
# no NAME x1 x2 x3
`grep -v "#" $OPTS | awk -v l=$LEN \
	'function abs(x){return (0>x)?-x:x}
	 BEGIN{i=0}
	 {if (abs($2)<l && abs($3)<l){i=i+1;print i,$1,$3,$2,0}}'`
# number of stress observation segments
4
# nb  x1   x2 x3 length width strike dip friction
   1  95   -26 0      5    10    315  90 0.4
   2  62   -50 0      5    10    315  90 0.4
   3  62  -100 0      5    10    310  90 0.4
   4  35     0 0      5    10    145  90 0.4
# number of prestress interfaces with depth
0
# number of linear viscous interfaces
1
# no depth gammadot0 cohesion
   1    70       0.9   0.0
# number of linear ductile zones
7
# no dgammadot0    x1     x2 x3 length width thickness strike dip
   1        0.9 102.5 -28.16 45    400    25       125    147  90
   2       -2.5 119.0  -3.00 45     80   7.5        65    147  90
   3       -2.5  90.0 -70.00 45     70    25        45    180  90
   4       -2.5 110.0  40.00 45     30    25       100      0  90
   5       -2.5  55.0 130.00 45    400    25       130    147  90
   6        0.4  80.0 -15.00 10    110  12.5        30    155  90
   7        0.4 -18.0  30.00 10    400  12.5        30    147  90
# number of powerlaw viscous interfaces
0
# number of fault creep interface
1
# nb depth gammadot0 (a-b)sig friction cohesion
   1     0       3.0        1      0.6        0
# number of afterslip planes
3
# nb      x1      x2     x3 length width strike dip rake
   1 -22.288  19.168 14.572     51     8  311.2  75 -150
   2  10.654 -19.751 11.674     40    13  305.0  75 -150
   3 -29.980  12.137 13.068     66     8  130.0  60 -150
# number of interseismic dislocation
0
# number of interseismic dykes
0
# number of coseismic events
1
# number of shear dislocations
`awk 'BEGIN{c=0}{if ($5 >= 0 && $5 < 25  && $2>0){c=c+1}}END{print c}' $FLT`
# index slip x1 x2 x3 length width strike dip rake
`awk 'BEGIN{c=1}{if ($5 >= 0 && $5 < 25  && $2>0){$1=c;print $0;c=c+1}}' $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface loads
0
EOF

