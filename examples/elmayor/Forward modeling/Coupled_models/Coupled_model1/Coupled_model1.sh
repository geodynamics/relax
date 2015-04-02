#!/bin/sh

# El Mayor-Cucapah earthquake: multiple-mechanism model combining:
# 1) rate-strengthening afterslip on downward extensions of planes F2 and F3 in Wei et al [2011] slip model, extending down approximately to Moho depth (~23-25 km depth), and on NW-striking vertical plane in Yuha Desert, extending from surface to 12 km depth; (a-b)sigma = 1 MPa
# 2) Newtonian viscoelastic relaxation in narrow ductile zone in lower crust aligned with locations of high heat flow and geothermal activity in Salton Trough; top depth 10 km; bottom depth 22.5 km (Moho depth)
# 3) Newtonian viscoelastic relaxation below 45 km depth in region of shallow Lekic et al [2011] lithosphere-asthenosphere boundary geometry in the Salton Trough and below 70 km depth outside the Salton Trough
# Note that the grid is 1024^3, which requires >50 GB of RAM on one node. If that exceeds capabilities, try setting sx3 to 512.

#PBS -N Coupled_model1
#PBS -l nodes=1:ppn=12
#PBS -l walltime=240:00:00
#PBS -V

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

FLT=../../../Wei_slipmodel/elmayor_wei_km_highslip.flt

WDIR=./Coupled_model1
OPTS=../../../GPS/gps_within200km.coord

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

SX=1024
DX=0.6
LEN=`echo $SX $DX | awk '{print $1*$2/3}'`

OMP_NUM_THREADS=12 time ../relax/build/relax --no-vtk-output --no-stress-output --no-proj-output $* <<EOF | tee $WDIR/in.param
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
3 3 1
# number of observation planes
0
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
0
# number of prestress interfaces with depth
0
# number of linear viscous interfaces
1
# no depth gammadot0 cohesion
   1    70       0.9   0.0
# number of linear ductile zones
7
# no dgammadot0 x1  x2 x3 length width thickness strike dip
   1       0.9 102.5 -28.16 45   400  25       125    147  90
   2      -2.5 119 -3 45   80    7.5     65    147  90
   3      -2.5 90  -70 45  70    25     45    180  90
   4      -2.5 110 40  45  30  25   100    0  90
   5      -2.5 55 130 45  400  25    130      147   90
   6       0.3  80 -15 10  110   12.5    30      155   90
   7       0.3 -18  30 10  400   12.5    30      147   90
# number of powerlaw viscous interfaces
0
# number of fault creep interface
1
# nb depth gammadot0 (a-b)sig friction cohesion
   1     0       3.0        1      0.6        0
# number of afterslip planes
3
# nb      x1       x2  x3 length width strike   dip rake
   1    -23.289   20.286   14.572     51  9    311.2    75 -150
   2    +10.941  -25.585   0.083     24  12    305    90  180
   3    -28.997   11.004   13.068    66   9    130      60 -150
# number of interseismic dislocation
0
# number of interseismic dykes
0
# number of coseismic events
1
# number of shear dislocations
`awk 'BEGIN{c=0}{if ($2 > 0.00){c=c+1}}END{print c}' $FLT`
# index slip x1 x2 x3 length width strike dip rake
`awk 'BEGIN{c=1}{if ($2 > 0.00){$1=c;print $0;c=c+1}}' $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface loads
0
EOF

