#!/bin/sh

# El Mayor-Cucapah earthquake: Newtonian viscoelastic relaxation in mantle lithosphere geometry between Tape et al [2012] Moho surface and approximate Lekic et al [2011] lithosphere-asthenosphere boundary surface
# Note that the grid is 1024^3, which requires >50 GB of RAM on one node. If that exceeds capabilities, try setting sx3 to 512.

#PBS -N Mantle_lithosphere_Newtonian_rheology_topdepth_Mohodepth
#PBS -q default
#PBS -l nodes=12
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

WDIR=./Mantle_lithosphere_Newtonian_rheology_topdepth_Mohodepth
OPTS=../../../GPS/gps_within200km.dat

if [ ! -e $WDIR ]; then
	echo adding directory $WDIR
	mkdir $WDIR
fi

SX=1024
DX=0.6
LEN=`echo $SX $DX | awk '{print ($1-1)*$2/2}'`

OMP_NUM_THREADS=12 time relax --no-vtk-output --no-stress-output --no-proj-output $* <<EOF | tee $WDIR/in.param
# grid size (sx1,sx2,sx3)
$SX $SX 512
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
3 -10 1
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
   1    70         0   0.0
# number of linear ductile zones
408
# no dgammadot0 x1  x2 x3 length width thickness strike dip
`grep -v "#" ../../../tape_mantlelith/mantlelith_moho_speed2.dat`
 401      -2.0 51.91 40.57 52.5   50   17.5     65    327  90
 402      -2.0 38.29 19.60 52.5   80   17.5     15    327  90
 403      -2.0 110 12.5  52.5  29  17.5   60    180  90
 404      -2.0 90  -35 45  75    25     25    180  90
 405      -2.0 35.41 15.41 45   400  25        125    147  90
 406      -2.0 28.76 4.92 45   80    25     20    327  90
 407      -2.0 12.42 -20.24 45   20    25     40    327  90
 408      -2.0 17.87 -11.85 45   40    25     20    327  90
# number of powerlaw viscous interfaces
0
# number of fault creep interface
0
# number of interseismic dislocation
0
# number of interseismic dykes
0
# number of coseismic events
1
# number of shear dislocations
`awk 'BEGIN{c=0}{if ($2 >= 0.00){c=c+1}}END{print c}' $FLT`
# index slip x1 x2 x3 length width strike dip rake
`awk 'BEGIN{c=1}{if ($2 >= 0.00){$1=c;print $0;c=c+1}}' $FLT`
# number of tensile cracks
0
# number of dilatation sources
0
# number of surface loads
0
EOF

# create postseismic time series with coseismic offset removed
obsrelax.sh $WDIR/????.txt | tee -a $WDIR/in.param

