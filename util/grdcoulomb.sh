#!/bin/bash

set -e
self=$(basename $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

while getopts "d:f:hs:r:" flag; do
	case "$flag" in
	d) dset=1;dip=$OPTARG;;
	f) fset=1;mu=$OPTARG;;
	h) hset=1;;
	r) rset=1;rak=$OPTARG;;
	s) sset=1;str=$OPTARG;;
	esac
done
for item in $dset $fset $rset $sset;do
	shift;shift
done
for item in $hset;do
	shift
done

if [ $# -lt "1" -o "$hset" == "1" ]; then
	less -R <<EOF
[1m$self[0m
  projects a stress tensor produced by Relax into surface tractions and Coulomb stress
  along a hypothetic fault plane defined by strike and dip orientation and rake angle.

[1musage:[0m
  $self [-s strike] [-d dip] [-r rake] [-f friction] index
  or
  $self [-s strike] [-d dip] [-r rake] [-f friction] index1 ... indexN

  where index is the prefix of the stress components created by Relax, for example
  out/000 for out/000-s11.grd, out/000-s12.grd, out/000-s13.grd, out/000-s22.grd and 
  out/000-s33.grd.

[1moptions:[0m
  -s strike   set the strike angle in degrees of the receiver fault [0]
  -d dip      set the dip angle in degrees of the receiver fault [90]
  -r rake     set the rake in degrees of the receiver fault for Coulomb stress calculation [0]
  -f friction set the friction coefficient of the receiver fault for Coulomb stress calculation [0]

[1mexamples:[0m
  to compute the Coulomb stress for right-lateral motion on a north-striking, vertical, fault:

     $self -r 180 output/000

  to compute the Coulomb stress for a NNW gently dipping thrust fault:

     $self -s -45 -d 10 -r 90 output/000

[1moutput files:[0m
  create 7 files in the index base directory.
  index-t1.grd (north component of traction vector)
  index-t2.grd (east component of traction vector)
  index-t3.grd (depth component of traction vector)
  index-ts.grd (strike-slip component of traction vector, positive for left-lateral moment)
  index-td.grd (dip-slip component of traction vector, positive for thrust when 0<dip<90)
  index-tn.grd (normal component of traction vector, positive for extension)
  index-tc.grd (Coulomb stress)

[1mdescription:[0m
  using the coordinate system

             (x1, north)
            /
           /) strike
          +-----------------------+      (x2, east)
          |\        p .            \
          :-\      i .              \
          |  \    l .                \
          :90 \  s .                  \
          |-dip\  .                    \
          :     \. ) rake               \
          |      +-----------------------+
          :
          Z (x3)

  the fault normal, dip and strike vectors are:

     n = [-sin(str) sin(dip), +cos(str) sin(dip), -cos(dip)]
     s = [cos(str), sin(str), 0]
     d = [+sin(str) cos(dip), -cos(str) cos(dip), -sin(dip)]

  and form a right-handed coordinate system.

  the local traction components are:

     t1 = S11 n(1) + S12 n(2) + S13 n(3)
     t2 = S12 n(1) + S22 n(2) + S23 n(3)
     t3 = S13 n(1) + S23 n(2) + S33 n(3)

  and the strike, normal and dip components of the traction vector are

     tn = t(1) n(1) + t(2) n(2) + t(3) n(3)
     ts = t(1) s(1) + t(2) s(2)
     td = t(1) d(1) + t(2) d(2) + t(3) d(3)

  the Coulomb stress is

     tc = t(1) b(1) + t(2) b(2) + t(3) b(3) + mu * tn

  by convention, ts is positive for left lateral strike shear. dip direction
  is up dip. the 1, 2 and 3 directions are north, east and down, respectively.
  the vector s, d, and n form a right-handed coordinate system.
EOF

	exit
fi

if [ "$dset" != "1" ]; then
	dip=90
fi
echo $self: using dip angle $dip

if [ "$sset" != "1" ]; then
	str=0
fi
echo $self: using strike angle $str

if [ "$rset" != "1" ]; then
	rak=0
fi
echo $self: using rake angle $rak

if [ "$fset" != "1" ]; then
	mu=0
fi
echo $self: using friction coefficient $mu

# x1, x2 and x3 components of the normal, strike and dip unit vectors ( s(3) = 0 )
n1=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", -sin(s/180*pi)*sin(d/180*pi) }'`
n2=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", +cos(s/180*pi)*sin(d/180*pi) }'`
n3=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", -cos(d/180*pi) }'`
s1=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", +cos(s/180*pi) }'`
s2=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", +sin(s/180*pi) }'`
d1=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", +sin(s/180*pi)*cos(d/180*pi) }'`
d2=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", -cos(s/180*pi)*cos(d/180*pi) }'`
d3=`echo "" | awk -v d=$dip -v s=$str 'BEGIN{pi=atan2(1,0)*2}{printf "%f", -sin(d/180*pi) }'`

# slip (Burgers) vector components
b1=`echo $s1 $d1 | awk -v r=$rak 'BEGIN{pi=atan2(1,0)*2}{print cos(r/180*pi)*$1 + sin(r/180*pi)*$2}'`
b2=`echo $s2 $d2 | awk -v r=$rak 'BEGIN{pi=atan2(1,0)*2}{print cos(r/180*pi)*$1 + sin(r/180*pi)*$2}'`
b3=`echo   0 $d3 | awk -v r=$rak 'BEGIN{pi=atan2(1,0)*2}{print cos(r/180*pi)*$1 + sin(r/180*pi)*$2}'`

echo "$self: n=($n1,$n2,$n3), s=($s1,$s2,0), d=($d1,$d2,$d3), b=($b1,$b2,$b3)"

# loop over index files for batch processing
while [ "$#" != "0" ];do

	INDEX=$(basename $1)

	WDIR=$(dirname $1)
	echo $self": projecting index $INDEX in $WDIR in traction"

	# input stress components
	GRD11=$WDIR/$INDEX-s11.grd
	GRD12=$WDIR/$INDEX-s12.grd
	GRD13=$WDIR/$INDEX-s13.grd
	GRD22=$WDIR/$INDEX-s22.grd
	GRD23=$WDIR/$INDEX-s23.grd
	GRD33=$WDIR/$INDEX-s33.grd

	# output traction components in the x1, x2, and x3 directions
	T1=$WDIR/$INDEX-t1.grd
	T2=$WDIR/$INDEX-t2.grd
	T3=$WDIR/$INDEX-t3.grd

	# output traction components in the strike, dip and normal directions
	GRDTS=$WDIR/$INDEX-ts.grd
	GRDTD=$WDIR/$INDEX-td.grd
	GRDTN=$WDIR/$INDEX-tn.grd

	# output Coulomb stress
	GRDTC=$WDIR/$INDEX-tc.grd

	# traction components in the x1, x2, and x3 directions
	grdmath $GRD11 $n1 MUL $GRD12 $n2 MUL ADD $GRD13 $n3 MUL ADD = $T1
	grdmath $GRD12 $n1 MUL $GRD22 $n2 MUL ADD $GRD23 $n3 MUL ADD = $T2
	grdmath $GRD13 $n1 MUL $GRD23 $n2 MUL ADD $GRD33 $n3 MUL ADD = $T3

	# projection in fault-aligned reference system

	# normal component of traction vectors 
	# Tn=t(1) n(1) + t(2) n(2) + t(3) n(3)
	grdmath $T1 $n1 MUL $T2 $n2 MUL ADD $T3 $n3 MUL ADD = $GRDTN

	# strike component of traction ( s(3)=0 )
	# ts = t(1) s(1) + t(2) s(2)
	grdmath $T1 $s1 MUL $T2 $s2 MUL ADD = $GRDTS

	# dip component of traction
	# td = t(1) d(1) + t(2) d(2) + t(3) d(3)
	grdmath $T1 $d1 MUL $T2 $d2 MUL ADD $T3 $d3 MUL ADD = $GRDTD

	# Coulomb stress \tau+\mu t_n (t_n is positive for extension)
	# tc = t(1) b(1) + t(2) b(2) + t(3) b(3) + mu * tn
	grdmath $T1 $b1 MUL $T2 $b2 MUL ADD $T3 $b3 MUL ADD $GRDTN $mu MUL ADD = $GRDTC

	shift
	if [ "$#" != "0" ];then
		echo ""
	fi
done


