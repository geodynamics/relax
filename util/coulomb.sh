#!/bin/bash -e

usage(){
        echo "usage: $self "
	echo "or"
	echo "       cat file.dat | $self"
	echo ""
	echo " compute normal, shear, and Coulomb stress for an arbitrarily oriented"
	echo " receiver fault based on the stress tensor components s11,s12,s13,s22,s23,s33 "
	echo "found in the index wdir/index."
	echo ""
	echo "usage: $self"
	echo ""
        echo "options:"
        echo "         -h display this message"
	echo "         -m reference friction coefficient of receiver fault"
	echo "         -s strike angle of receiver fault"
	echo "         -d dip angle of receiver fault"
	echo "         -r rake angle of receiver fault"
	echo ""
        echo "description:"
	echo ""
	echo "                   North (x1)"
	echo "                  / "
	echo "                 /) Strike"
	echo "     x1,x2,x3 ->@------------------------  (x2)"
	echo "                :\                       \"
	echo "                :-\                       \"
	echo "                :  \                       \"
	echo "                :90 \                       \"
	echo "                :-Dip\                       \"
	echo "                :     \                       \"
	echo "                :      -------------------------"
	echo "                :        "
	echo "                Z (x3)"
	echo ""
	echo ""
	
        exit
}

while getopts "d:hm:r:s:" flag
do
  case "$flag" in
    d) dset=1;dip=$(echo $OPTARG | awk '{print $1*3.1415926535897932385/180.0}');;
    h) hset=1;;
    m) mset=1;mu0=$OPTARG;;
    r) rset=1;rak=$(echo $OPTARG | awk '{print $1*3.1415926535897932385/180.0}');;
    s) sset=1;str=$(echo $OPTARG | awk '{print $1*3.1415926535897932385/180.0}');;
  esac
done
shift $((OPTIND -1))

# usage
if [ "1" == "$hset" ]; then
	usage
fi

# default values
if [ "" == "$mset" ]; then
	mu0=0.6
fi
if [ "" == "$rset" ]; then
	rak=0
fi
if [ "" == "$sset" ]; then
	str=0
fi
if [ "" == "$dset" ]; then
	dip=1.570796326794897
fi

BC=$(tput bold)
NC=$(tput sgr0)
RC="\033[0;31m"

while [ "$#" != "0" ];do
	WDIR=$(dirname $1)
	index=$(basename $1)
	echo -e "# processing Coulomb stress in ${RC}${BC}$WDIR/${index}${NC}"

	if [ ! -e "$WDIR/${index}-s12.grd" ]; then
		echo -e "# ${RC}error: stress components ${index}-s12.grd does not exit in $WDIR/"
		exit -1
	fi

	# strike vector
	s1=$(echo $str | awk 'BEGIN{pi=atan2(0,1)*2}{print cos($1)}')
	s2=$(echo $str | awk 'BEGIN{pi=atan2(0,1)*2}{print sin($1)}')
	s3=0

	# dip vector (positive for thrust motion)
	d1=$(echo $str $dip | awk 'BEGIN{pi=atan2(0,1)*2}{print +sin($1)*cos($2)}')
	d2=$(echo $str $dip | awk 'BEGIN{pi=atan2(0,1)*2}{print -cos($1)*cos($2)}')
	d3=$(echo $str $dip | awk 'BEGIN{pi=atan2(0,1)*2}{print -sin($2)}')

	# normal vector
	n1=$(echo $str $dip | awk 'BEGIN{pi=atan2(0,1)*2}{print -sin($1)*sin($2)}')
	n2=$(echo $str $dip | awk 'BEGIN{pi=atan2(0,1)*2}{print +cos($1)*sin($2)}')
	n3=$(echo $str $dip | awk 'BEGIN{pi=atan2(0,1)*2}{print -cos($2)}')

	# rake vector
	r1=$(echo $rak $s1 $d1 | awk 'BEGIN{pi=atan2(0,1)*2}{print cos($1)*$2+sin($1)*$3}')
	r2=$(echo $rak $s2 $d2 | awk 'BEGIN{pi=atan2(0,1)*2}{print cos($1)*$2+sin($1)*$3}')
	r3=$(echo $rak $s3 $d3 | awk 'BEGIN{pi=atan2(0,1)*2}{print cos($1)*$2+sin($1)*$3}')

	# traction vector
	T1="$WDIR/$index-receiver-t1.grd"
	T2="$WDIR/$index-receiver-t2.grd"
	T3="$WDIR/$index-receiver-t3.grd"

	# compute traction vector
	gmt grdmath $WDIR/${index}-s11.grd $n1 MUL $WDIR/${index}-s12.grd $n2 MUL ADD $WDIR/${index}-s13.grd $n3 MUL ADD = $T1
	gmt grdmath $WDIR/${index}-s12.grd $n1 MUL $WDIR/${index}-s22.grd $n2 MUL ADD $WDIR/${index}-s23.grd $n3 MUL ADD = $T2
	gmt grdmath $WDIR/${index}-s13.grd $n1 MUL $WDIR/${index}-s23.grd $n2 MUL ADD $WDIR/${index}-s33.grd $n3 MUL ADD = $T3
	echo "# ${BC}traction components${NC}: $(basename $T1) $(basename $T2) $(basename $T3)"

	# normal component
	TN="$WDIR/$index-receiver-tn.grd"

	# normal component (scalar, positive for extension)
	gmt grdmath $T1 $n1 MUL $T2 $n2 MUL ADD $T3 $n3 MUL ADD = $TN
	echo "# ${BC}normal stress${NC}: $(basename $TN) (positive for extension)"

	# normal component
	TS="$WDIR/$index-receiver-ts.grd"

	# shear component in the rake direction (scalar, positive in the rake direction)
	gmt grdmath $T1 $r1 MUL $T2 $r2 MUL ADD $T3 $r3 MUL ADD = $TS
	echo "# ${BC}shear component${NC}: $(basename $TS) (in rake direction)"

	# Coulomb stress
	CS=$WDIR/${index}-receiver-coulomb.grd

	# Coulomb stress
	gmt grdmath "$TN" "$mu0" MUL "$TS" ADD = "$CS"
	echo "# ${BC}Coulomb stress${NC}: $(basename $CS)"

	shift
	if [ "$#" -gt 0 ]; then
		echo ""
	fi
done


