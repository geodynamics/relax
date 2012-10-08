#!/bin/bash

set -e
self=$(basename $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if [ "$#" == "0" ]; then
	echo usage $self [-r gpsr.txt] gps1.txt [... gpsN.txt]
	echo  
	echo removes the coseismic component from the time series gps1.txt ... gpsN.txt
	echo and optionally outputs the displacement relative to station gpsr
	echo 
	echo example:  obsrelax.sh "output{1,2,3,4,5,6,7}/{lae1,lae2,lae3,lae4,law1,law2,law3,law4,oldw,oldd,meek,rich,sanh}"
	exit 1
fi

while getopts "r:" flag
do
  case "$flag" in
    r) rset=1;RFILE=$OPTARG;;
  esac
done
for i in $rset;do
        shift;shift
done

if [ "$rset" == "1" ];then
	TMP=temp
	REF=$(basename $RFILE .txt)
	grep -v "#" $RFILE | \
		awk '{if (1==NR){i2=$2;i3=$3;i4=$4;i5=$5;i6=$6;i7=$7;i8=$8;i9=$9;i10=$10}; \
			print $2-i2,$3-i3,$4-i4,$5-i5,$6-i6,$7-i7,$8-i8,$9-i9,$10-i10}' > $TMP
fi

while [ "$#" != "0" ]; do
	WDIR=$(dirname $1)
	IFILE=$WDIR/$(basename $1 .txt).txt
	if [ "$rset" != "1" ];then
		# postseismic component
		OFILE=$WDIR/$(basename $1 .txt)-relax.txt
		echo $self: changing $IFILE to $OFILE
		grep -v "#" $IFILE | 
			awk 'BEGIN{print "#         t         u1         u2         u3        s11        s12        s13        s22        s23        s33"}{if (1==NR){i2=$2;i3=$3;i4=$4;i5=$5;i6=$6;i7=$7;i8=$8;i9=$9;i10=$10}; \
			$2=$2-i2;$3=$3-i3;$4=$4-i4;$5=$5-i5;$6=$6-i6;$7=$7-i7;$8=$8-i8;$9=$9-i9;$10=$10-i10;print $0}' > $OFILE
	else
		# postseismic displacement relative to reference station
		OFILE=$WDIR/$(basename $1 .txt)-$REF-relax.txt
		echo $self: change $IFILE to $OFILE
		grep -v "#" $IFILE | awk '{if (1==NR){i2=$2;i3=$3;i4=$4;i5=$5;i6=$6;i7=$7;i8=$8;i9=$9;i10=$10}; \
			$2=$2-i2;$3=$3-i3;$4=$4-i4;$5=$5-i5;$6=$6-i6;$7=$7-i7;$8=$8-i8;$9=$9-i9;$10=$10ii10;print $0}' | \
			paste - $TMP | \
			awk 'BEGIN{print "#         t         u1         u2         u3        s11        s12        s13        s22        s23        s33"}{$2=$2-$11;$3=$3-$12;$4=$4-$13;$5=$5-$14;$6=$6-$15;$7=$7-$16;$8=$8-$18;$9=$9-$18;$10=$10-$19; \
				$11="";$12="";$13="";print $0}' > $OFILE
	fi
	shift
done

if [ "$rset" == "1" ];then
	rm -f $TMP
fi
