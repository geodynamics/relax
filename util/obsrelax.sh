#!/bin/bash

set -e
self=$(basename $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

usage(){
	echo "usage: $self [-r gpsr.txt] [-t time0/../time1] gps1.txt [... gpsN.txt]"
	echo or
	echo "       $self -t 0/7 gp*.txt"
	echo  
	echo removes the coseismic component from the time series gps1.txt ... gpsN.txt
	echo and optionally outputs the displacement relative to station gpsr
	echo ""
	echo options:
	echo    -r  gps1.txt  compute the displacements relative to station gps1.txt
	echo    -t  time      remove coseismic displacements at given times t1/../t3 [0].
	echo 
	echo example:  obsrelax.sh "output{1,2,3,4,5,6,7}/{lae1,lae2,lae3,lae4,law1,law2,law3,law4,oldw,oldd,meek,rich,sanh}"
	echo 
}

stripoffset(){
	WDIR=$(dirname $1)
	IFILE=$WDIR/$(basename $1 .txt).txt
	if [ "$rset" != "1" ];then
		EVENTS=`grep -v "#" $IFILE | \
			awk -v t="$EPOCHS" 'BEGIN{ \
			split(t,time,"/"); \
			}{ \
			for (t in time){ \
				if (time[t]<=$1){ \
					printf "%d ",NR; \
					delete time[t]; \
				}; \
			}; \
			if (0==length(time)){ \
				exit \
			}; \
			}'`

		# postseismic component
		if [ ! -t 1 ]; then
			# if output (1) is redirected, print to stdout
			OFILE=/dev/stdout
		else
			OFILE=$WDIR/$(basename $1 .txt)-relax.txt
			echo $self: changing $IFILE to $OFILE
		fi

		grep -v "#" $IFILE |
			awk -v e="$EVENTS" \
				'BEGIN{ \
					print "#         t         u1         u2         u3        s11        s12        s13        s22        s23        s33"; \
					split(e,events," "); \
					p2=0;p3=0;p4=0;p5=0;p6=0;p7=0;p8=0;p9=0;p10=0; \
					o2=0;o3=0;o4=0;o5=0;o6=0;o7=0;o8=0;o9=0;o10=0; \
				}{ \
				for (e in events){ \
					if (events[e]==NR){ \
						o2=o2+$2-p2;o3=o3+$3-p3;o4=o4+$4-p4;o5=o5+$5-p5; \
						o6=o6+$6-p6;o7=o7+$7-p7;o8=o8+$8-p8;o9=o9+$9-p9;o10=o10+$10-p10; \
					} \
				}; \
				p2=$2;p3=$3;p4=$4;p5=$5;p6=$6;p7=$7;p8=$8;p9=$9;p10=$10;
				$2=$2-o2;$3=$3-o3;$4=$4-o4;$5=$5-o5;$6=$6-o6;$7=$7-o7;$8=$8-o8;$9=$9-o9;$10=$10-o10; \
				printf "%11.4e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n", \
					$1,$2,$3,$4,$5,$6,$7,$8,$9,$10; \
			}' > $OFILE
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
}

if [ "$#" == "0" ]; then
	usage
	exit 1
fi

while getopts "r:t:" flag
do
  case "$flag" in
    r) rset=1;RFILE=$OPTARG;;
    t) tset=1;EPOCHS="$OPTARG";;
  esac
done
for i in $rset $tset; do
        shift;shift
done

if [ "$rset" == "1" ];then
	TMP=temp
	REF=$(basename $RFILE .txt)
	grep -v "#" $RFILE | \
		awk '{if (1==NR){i2=$2;i3=$3;i4=$4;i5=$5;i6=$6;i7=$7;i8=$8;i9=$9;i10=$10}; \
			print $2-i2,$3-i3,$4-i4,$5-i5,$6-i6,$7-i7,$8-i8,$9-i9,$10-i10}' > $TMP
fi

if [ "$tset" != "1" ];then
	EPOCHS="0"
fi

while [ "$#" != "0" ]; do
	stripoffset $1
	shift
done

if [ "$rset" == "1" ];then
	rm -f $TMP
fi

