#!/bin/bash

set -e
self=$(basename $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if (test $# -lt "1"); then
        echo "usage: $self [-s 45] [-d 90] [...] file1.txt ... fileN.txt"
	echo ""
        echo "options:"
        echo "         -d phi sets the dip angle to phi [default is 90]"
        echo "         -l L sets the lengths of slip patches" 
        echo "         -m m multiplies the slip by a factor m"
	echo "         -s theta sets the strike angle to theta [default is 0]"
        echo "         -w W sets the width of slip patches"
	echo ""
	echo "Creates source input files from slip distributions files"
	echo "exported in txt format from an observation plane. Files"
	echo "dir1/file1.flt ... dirN/indexN.flt may be used as input"
	echo "slip distribution"
	echo ""
	echo "Export of observation planes provides the amplitude of slip"
	echo "at a center location (x0,y0,z0). The corner location (x1,y1,z0)"
	echo "is given by"
	echo ""
	echo "   x0 = x1 + L/2 cos(theta) - W/2 sin(theta)*cos(phi)"
	echo "   y0 = y1 + L/2 sin(theta) + W/2 cos(theta)*cos(phi)"
	echo "   z0 = z1 + W/2 sin(phi)"
	echo ""
	echo "The first coordinates in files index.s00001.slip.txt"
	echo "contain (x0,y0,z0). Then follows along-strike and down-dip"
	echo "distance (x',z'). The last 3 columns are total slip, "
	echo "strike slip and dip slip."
	echo ""
        
        exit
fi

pi=3.14159265359

while getopts "d:l:m:s:w:" flag
do
  case "$flag" in
    d) dset=1;dip=$OPTARG;;
    l) lset=1;length=$OPTARG;;
    m) mset=1;m=$OPTARG;;
    s) sset=1;strike=$OPTARG;;
    w) wset=1;width=$OPTARG;;
  esac
done
for item in $dset $lset $mset $rset $sset $wset;do
	shift;shift
done

while [ "$#" != "0" ];do

WDIR=$(dirname $1)
IN=$WDIR/$(basename $1 .txt).txt
OUT=$WDIR/$(basename $1 .txt).flt

# default value for dip angle phi: 90
if [ "$dset" != "1" ]; then
	dip=90;
fi
# default value for strike angle theta: 0
if [ "$sset" != "1" ]; then
	strike=0;
fi

# default value for strike angle theta: 0
dlength=`head -n 2 $IN | awk '{if (NR==1){p=$4}else{print $4-p}}'`
if [ "$lset" != "1" ]; then
	length=$dlength;
fi
# default width is length
if [ "$wset" != "1" ]; then
	width=$length;
fi
# default scaling factor is 1
if [ "$mset" != "1" ]; then
	m=1;
fi

echo ${self}: using length=${length}, width=${width}, strike=${strike}, dip=${dip} 

cs=`echo $strike | awk '{print cos($1*3.14159265359/180)}'`
ss=`echo $strike | awk '{print sin($1*3.14159265359/180)}'`
cd=`echo $dip | awk '{print cos($1*3.14159265359/180)}'`
sd=`echo $dip | awk '{print sin($1*3.14159265359/180)}'`

awk -v s=$strike -v d=$dip \
    -v cs=$cs -v ss=$ss -v cd=$cd -v sd=$sd \
    -v m=$m \
    -v l=$length -v w=$width \
	'BEGIN{i=1}{x1=$1-l/2*cs+w/2*ss*cd;x2=$2-l/2*ss-w/2*cs*cd;x3=$3-l/2*sd;x1=x1+1e4*l;x2=x2+1e4*l;x3=x3+1e4*l; \
          if ($6^2>0){print i,-$6*m,x1-1e4*l,x2-l*1e4,x3-1e4*l,l,w,s,d,atan2($8,$7)*180/3.14159265359;i=i+1}}' \
 	$IN > $OUT

echo ${self}: $OUT created successfully

shift
if [ "$#" != "0" ];then
	echo ""
fi
done

