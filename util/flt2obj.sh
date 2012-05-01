#!/bin/sh


#  This is a cube with negative vertex reference numbers. Each element
#  references the vertices stored immediately above it in the file. Note
#  that vertices are not shared.
#
#  v 0.000000 2.000000 2.000000
#  v 0.000000 0.000000 2.000000
#  v 2.000000 0.000000 2.000000
#  v 2.000000 2.000000 2.000000
#  f -4 -3 -2 -1

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if [ $# -eq 0 ]; then
        echo "$self converts an .flt file to a obj polygon"
	echo "file in ASCII format for visualization in LightWave"
	echo ""
	echo "usage: $self [-c com] file.flt"
	echo "       $self file1.flt file2.flt file3.flt "
	echo ""
	echo "file foo.flt is converted to foo.obj"
	echo "files foo.ext are converted to foo.ext.obj"
	echo ""

	exit
fi

# optional command-line parameters
while getopts "c:" flag
do
  case "$flag" in
    c) cset=1;com=$OPTARG;;
  esac
done
for item in $cset ;do
	shift;shift
done

if [ "$cset" != "1" ]; then
	com=">"
else
	echo $self: using separator $com
fi

# loop over list of files to convert
while [ $# -ne 0 ];do

FLTFILE=$1

if [ ! -e $FLTFILE ]; then
	echo $self: could not find $FLTFILE. exiting.
	exit 2
fi

# define output file name (file.flt is converted to file.obj, file.ext to file.ext.obj)
OBJFILE=$(dirname $1)/$(basename $1 .flt).obj

echo $self: converting $1 to $OBJFILE

grep -v "#" $FLTFILE | awk '
	BEGIN{
		pi=atan2(1,0)*2;
		printf("# flt2obj.sh %s\n",$FLTFILE);
	}
	{
		i=$1;
		x1=$2;
		x2=$3;
		x3=$4;
		L=$5;
		W=$6;
		strike=$7*pi/180;
		dip=$8*pi/180;
		s[0]=cos(strike);
		s[1]=sin(strike);
		s[2]=0;
		d[0]=+sin(strike)*cos(dip);
		d[1]=-cos(strike)*cos(dip);
		d[2]=-sin(dip);
                printf("v %f %f %f\n",x1              ,x2              ,x3              );
                printf("v %f %f %f\n",x1       +s[0]*L,x2       +s[1]*L,x3       +s[2]*L);
		printf("v %f %f %f\n",x1-d[0]*W+s[0]*L,x2-d[1]*W+s[1]*L,x3-d[2]*W+s[2]*L);
                printf("v %f %f %f\n",x1-d[0]*W       ,x2-d[1]*W       ,x3-d[2]*W       );
	}
	END{
		printf("\n");
	}' > $OBJFILE

	grep -v "#" $FLTFILE | awk 'BEGIN{o=1}{printf("f %d %d %d %d\n",o+0,o+1,o+2,o+3);o=o+4}' >> $OBJFILE

	shift
done
