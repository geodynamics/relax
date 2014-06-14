#!/bin/bash

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

flt2xyz(){
	grep -v "#" $FLTFILE | awk '
	BEGIN{
		pi=atan2(1,0)*2;
	}
	{
		i=$1;
		if (NF==9){
		x1=$2;
		x2=$3;
		x3=$4;
		L=$5;
		W=$6;
		strike=$7*pi/180;
		dip=$8*pi/180;
		rake=$9*pi/180;
		nstrike=7;
		ndip=8;
		} else {
		slip=$2
		x1=$3;
		x2=$4;
		x3=$5;
		L=$6;
		W=$7;
		strike=$8*pi/180;
		dip=$9*pi/180;
		rake=$10*pi/180;
		nstrike=8;
		ndip=9;
		}
		s[0]=cos(strike);
		s[1]=sin(strike);
		s[2]=0;
		d[0]=+sin(strike)*cos(dip);
		d[1]=-cos(strike)*cos(dip);
		d[2]=-sin(dip);
		b[0]=s[0]*cos(rake)+d[0]*sin(rake);
		b[1]=s[1]*cos(rake)+d[1]*sin(rake);
		b[2]=d[2]*sin(rake);
		if (10==NF){
			printf("> -Z%f\n",slip);
		} else {
			printf(">\n");
		}
                printf("%f %f %f\n",x2              ,x1              ,-x3       );
                printf("%f %f %f\n",x2       +s[1]*L,x1       +s[0]*L,-x3       );
		printf("%f %f %f\n",x2-d[1]*W+s[1]*L,x1-d[0]*W+s[0]*L,-x3+d[2]*W);
                printf("%f %f %f\n",x2-d[1]*W       ,x1-d[0]*W       ,-x3+d[2]*W);
	}
	END{
	}' > $XYZFILE
}

usage(){
        echo "$self converts an .flt file to a .xyz polygon"
	echo "file for visualization in GMT with psxyz"
	echo ""
	echo "usage: $self file.flt"
	echo "       $self file1.flt file2.flt file3.flt "
	echo ""
	echo "or from the standard input"
	echo ""
	echo "       cat file.flt | $self > file.xyz"
	echo ""
	echo "file foo.flt is converted to foo.xyz"
	echo "files foo.ext are converted to foo.ext.xyz"
	echo ""
}

if [ -t 0 ] && [ $# -eq 0 ]; then
	usage
	exit
fi

if [ ! -t 0 ]; then
	FLTFILE=-
	XYZFILE=/dev/stdout
	flt2xyz
else
	# loop over list of files to convert
	while [ $# -ne 0 ];do
		# define output file name (file.flt is converted to file.xyz, file.ext to file.ext.xyz)
		XYZFILE=$(dirname $1)/$(basename $1 .flt).xyz
		FLTFILE=$1
		if [ ! -e $FLTFILE ]; then
			echo $self: could not find $FLTFILE. exiting.
			exit 2
		fi
		echo $self: converting $1 to $XYZFILE
		flt2xyz
		shift
	done
fi

