#!/bin/bash

set -e
self=$(basename $0)
selfdir=$(dirname $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

tri2xyz(){
	grep -v "#" $NEDFILE | \
		awk 'NR==FNR {for (j=2;j<=4;j++){a[j-1,0+$1]=$j}; next} 
		{
		if ("#" != substr($0,0,1)) {
			if (NF==6){
			        print "> -Z ",$2
				print a[2,0+$3],a[1,0+$3],a[3,0+$3];
				print a[2,0+$4],a[1,0+$4],a[3,0+$4];
				print a[2,0+$5],a[1,0+$5],a[3,0+$5];
				print a[2,0+$3],a[1,0+$3],a[3,0+$3];
			} else {
			        print ">" 
				print a[2,0+$2],a[1,0+$2],a[3,0+$2];
				print a[2,0+$3],a[1,0+$3],a[3,0+$3];
				print a[2,0+$4],a[1,0+$4],a[3,0+$4];
				print a[2,0+$2],a[1,0+$2],a[3,0+$2];
			}
		}
	} END {
	}' - $TRIFILE > $XYZFILE

	if [ "1" == "$sset" ]; then
		grep -v "#" $NEDFILE | \
		awk -v s=$segment 'function max(x,y){return (x>=y)?x:y} function abs(x){return (x>=0)?x:-x}; NR==FNR {for (j=2;j<=4;j++){a[j-1,0+$1]=$j}; next} 
			{
			if ("#" != substr($0,0,1)) {
				if (NF==6){
					l1=sqrt((a[1,0+$3]-a[1,0+$4])^2+(a[2,0+$3]-a[2,0+$4])^2+(a[3,0+$3]-a[3,0+$4])^2);
					l2=sqrt((a[1,0+$3]-a[1,0+$5])^2+(a[2,0+$3]-a[2,0+$5])^2+(a[3,0+$3]-a[3,0+$5])^2);
					l3=sqrt((a[1,0+$5]-a[1,0+$4])^2+(a[2,0+$5]-a[2,0+$4])^2+(a[3,0+$5]-a[3,0+$4])^2);
				} else {
					l1=sqrt((a[1,0+$2]-a[1,0+$3])^2+(a[2,0+$2]-a[2,0+$3])^2+(a[3,0+$2]-a[3,0+$3])^2);
					l2=sqrt((a[1,0+$2]-a[1,0+$4])^2+(a[2,0+$2]-a[2,0+$4])^2+(a[3,0+$2]-a[3,0+$4])^2);
					l3=sqrt((a[1,0+$4]-a[1,0+$3])^2+(a[2,0+$4]-a[2,0+$3])^2+(a[3,0+$4]-a[3,0+$3])^2);
				}
				l=max(max(l1,l2),l3);
				if (l<s){
					print $0;
				} else {
					print "# ",$0
				}
			}
		}' - $TRIFILE 
fi

}

usage(){
        echo "$self converts an .tri file to a .xyz file"
	echo "compatible with GMT"
	echo ""
	echo "usage: $self file.tri"
	echo "       $self file1.tri file2.tri file3.tri "
	echo ""

	exit
}

if [ $# -eq 0 ]; then
	usage
fi

while getopts "hs:" flag
do
	case "$flag" in
	h) hset=1;;
	s) sset=1;segment=$OPTARG;;
	esac
done
for item in $sset; do
	shift; shift;
done
for item in $hset; do
	shift;
done

# loop over list of files to convert
while [ $# -ne 0 ]; do
	# define output file name (file.flt is converted to file.xyz, file.ext to file.ext.xyz)
	XYZFILE=$(dirname $1)/$(basename $1 .tri).xyz
	TRIFILE=$(dirname $1)/$(basename $1 .tri).tri
	NEDFILE=$(dirname $1)/$(basename $1 .tri).ned
	if [ ! -e $TRIFILE ]; then
		echo "# $self: could not find $TRIFILE. exiting."
		exit 2
	fi
	if [ ! -e $NEDFILE ]; then
		echo "# $self: could not find $NEDFILE. exiting."
		exit 2
	fi

	echo "# $self: converting $1 to $XYZFILE"
	echo "# $self $cmdline"
	tri2xyz
	shift
done

