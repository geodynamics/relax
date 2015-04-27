#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# script cpt2rgb.sh
# converts a z-value to red-green-blue using a GMT color scale.
#
# Sylvain Barbot (04/27/15) - original form
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set -e
self=$(basename $0)
selfdir=$(dirname $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

usage(){
        echo "usage: $self -C cpt"
	echo ""
	echo "converts a z-value to hexadecimal rgb using a GMT color palette."
	echo ""
        echo "options:"
        echo "         -C cpt "
        echo "         -h display this message"

	exit
}

while getopts "C:h" flag
do
	case "$flag" in
	C) Cset=1;CPT=$OPTARG;;
	h) hset=1;;
	esac
done
for item in $Cset ; do
	shift;shift
done
for item in $hset ; do
	shift
done

if [ "1" == "$hset" ]; then
	usage
fi

if [ "1" != "$Cset" ]; then
	CPT=$selfdir/../share/jet.cpt
fi

cat - | while read -r VALUE; do
	grep -v "#" $CPT | awk -v value=$VALUE ' \
		{
			if ($1<=value && value<=$5){
				c=(value-$1)/($5-$1);
				r=int($2+c*($6-$2));
				g=int($3+c*($7-$3));
				b=int($4+c*($8-$4));
			};
		} \
		END{print r,g,b}'
done




