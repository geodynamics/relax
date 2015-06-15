#!/bin/bash

set -e
self=$(basename $0)
selfdir=$(dirname $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

tri2vtp(){
	cat <<EOF > $VTPFILE
<?xml version="1.0"?>
<VTKFile type="PolyData" version="0.1">
<PolyData>
EOF

	grep -v "#" $NEDFILE | \
		awk 'NR==FNR {for (j=2;j<=4;j++){a[j-1,0+$1]=$j}; next} 
		{
		if ("#" != substr($0,0,1)) {
			printf("\t<Piece NumberOfPoints=\"3\" NumberOfPolys=\"1\">\n");
			printf("\t\t<Points>\n");
			printf("\t\t\t<DataArray type=\"Float32\" Name=\"Fault Patch\" NumberOfComponents=\"3\" format=\"ascii\">\n");
			print a[1,0+$3],a[2,0+$3],a[3,0+$3];
			print a[1,0+$4],a[2,0+$4],a[3,0+$4];
			print a[1,0+$5],a[2,0+$5],a[3,0+$5];
		}
			printf("\t\t\t</DataArray>\n");
			printf("\t\t</Points>\n");
			printf("\t\t<Polys>\n");
			printf("\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"2\">\n");
			printf("0 1 2\n");
			printf("\t\t\t</DataArray>\n");
			printf("\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"3\" RangeMax=\"3\">\n");
			printf("3\n");
			printf("\t\t\t</DataArray>\n");
			printf("\t\t</Polys>\n");
			printf("\t\t<CellData Scalar=\"geometry\">\n");
			printf("\t\t\t<DataArray type=\"Float32\" Name=\"slip amplitude\" NumberOfComponents=\"1\" format=\"ascii\">\n");
			printf("%f\n",$2);
			printf("\t\t\t</DataArray>\n");
			printf("\t\t\t<DataArray type=\"Float32\" Name=\"rake\" NumberOfComponents=\"1\" format=\"ascii\">\n");
			printf("%f\n",$6);
			printf("\t\t\t</DataArray>\n");
			printf("\t\t</CellData>\n");
			printf("\t</Piece>\n");
	} END {
		printf("</PolyData>\n");
		printf("</VTKFile>\n");
	}' - $TRIFILE >> $VTPFILE

}

usage(){
        echo "$self converts an .tri file to a vtp polygon"
	echo "file in xml format for visualization in Paraview"
	echo ""
	echo "usage: $self file.tri"
	echo "       $self file1.tri file2.tri file3.tri "
	echo ""

	exit
}

if [ $# -eq 0 ]; then
	usage
fi

while getopts "C:hs:x:y:z:" flag
do
	case "$flag" in
	C) Cset=1;CPT="$OPTARG";;
	h) hset=1;;
	s) sset=1;SCALE=$OPTARG;;
	x) xset=1;x=$OPTARG;;
	y) yset=1;y=$OPTARG;;
	z) zset=1;UTMZONE=$OPTARG;;
	esac
done
for item in $Cset $sset $xset $yset $zset; do
	shift;shift
done
for item in $hset; do
	shift;
done

if [ "1" != "$sset" ]; then
	SCALE=1
fi

if [ "1" != "$Cset" ]; then
	CPT=$selfdir/../share/jet.cpt
fi
echo $self: using color palette $(basename $CPT)

# loop over list of files to convert
while [ $# -ne 0 ]; do
	# define output file name (file.flt is converted to file.vtp, file.ext to file.ext.vtp)
	VTPFILE=$(dirname $1)/$(basename $1 .tri).vtp
	TRIFILE=$(dirname $1)/$(basename $1 .tri).tri
	NEDFILE=$(dirname $1)/$(basename $1 .tri).ned
	if [ ! -e $TRIFILE ]; then
		echo $self: could not find $TRIFILE. exiting.
		exit 2
	fi
	if [ ! -e $NEDFILE ]; then
		echo $self: could not find $NEDFILE. exiting.
		exit 2
	fi

	echo $self: converting $1 to $VTPFILE
	tri2vtp
	shift
done

