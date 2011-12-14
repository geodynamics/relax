#!/bin/sh

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if [ $# -eq 0 ]; then
        echo "$self converts an .flt file to a vtp polygon"
	echo "file in xml format for visualization in Paraview"
	echo ""
	echo "usage: $self [-c com] file.flt"
	echo "       $self file1.flt file2.flt file3.flt "
	echo ""
	echo "file foo.flt is converted to foo.vtp"
	echo "files foo.ext are converted to foo.ext.vtp"
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

# define output file name (file.flt is converted to file.vtp, file.ext to file.ext.vtp)
VTKFILE=$(dirname $1)/$(basename $1 .flt).vtp

echo $self: converting $1 to $VTKFILE

grep -v "#" $FLTFILE | awk '
	BEGIN{
		pi=atan2(1,0)*2;
		printf("<?xml version=\"1.0\"?>\n");
		printf("<VTKFile type=\"PolyData\" version=\"0.1\">\n");
		printf("  <PolyData>\n");
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
		printf("    <Piece NumberOfPoints=\"4\" NumberOfPolys=\"1\">\n");
		printf("      <Points>\n");
		printf("        <DataArray type=\"Float32\" Name=\"Fault Patch\" NumberOfComponents=\"3\" format=\"ascii\">\n");
                printf("%f %f %f\n",x1              ,x2              ,x3              );
                printf("%f %f %f\n",x1       +s[0]*L,x2       +s[1]*L,x3       +s[2]*L);
		printf("%f %f %f\n",x1-d[0]*W+s[0]*L,x2-d[1]*W+s[1]*L,x3-d[2]*W+s[2]*L);
                printf("%f %f %f\n",x1-d[0]*W       ,x2-d[1]*W       ,x3-d[2]*W       );
		printf("        </DataArray>\n");
		printf("      </Points>\n");
		printf("      <Polys>\n");
		printf("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%i\">\n",n-1);
		printf("0 1 2 3\n");
		printf("        </DataArray>\n");
		printf("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"4\" RangeMax=\"4\">\n");
		printf("          4\n");
		printf("        </DataArray>\n");
		printf("      </Polys>\n");
		printf("    </Piece>\n");
	}
	END{
		printf("  </PolyData>\n");
		printf("</VTKFile>\n");
	}' > $VTKFILE

	shift
done
