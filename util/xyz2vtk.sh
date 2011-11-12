#!/bin/sh

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

if [ $# -eq 0 ]; then
        echo "$self converts an xyz file to a vtp line polydata"
	echo "file in xml format for visualization in Paraview"
	echo ""
	echo "usage: $self [-c com] file.xyz "
	echo "       $self file1.xyz file2.xyz file3.xyz"
	echo ""
	echo "options:"
	echo "         -c com sets the segment separator [default >]"
	echo ""
	echo "file foo.xyz is converted to foo.vtp"
	echo "files foo.ext are converted to foo.ext.vtp"
	echo "note that the 3rd column of the .xyz file is"
	echo "interpreted as height in Paraview."
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

XYZFILE=$1

if [ ! -e $XYZFILE ]; then
	echo $self: could not find $XYZFILE. exiting.
	exit 2
fi

EXT=`echo $1 | awk '{print substr($0,length($0)-2,3)}'`

# define output file name (file.xyz is converted to file.vtp, file.ext to file.ext.vtp)
VTKFILE=$(dirname $1)/$(basename $1 .xyz).vtp

echo $self: converting $1 to $VTKFILE

cat $XYZFILE | awk -v c=$com '
	BEGIN{
	n=0;line="";m=0;
	printf("<?xml version=\"1.0\"?>\n");
	printf("<VTKFile type=\"PolyData\" version=\"0.1\">\n");
	printf("  <PolyData>\n");
	}
	{
	if (substr($1,0,1)==c){
		if (0!=n){
			printf("    <Piece NumberOfPoints=\"%i\" NumberOfLines=\"1\">\n",n);
			printf("      <Points>\n");
			printf("        <DataArray type=\"Float32\" Name=\"Fault Segment %i\" NumberOfComponents=\"3\" format=\"ascii\">\n",m);
			print line; 
			printf("        </DataArray>\n");
			printf("      </Points>\n");
			printf("      <Lines>\n");
			printf("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%i\">\n",n-1);
			for (i=0;i<n;i++){printf("%i ",i)};
			printf("\n");
			printf("        </DataArray>\n");
			printf("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"%i\" RangeMax=\"%i\">\n",n,n);
			printf("          %i\n",n);
			printf("        </DataArray>\n");
			printf("      </Lines>\n");
			printf("    </Piece>\n");
			n=0;
			line=""
			m=m+1;
		}
	}
	else{
		line=line " " $2" "$1" "$3;
		n=n+1;
	}
	}
	END{
	if (n!=0){
		printf("    <Piece NumberOfPoints=\"%i\" NumberOfLines=\"1\">\n",n);
		printf("      <Points>\n");
		printf("        <DataArray type=\"Float32\" Name=\"Fault Segment %i\" NumberOfComponents=\"3\" format=\"ascii\">\n",m);
		print line; 
		printf("        </DataArray>\n");
		printf("      </Points>\n");
		printf("      <Lines>\n");
		printf("        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\" RangeMin=\"0\" RangeMax=\"%i\">\n",n-1);
		for (i=0;i<n;i++){printf("%i ",i)};
		printf("\n");
		printf("        </DataArray>\n");
		printf("        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" RangeMin=\"%i\" RangeMax=\"%i\">\n",n,n);
		printf("          %i\n",n);
		printf("        </DataArray>\n");
		printf("      </Lines>\n");
		printf("    </Piece>\n");
	}
	printf("  </PolyData>\n");
	printf("</VTKFile>\n");
	}' > $VTKFILE

	shift
done
