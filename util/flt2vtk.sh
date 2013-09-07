#!/bin/bash

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

flt2vtk(){
	grep -v "#" $FLTFILE | awk '
	BEGIN{
		pi=atan2(1,0)*2;
		printf("<?xml version=\"1.0\"?>\n");
		printf("<VTKFile type=\"PolyData\" version=\"0.1\">\n");
		printf("  <PolyData>\n");
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
		b[0]=+cos(strike)*cos(rake)+sin(strike)*cos(dip)*sin(rake);
		b[1]=+sin(strike)*cos(rake)-cos(strike)*cos(dip)*sin(rake);
		b[2]=-sin(dip)*sin(rake);
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
		printf("      <CellData Scalar=\"geometry\">\n");
		printf("        <DataArray type=\"Float32\" Name=\"strike\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		printf("%f\n",((($nstrike+180)%360-360)%360+180));
		printf("        </DataArray>\n");
		printf("        <DataArray type=\"Float32\" Name=\"dip\" NumberOfComponents=\"1\" format=\"ascii\">\n");
		printf("%f\n",$ndip);
		printf("        </DataArray>\n");
		if (10==NF){
		printf("        <DataArray type=\"Float32\" Name=\"slip\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		printf("%f %f %f\n",slip*b[0],slip*b[1],slip*b[2]);
		printf("        </DataArray>\n");
		} else {
		printf("        <DataArray type=\"Float32\" Name=\"unit slip\" NumberOfComponents=\"3\" format=\"ascii\">\n");
		printf("%f %f %f\n",b[0],b[1],b[2]);
		printf("        </DataArray>\n");
		}
		printf("      </CellData>\n");
		printf("    </Piece>\n");
	}
	END{
		printf("  </PolyData>\n");
		printf("</VTKFile>\n");
	}' > $VTKFILE
}

usage(){
        echo "$self converts an .flt file to a vtp polygon"
	echo "file in xml format for visualization in Paraview"
	echo ""
	echo "usage: $self file.flt"
	echo "       $self file1.flt file2.flt file3.flt "
	echo ""
	echo "or from the standard input"
	echo ""
	echo "       cat file.flt | $self > file.vtp"
	echo ""
	echo "file foo.flt is converted to foo.vtp"
	echo "files foo.ext are converted to foo.ext.vtp"
	echo ""
}

if [ -t 0 ] && [ $# -eq 0 ]; then
	usage
	exit
fi

if [ ! -t 0 ]; then
	FLTFILE=-
	VTKFILE=/dev/stdout
	flt2vtk
else
	# loop over list of files to convert
	while [ $# -ne 0 ];do
		# define output file name (file.flt is converted to file.vtp, file.ext to file.ext.vtp)
		VTKFILE=$(dirname $1)/$(basename $1 .flt).vtp
		FLTFILE=$1
		if [ ! -e $FLTFILE ]; then
			echo $self: could not find $FLTFILE. exiting.
			exit 2
		fi
		echo $self: converting $1 to $VTKFILE
		flt2vtk
		shift
	done
fi

