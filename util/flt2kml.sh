#!/bin/bash

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

flt2kml(){
	# save input to variable for multiple access
	IN=`cat $FLTFILE`
	
	# first convert slip to color using cpt2rgb.sh
	echo "$IN" | grep -v "#" | awk '{if (10==NF){print $2}else{print 1}}' | cpt2rgb.sh $CPT | awk -v name=$KMLFILE '
		BEGIN{
			print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
			print "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">";
			print "<Document>";
			printf "\t<name>%s</name>\n",name;
		}
		{
		printf "\t<Style id=\"flt2kml_%05.5d\">\n",NR;
		print "\t\t<ListStyle>";
		print "\t\t\t<listItemType>checkHideChildren</listItemType>";
	        print "\t\t</ListStyle>";
		print "\t\t<LineStyle>";
		printf "\t\t\t<color>bb%02x%02x%02x</color>\n",$3,$2,$1;
		print "\t\t</LineStyle>";
		print "\t\t<PolyStyle>";
		# aabbggrr
		printf "\t\t\t<color>cc%02x%02x%02x</color>\n",$3,$2,$1;
		print "\t\t</PolyStyle>";
		print "\t</Style>";
		}' > $KMLFILE

	echo "$IN" | grep -v "#" | awk '
		BEGIN{
			pi=atan2(1,0)*2;
		} {
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
			} else {
			x1=$3;
			x2=$4;
			x3=$5;
			L=$6;
			W=$7;
			strike=$8*pi/180;
			dip=$9*pi/180;
			rake=$10*pi/180;
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
        	        printf("%f %f\n",x2              ,x1              );
        	        printf("%f %f\n",x2       +s[1]*L,x1       +s[0]*L);
			printf("%f %f\n",x2-d[1]*W+s[1]*L,x1-d[0]*W+s[0]*L);
        	        printf("%f %f\n",x2-d[1]*W       ,x1-d[0]*W       );
        	        printf("%f %f\n",x2              ,x1              );
		}' | \
		awk -v x=$x -v y=$y -v scale=$SCALE '{$1=$1*scale+x;$2=$2*scale+y;printf "%f %f\n", $1,$2}' | \
		invproj +proj=utm +zone=$UTMZONE -f "%f" | awk '{
			printf "%f,%f,0 ",$1,$2
			if (4==((NR-1)%5)){
				print "";
			}
		}' | awk ' 
		{
			print "\t<Placemark>";
			print "\t\t<name>TEST1</name>"
			printf "\t\t<styleUrl>#flt2kml_%05.5d</styleUrl>\n",NR;
			print "\t\t<Polygon>"
			print "\t\t\t<tessellate>1</tessellate>"
			print "\t\t\t<outerBoundaryIs>"
			print "\t\t\t\t<LinearRing>"
			print "\t\t\t\t\t<coordinates>"
			print "\t\t\t\t\t\t",$0
			print "\t\t\t\t\t</coordinates>"
			print "\t\t\t\t</LinearRing>"
			print "\t\t\t</outerBoundaryIs>"
			print "\t\t</Polygon>"
			print "\t</Placemark>";
		} END {
			print "</Document>";
			print "</kml>";
		}' >> $KMLFILE

}

usage(){
        echo "$self converts an .flt file to a kml polygon"
	echo "file in xml format for visualization in Google Earth"
	echo ""
	echo "usage: $self file.flt"
	echo "       $self file1.flt file2.flt file3.flt "
	echo ""
	echo "or from the standard input"
	echo ""
	echo "       cat file.flt | $self > file.kml"
	echo ""
	echo "file foo.flt is converted to foo.kml"
	echo "files foo.ext are converted to foo.ext.kml"
	echo ""

	exit
}

if [ -t 0 ] && [ $# -eq 0 ]; then
	usage
fi

while getopts "C:hs:x:y:z:" flag
do
	case "$flag" in
	C) Cset=1;CPT="-C $OPTARG";;
	h) hset=1;;
	s) sset=1;scale=$OPTARG;;
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

if [ ! -t 0 ]; then
	FLTFILE=-
	KMLFILE=/dev/stdout
	flt2kml
else
	# loop over list of files to convert
	while [ $# -ne 0 ]; do
		# define output file name (file.flt is converted to file.kml, file.ext to file.ext.kml)
		KMLFILE=$(dirname $1)/$(basename $1 .flt).kml
		FLTFILE=$1
		if [ ! -e $FLTFILE ]; then
			echo $self: could not find $FLTFILE. exiting.
			exit 2
		fi
		echo $self: converting $1 to $KMLFILE
		flt2kml
		shift
	done
fi

