#!/bin/bash

set -e
self=$(basename $0)
selfdir=$(dirname $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

tri2kml(){
	# first convert slip to color using cpt2rgb.sh
	grep -v "#" $TRIFILE | awk '{print $2}' | cpt2rgb.sh -C $CPTFILE | awk -v name=$KMLFILE '
		BEGIN{
			print "<?xml version=\"1.0\" encoding=\"UTF-8\"?>";
			print "<kml xmlns=\"http://www.opengis.net/kml/2.2\" xmlns:gx=\"http://www.google.com/kml/ext/2.2\" xmlns:kml=\"http://www.opengis.net/kml/2.2\" xmlns:atom=\"http://www.w3.org/2005/Atom\">";
			print "<Document>";
			printf "\t<name>%s</name>\n",name;
			print "\t<Style>";
			print "\t\t<ListStyle>";
			print "\t\t\t<listItemType>checkHideChildren</listItemType>";
	        	print "\t\t</ListStyle>";
			print "\t</Style>";
		}
		{
		printf "\t<Style id=\"tri2kml_%05.5d\">\n",NR;
		print "\t\t<LineStyle>";
		printf "\t\t\t<color>bb%02x%02x%02x</color>\n",$3,$2,$1;
		print "\t\t</LineStyle>";
		print "\t\t<PolyStyle>";
		# aabbggrr
		printf "\t\t\t<color>cc%02x%02x%02x</color>\n",$3,$2,$1;
		print "\t\t</PolyStyle>";
		print "\t</Style>";
		}' > $KMLFILE


		grep -v "#" $NEDFILE | \
			awk -v x=$x -v y=$y -v scale=$SCALE '{$3=$3*scale+x;$2=$2*scale+y;printf "%f %f %d\n", $3,$2,$1}' | \
			invproj +proj=utm +zone=$UTMZONE -f "%f" | awk '{printf "%f %f %d\n",$1,$2,$3}' | \
			awk 'NR==FNR {for (j=1;j<=2;j++){a[j+0,0+$3]=$j}; next} { 
			if ("#" != substr($0,0,1)) {
				print "\t<Placemark>";
				print "\t\t<name>TEST1</name>"
				printf "\t\t<styleUrl>#tri2kml_%05.5d</styleUrl>\n",$1;
				print "\t\t<Polygon>"
				print "\t\t\t<tessellate>1</tessellate>"
				print "\t\t\t<outerBoundaryIs>"
				print "\t\t\t\t<LinearRing>"
				print "\t\t\t\t\t<coordinates>"
				print "\t\t\t\t\t\t",a[1,0+$3]","a[2,0+$3]",0 ",a[1,0+$4]","a[2,0+$4]",0 ",a[1,0+$5]","a[2,0+$5]",0 ",a[1,0+$3]","a[2,0+$3]",0"
				print "\t\t\t\t\t</coordinates>"
				print "\t\t\t\t</LinearRing>"
				print "\t\t\t</outerBoundaryIs>"
				print "\t\t</Polygon>"
				print "\t</Placemark>";
			}
		} END {
			print "</Document>";
			print "</kml>";
		}' - $TRIFILE >> $KMLFILE

}

usage(){
        echo "$self converts an .tri file to a kml polygon"
	echo "file in xml format for visualization in Google Earth"
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
	# define output file name (file.flt is converted to file.kml, file.ext to file.ext.kml)
	KMLFILE=$(dirname $1)/$(basename $1 .tri).kml
	TRIFILE=$(dirname $1)/$(basename $1 .tri).tri
	NEDFILE=$(dirname $1)/$(basename $1 .tri).ned
	CPTFILE=$(dirname $1)/palette.cpt
	if [ ! -e $TRIFILE ]; then
		echo $self: could not find $TRIFILE. exiting.
		exit 2
	fi
	if [ ! -e $NEDFILE ]; then
		echo $self: could not find $NEDFILE. exiting.
		exit 2
	fi

	# create color palette $CPTFILE
	T=`minmax -C $TRIFILE | awk '{print $3"/"$4"/"(($4-$3)/30)}'`
	makecpt -C$CPT -T$T -Z -N > $CPTFILE

	echo $self: converting $1 to $KMLFILE
	tri2kml
	rm -f $CPTFILE
	shift
done

