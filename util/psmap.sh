#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# script grdmap.sh
# display a map view of GMT .grd files with
# overlay from user-provided GMT scripts.
#
# Sylvain Barbot (01/01/07) - original form
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set -e
self=$(basename $0)
selfdir=$(dirname $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

usage(){
        echo "usage: $self [-b -5/5/-5/5] [-p -1.5/1.5/0.1] [...] file1.grd ... fileN.grd"
	echo "or"
        echo "       $self [-b -5/5/-5/5] [-p -1.5/1.5/0.1] [...] dir1/index1 ... dirN/indexN"
	echo ""
        echo "options:"
        echo "         -b bds sets the map bound to bds"
	echo "         -c palette_name [default jet]"
        echo "         -e file.sh runs file.sh to add content to map"
        echo "         -g switch to geographic projection (longitude/latitude)"
        echo "         -h display this error message and exit"
	echo "         -i file.grd illuminates the grd image with file.grd"
        echo "         -o output directory"
        echo "         -p palette sets the color scale bounds to palette"
        echo "         -r reverse the y-axis"
        echo "         -s step sets the distance between vectors"
        echo "         -u defines the color scale unit"
        echo "         -v vector sets the vector scale to vector"
	echo "         -t tick interval"
	echo "         -T title header"
	echo "         -x do not display map (only create .ps file)"
	echo "         -C interval plots contours every interval distance"
	echo "         -J overwrites the geographic projections (for -J o the -b option is relative)"
	echo "         -L draws a distance scale at coordinates lon/lat"
        echo "         -Y shift the plot vertically on the page"
	echo ""
	echo "Creates N maps based on grd files file1.grd ... fileN.grd, or based"
	echo "on files dir1/index1-{up,north.east}.grd ... dirN/indexN-{up,north,east}.grd"
	echo ""
        
        exit
}

my_gmt(){

	psbasemap -R$bds -J${PROJ} \
	    $AXIS \
	     -K -P -X1.2i -Y${YSHIFT}i \
	    > $PSFILE
	echo "psbasemap -R$bds -J$PROJ $AXIS -P -Y$YSHIFT"

	# running all required subprograms
	for subprog in $EXTRA; do
		OPTIONP=""
		if [ -e "$selfdir/$subprog" ]; then
			eval "$selfdir/$subprog $gset -b $bds -v $VECTOR $OPTIONP -H $HEIGHT $Jset $PSFILE"
		else
			if [ -e "$subprog" ]; then
				eval "$subprog $gset -b $bds -v $VECTOR $OPTIONP -H $HEIGHT $Jset $PROJ $PSFILE"
			fi
		fi
	done

	if [ "$VECTOR" != "-" ]; then
		UL=`echo $bds | awk -F "/" '{print $1,$4}' `
		pstext -O -K -J${PROJ} -N -R$bds \
			-G0/0/0 -Yr0.3i \
			<<EOF >> $PSFILE
$UL 14 0 4 LM $SIZE $unit
EOF
		psxy -O -K -J${PROJ} -R$bds -N \
			-W0.5p/0/0/0 -Xr0.9i \
			-Sv0.05/0.2/0.08n0.3c \
			<<EOF >> $PSFILE
$UL 0 1
EOF
		REVERT="-Xr-0.9i -Yr-0.3i"
	fi

	if [ "$Lset" != "" ]; then
		psbasemap -O -K -J$PROJ -R$bds $Lset >> $PSFILE
	fi

	if [ -e "$colors" ]; then
		#-Q0.20c/1.0c/0.4cn1.0c \
		psscale -O -K -B$PSSCALE/:$unit: -D3.5i/-0.8i/7.1i/0.2ih \
			-C$colors $REVERT $illuminationopt \
			>> $PSFILE 
		
		rm -f $colors
	fi

	echo 0 0 | psxy -O -J${PROJ} -R$bds -Sp0.001c >> $PSFILE

}

gmtset LABEL_FONT_SIZE 12p
gmtset HEADER_FONT_SIZE 12p
gmtset ANOT_FONT_SIZE 12p 
gmtset COLOR_BACKGROUND 255/255/255
gmtset COLOR_FOREGROUND 0/0/0
gmtset COLOR_NAN 255/255/255
gmtset PAPER_MEDIA archA

libdir=$(dirname $0)/../share
EXTRA=""

while getopts "b:c:e:ghi:o:p:v:s:t:T:u:xrC:J:L:Y:" flag
do
	case "$flag" in
	b) bset=1;bds=$OPTARG;;
	c) cset="-c";carg=$OPTARG;;
	e) eset=1;EXTRA="$EXTRA $OPTARG";;
	g) gset="-g";;
	h) hset=1;;
	i) iset=1;illumination="-I$OPTARG";illuminationopt="-I";;
	o) oset=1;ODIR=$OPTARG;;
	p) pset=1;P=$OPTARG;PSSCALE=`echo $OPTARG | awk -F"/" 'function abs(x){return x<0?-x:x}{print abs($2-$1)/6}'`;;
	r) rset=1;;
	s) sset=1;ADX=$OPTARG;;
	t) tset=1;tick=$OPTARG;;
	T) Tset=1;title=$OPTARG;;
	u) uset=1;unit=$OPTARG;;
	v) vset=1;SIZE=$OPTARG;VECTOR=$OPTARG"c";;
	x) xset=1;;
	C) Cset="-C";contour=$OPTARG;;
	J) Jset="-J";PROJ=$OPTARG;;
	L) Lset="-L$OPTARG";;
	Y) Yset=1;Yshift=$OPTARG;;
	esac
done
for item in $bset $cset $iset $oset $pset $vset $sset $tset $Tset $uset $Yset $Cset $Oset $Jset $Lset $EXTRA;do
	shift;shift
done
for item in $gset $hset $xset $rset;do
	shift;
done

#-------------------------------------------------------- 
# DEFAULTS
#-------------------------------------------------------- 

if [ "$vset" != "1" ]; then
	# unused value but preserve the number of elements in call to subroutine
	VECTOR="-"
fi

# vertical shift of plot
if [ "$Yset" != "1" ]; then
	YSHIFT=2.0
else
	YSHIFT=`echo $Yshift | awk '{print $1+2}'`
fi

# color scale
if [ "$cset" != "-c" ]; then
	cptfile=$libdir/jet
else
	if [ -e "$libdir/$carg" ]; then 
		cptfile=$libdir/$carg
	else 
		cptfile=$carg;
	fi
fi

# usage
if [ $# -lt "1" -a "$Oset" != "1" ] || [ "$hset" == "1" ] ; then
	usage
else
	echo $self: using colorfile $cptfile
fi

if [ "$1" != "" ]; then
	PSFILE=$(dirname $1)/$(basename $1 .ps).ps;
	WDIR=$(dirname $PSFILE)
	ODIR=$WDIR
	TITLE=$(basename $PSFILE .ps)
	#PSFILE=$ODIR/plot.ps
	PDFFILE=$ODIR/$(basename $PSFILE .ps).pdf
	
	if [ "$pset" == "1" ]; then
		colors=$WDIR/palette.cpt
	fi
	if [ "$bset" != "1" ]; then
		echo "$self: -b option should be set with the -E option. exiting."
		echo ""
		usage
	else

		# tick marks
		if [ "$tset" != "1" ]; then
			tick=`echo $bds | awk -F "/" '{s=1;print ($2-$1)/s/7}'`
		fi
		tickf=`echo $tick | awk '{print $1/2}'`

		# Cartesian vs geographic coordinates
		if [ "$gset" != "-g" ]; then
			if [ "$Jset" == "" ]; then
				HEIGHT=`echo $bds | awk -F "/" '{printf("%fi\n",($4-$3)/($2-$1)*7)}'`
				if [ "$rset" != "1" ]; then
					PROJ="X7i/"$HEIGHT
				else
					PROJ="X7i/"-$HEIGHT
				fi
			else
				HEIGHT="-"
			fi
		        AXIS=-Ba${tick}:"":/a${tick}:""::."$TITLE":WSne
		else
			HEIGHT=7i
			PROJ="M0/0/$HEIGHT"
		        AXIS=-B${tick}:"":/${tick}:""::."$TITLE":WSne
		fi
	fi
fi

# color bounds
if [ "$pset" == "1" ]; then
	makecpt -C$cptfile -T$P -D > $colors;
	m=`echo $P | awk -F"/" 'function max(x,y){return (x>y)?x:y};function abs(x){return (0<x)?x:-x}{print max(abs($1),$2)}'`
fi

my_gmt

# add trailer information
echo %% Postscript created with >> $PSFILE
echo %% $(basename $0) $cmdline >> $PSFILE

# open file for display
echo $self": created map "$PSFILE

if [ "$xset" != "1" ]; then
	#display -trim $PSFILE &
	#gv -geometry +0+0 -spartan -scale=0.5 $PSFILE &
	ps2pdf -sPAPERSIZE="archA" -dPDFSETTINGS=/prepress $PSFILE $PDFFILE
	echo $self": converted to pdf file "$PDFFILE
	xpdf -geometry +0+0 -paper "archA" $PDFFILE -z 100 -g 565x655 >& /dev/null &
fi

if [ "$#" != "0" ];then
	shift
	if [ "$#" != "0" ];then
		echo ""
	fi
fi

