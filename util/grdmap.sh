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
	echo "         -c palette_name [default my_jet]"
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
	echo "         -E file.ps only plot base map with extra scripts"
	echo "         -J overwrites the geographic projections"
        echo "         -Y shift the plot vertically on the page"
	echo ""
	echo "Creates N maps based on grd files file1.grd ... fileN.grd, or based"
	echo "on files dir1/index1-{up,north.east}.grd ... dirN/indexN-{up,north,east}.grd"
	echo ""
        
        exit
}

my_gmt(){

	if [ -e "$U3" ]; then
		grdimage $U3 -R$bds -J${PROJ} -Q \
		    $AXIS \
		    -K -C$colors -P -X1.2i -Y${YSHIFT}i $illumination \
		    > $PSFILE
		if [ "$Cset" == "-C" ]; then
			grdcontour -K -O $U3 -W0.3p/100/100/100 -R$bds -J${PROJ} -S -C$contour -P >> $PSFILE
		fi
	else
		psbasemap -R$bds -J${PROJ} \
		    $AXIS \
		     -K -P -X1.2i -Y${YSHIFT}i \
		    > $PSFILE
	fi

	# running all required subprograms
	for subprog in $EXTRA; do
		if [ "" != "$U3" ]; then
			OPTIONP="-p $U3"
		else
			OPTIONP=""
		fi
		if [ -e "$selfdir/$subprog" ]; then
			#echo $self: running $subprog $PSFILE $bds $VECTOR $U3 $HEIGHT
			eval "$subprog $gset -b $bds -v $VECTOR $OPTIONP -H $HEIGHT $PSFILE"
		else
			if [ -e "$subprog" ]; then
				#echo $self: running $subprog $PSFILE $bds $VECTOR $U3 $HEIGHT
				eval "$subprog $gset -b $bds -v $VECTOR $OPTIONP -H $HEIGHT $PSFILE"
			fi
		fi
	done

	# plotting vectors
	if [ -e "$U1" ]; then

		#echo $self": found "$U1": plotting vectors"
		echo $self": Using VECTOR="$VECTOR", STEP="$ADX

		# arrowwidth/headlength/headwidth
		# grdvector crashes with wrong sampling
		#	-Q0.51c/0.8c/0.7cn1.0c \
		grdvector $U1 $U2 -K -J${PROJ} -O -R$bds -I$ADX/$ADX \
			-Q0.3c/0.5c/0.4cn1.0c \
			-S$VECTOR -W0.1p/0/0/0 \
			-G255/255/255 \
			 >> $PSFILE
	fi

	# plot the vector legend if $VECTOR is set
	if [ "$VECTOR" != "-" ]; then
		UL=`echo $bds | awk -F "/" '{print $1,$4}' `
		pstext -O -K -J${PROJ} -N -R$bds \
			-G0/0/0 -Yr0.3i \
			<<EOF >> $PSFILE
$UL 14 0 4 LM $SIZE $unit
EOF
		psxy -O -K -J${PROJ} -R$bds -N \
			-W0.5p/0/0/0 -Xr0.9i \
			-Sv0.2c/1.0c/0.4cn1.0c \
			<<EOF >> $PSFILE
$UL 0 1
EOF
		REVERT="-Xr-0.9i -Yr-0.3i"
	fi

	if [ -e "$colors" ]; then
		#-Q0.20c/1.0c/0.4cn1.0c \
		psscale -O -K -B$PSSCALE/:$unit: -D2.0i/-0.8i/3.1i/0.2ih \
			-C$colors $REVERT \
			>> $PSFILE 
		
		rm -f $colors
	fi

	echo 0 0 | psxy -O -J${PROJ} -R$bds -Sp0.001c >> $PSFILE

}

gmtset LABEL_FONT_SIZE 12p
gmtset HEADER_FONT_SIZE 12p
gmtset ANOT_FONT_SIZE 12p 
gmtset COLOR_BACKGROUND 0/0/255
gmtset COLOR_FOREGROUND 255/0/0
gmtset COLOR_NAN 255/255/255
gmtset PAPER_MEDIA a5

libdir=$(dirname $0)
EXTRA=""

while getopts "b:c:e:ghi:o:p:v:s:t:T:u:xrC:E:J:Y:" flag
do
	case "$flag" in
	b) bset=1;bds=$OPTARG;;
	c) cset="-c";carg=$OPTARG;;
	e) eset=1;EXTRA="$EXTRA $OPTARG";;
	g) gset="-g";;
	h) hset=1;;
	i) iset=1;illumination="-I$OPTARG";;
	o) oset=1;ODIR=$OPTARG;;
	p) pset=1;P=$OPTARG;PSSCALE=`echo $OPTARG | awk -F"/" 'function abs(x){return x<0?-x:x}{print abs($2-$1)/4}'`;;
	r) rset=1;;
	s) sset=1;ADX=$OPTARG;;
	t) tset=1;tick=$OPTARG;;
	T) Tset=1;title=$OPTARG;;
	u) uset=1;unit=$OPTARG;;
	v) vset=1;SIZE=$OPTARG;VECTOR=$OPTARG"c";;
	x) xset=1;;
	C) Cset="-C";contour=$OPTARG;;
	E) Eset=1;PSFILE=$(dirname $OPTARG)/$(basename $OPTARG .ps).ps;;
	J) Jset=1;PROJ=$OPTARG;;
	Y) Yset=1;Yshift=$OPTARG;;
	esac
done
for item in $bset $cset $iset $oset $pset $vset $sset $tset $Tset $uset $Yset $Cset $Eset $Jset $EXTRA;do
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
	#cptfile=hot
	cptfile=$libdir/my_jet
	#cptfile=$libdir/my_hot_inv
else
	if [ -e "$libdir/$carg" ]; then 
		cptfile=$libdir/$carg
	else 
		cptfile=$carg;
	fi
fi

# usage
if [ $# -lt "1" -a "$Eset" != "1" ] || [ "$hset" == "1" ] ; then
	usage
else
	echo $self: using colorfile $cptfile
fi


# loop over grd files
while [ "$#" != "0" -o "$Eset" == "1" ];do

	if [ "$1" != "" ]; then
		WDIR=`dirname $1`
		INDEX=`basename $1`

		# default names
		U1=$WDIR/$INDEX-east.grd
		U2=$WDIR/$INDEX-north.grd
		U3=$WDIR/$INDEX-up.grd

		if [ -e "$U3" ]; then
			U3=$U3;
		else
			U3=$WDIR/$INDEX
		fi
		if [ "$bset" != "1" ]; then
			bds=`grdinfo -C $U3 | awk '{s=1;print $2/s"/"$3/s"/"$4/s"/"$5/s}'`
		fi

		# tick marks
		if [ "$tset" != "1" ]; then
			tick=`echo $bds | awk -F "/" '{s=1;print ($2-$1)/s/4}'`
		fi
		tickf=`echo $tick | awk '{print $1/2}'`

		echo $self": Using directory "$WDIR", plotting index "$INDEX

		# defaults

		# output directory
		if [ "$Tset" != "1" ]; then
			title=$U3
		fi

		# output directory
		if [ "$oset" != "1" ]; then
			ODIR=$WDIR
		fi

		# units
		if [ "$uset" != "1" ]; then
			# retrieve the "Remarks"
			# use grdedit -D:=:=:=:=:=:=:cm: file.grd to update the value
			unit=`grdinfo $U3 | grep "Remark" | awk -F ": " '{print $3}'`
		fi

		colors=$ODIR/palette.cpt
		if [ "$pset" != "1" ]; then
			# cool, copper, gebco, globe, gray, haxby, hot, jet, no_green, ocean
			# polar, rainbow, red2green, relief, topo, sealand, seis, split, wysiwyg  
			PSSCALE=`grdinfo $U3 -C | awk 'function abs(x){return x<0?-x:x}{if (abs($6) >= abs($7)) print abs($6)/2; else print abs($7)/2}'`
			if [ "0" == $PSSCALE ]; then
				grd2cpt $U3 -C$cptfile -Z -T= -L-1/1 > $colors
				PSSCALE=0.5
			else
				grd2cpt $U3 -C$cptfile -Z -T= > $colors
			fi
		fi

		if [ -e "$U1" ]; then
			if [ "$vset" != "1" ]; then
		                MAX1=`grdinfo $U1 -C | awk '{if (sqrt($7^2)>sqrt($6^2)){print sqrt($7^2)}else{print sqrt($6^2)}}'`
		                MAX2=`grdinfo $U2 -C | awk '{if (sqrt($7^2)>sqrt($6^2)){print sqrt($7^2)}else{print sqrt($6^2)}}'`
				SIZE=`echo "$MAX1 $MAX2"| awk '{if (0==$1 && 0==$2){print 1}else{print ($1+$2)*0.95}}'`
				VECTOR=$SIZE"c"
			fi
			if [ "$sset" != "1" ]; then
				ADX=`grdinfo $U2 -C | awk '{printf "5*%11.9f\n", $8}' | bc`
			fi
		fi
	
		if [ "$gset" != "-g" ]; then
			if [ "$Jset" == "" ]; then
				# Cartesian coordinates
				HEIGHT=`echo $bds | awk -F "/" '{printf("%fi\n",($4-$3)/($2-$1)*4)}'`
				if [ "$rset" != "1" ]; then
					PROJ="X4i/"$HEIGHT
				else
					PROJ="X4i/"-$HEIGHT
				fi
			else
				HEIGHT="-"
			fi
			AXIS=-Bf${tickf}a${tick}:"":/f${tickf}a${tick}:""::."$title":WSne
		else
			# geographic coordinates
			HEIGHT=4i
			PROJ="M$HEIGHT"
		        AXIS=-B${tick}:"":/${tick}:""::."$title":WSne
		fi

		echo $self": z-min/z-max for "$U3": "`grdinfo -C $U3 | awk '{print $6,$7}'`
	
		PSFILE=$ODIR/$INDEX-plot.ps
		PDFFILE=$ODIR/$INDEX-plot.pdf
	
	else
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
				tick=`echo $bds | awk -F "/" '{s=1;print ($2-$1)/s/4}'`
			fi
			tickf=`echo $tick | awk '{print $1/2}'`

			# Cartesian vs geographic coordinates
			if [ "$gset" != "-g" ]; then
				HEIGHT=`echo $bds | awk -F "/" '{printf("%fi\n",($4-$3)/($2-$1)*4)}'`
				if [ "$rset" != "1" ]; then
					PROJ="X4i/"$HEIGHT
				else
					PROJ="X4i/"-$HEIGHT
				fi
			        AXIS=-Ba${tick}:"":/a${tick}:""::."$TITLE":WSne
			else
				HEIGHT=4i
				PROJ="M0/0/$HEIGHT"
			        AXIS=-B${tick}:"":/${tick}:""::."$TITLE":WSne
			fi
		fi
	fi

	# color bounds
	if [ "$pset" == "1" ]; then
		makecpt -C$cptfile -T$P -D > $colors;
	fi

	my_gmt $INDEX
	
	# add trailer information
	echo %% Postscript created with >> $PSFILE
	echo %% $(basename $0) $cmdline >> $PSFILE
	
	# open file for display
	echo $self": Created map "$PSFILE
	
	if [ "$xset" != "1" ]; then
		#display -trim $PSFILE &
		#gv -geometry +0+0 -spartan $PSFILE &
		ps2pdf -sPAPERSIZE="a4" $PSFILE $PDFFILE
		echo $self": Converted to pdf file "$PDFFILE
		xpdf -geometry +0+0 -paper "A5" $PDFFILE -z 100 -g 565x655 >& /dev/null &
	fi
	
	if [ "$#" != "0" ];then
		shift
		if [ "$#" != "0" ];then
			echo ""
		fi
	fi

	# prevent more empty plots (cancel -E option)
	Eset=""
done
