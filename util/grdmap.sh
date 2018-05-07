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
	echo "         -l clip.xy clip the background color and vector plot"
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
	echo "         -H fixed height"
	echo "         -L draws a distance scale at coordinates lon/lat"
	echo "         -O file.ps only plot base map with extra scripts"
	echo "         -J overwrites the geographic projections (for -J o the -b option is relative)"
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

	if [ "" != "$lset" ]; then
		psclip -O -K -J${PROJ} -R$bds -m $clip >> $PSFILE
	fi

	if [ -e "$U3" ]; then
		grdimage $U3 -R$bds -J${PROJ} -Q \
		    -O -K -C$colors -P $illumination \
		    >> $PSFILE
		if [ "$Cset" == "-C" ]; then
			grdcontour -K -O $U3 -W0.3p/100/100/100 -R$bds -J${PROJ} -S4 -C$contour -P >> $PSFILE
		fi
	fi

	if [ "" != "$lset" ]; then
		psclip -O -K -J${PROJ} -R$bds -C >> $PSFILE
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
			eval "$selfdir/$subprog $gset -b $bds -v $VECTOR $OPTIONP -H $HEIGHT $Jset $PROJ $PSFILE"
		else
			if [ -e "$subprog" ]; then
				#echo $self: running $subprog $gset -b $bds -v $VECTOR $OPTIONP -H $HEIGHT $Jset $PROJ $PSFILE
				eval "$subprog $gset -b $bds -v $VECTOR $OPTIONP -H $HEIGHT $Jset $PROJ $PSFILE"
			fi
		fi
	done

	# plotting vectors
	if [ -e "$U1" ]; then

		#echo $self": found "$U1": plotting vectors"
		echo $self": Using VECTOR="$VECTOR", STEP="$ADX

		if [ "" != "$lset" ]; then
			psclip -O -K -J${PROJ} -R$bds -m $clip >> $PSFILE
		fi

		# arrowwidth/headlength/headwidth
		# grdvector crashes with wrong sampling
		#	-Q0.51c/0.8c/0.7cn1.0c \
	        #	-Q0.3c/0.5c/0.4cn1.0c \
		grdvector $U1 $U2 -K -J${PROJ} -O -R$bds -I$ADX/$ADX \
			-Q0.05/0.2/0.08n0.3c \
			-S$VECTOR -W0.5p/0/0/0 \
			-G255/255/255 \
			 >> $PSFILE

		if [ "" != "$lset" ]; then
			psclip -O -K -J${PROJ} -R$bds -C >> $PSFILE
		fi

	fi

	# plot the vector legend if $VECTOR is set
	if [ "$VECTOR" != "-" ]; then
		UL=`echo $bds | awk -F "/" '{print $1,$4}' `
		pstext -O -K -J${PROJ} -N -R$bds \
			-G0/0/0 -Yr+0.3i \
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
			$TRIANGLES -C$colors $REVERT $illuminationopt \
			>> $PSFILE 
		
		rm -f $colors
	fi

	echo 0 0 | psxy -O -J${PROJ} -R$bds -Sp0.001c >> $PSFILE

}

gmtset LABEL_FONT_SIZE 12p
gmtset HEADER_FONT_SIZE 12p
gmtset ANOT_FONT_SIZE 12p 
gmtset COLOR_BACKGROUND 0/0/0
gmtset COLOR_FOREGROUND 255/255/255
gmtset COLOR_NAN 255/255/255
gmtset PAPER_MEDIA archA

libdir=$(dirname $0)/../share
EXTRA=""

while getopts "b:c:e:ghi:l:o:p:v:s:t:u:xrC:HJ:L:O:T:Y:" flag
do
	case "$flag" in
	b) bset=1;bds=$OPTARG;;
	c) cset="-c";carg=$OPTARG;;
	e) eset=1;EXTRA="$EXTRA $OPTARG";;
	g) gset="-g";;
	h) hset=1;;
	i) iset=1;illumination="-I$OPTARG";illuminationopt="-I";;
	l) lset=1;clip=$OPTARG;;
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
	H) Hset="-H";;
	L) Lset="-L$OPTARG";;
	O) Oset=1;PSFILE=$(dirname $OPTARG)/$(basename $OPTARG .ps).ps;;
	J) Jset="-J";PROJ="$OPTARG";;
	Y) Yset=1;Yshift=$OPTARG;;
	esac
done
for item in $bset $cset $iset $lset $oset $pset $vset $sset $tset $Tset $uset $Yset $Cset $Lset $Oset $Jset $EXTRA;do
	shift;shift
done
for item in $gset $hset $xset $rset $Hset;do
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


# loop over grd files
while [ "$#" != "0" -o "$Oset" == "1" ];do

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
			xtick=`echo $bds | awk -F "/" '{s=1;print ($2-$1)/s/7}'`
			ytick=`echo $bds | awk -F "/" '{s=1;print ($4-$3)/s/7}'`
		else
			xtick=`echo $tick | awk -F "/" '{print $1}'`
			ytick=`echo $tick | awk -F "/" '{if (1==NF){print $1}else{print $2}}'`
		fi
		xtickf=`echo $xtick | awk '{print $1/2}'`
		ytickf=`echo $ytick | awk '{print $1/2}'`

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
				if [ "$Hset" == "" ]; then
					HEIGHT=`echo $bds | awk -F "/" '{printf("%fi\n",($4-$3)/($2-$1)*7)}'`
				else
					HEIGHT="7i"
				fi
				if [ "$rset" != "1" ]; then
					PROJ="X7i/"$HEIGHT
				else
					PROJ="X7i/"-$HEIGHT
				fi
			else
				HEIGHT="-"
			fi
			AXIS=-Bf${xtickf}a${xtick}:"":/f${ytickf}a${ytick}:""::."$title":WSne
		else
			# geographic coordinates
			HEIGHT=7i
			PROJ="M$HEIGHT"
		        AXIS=-B${xtick}:"":/${ytick}:""::."$title":WSne
		fi

		Jset="-J"
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
				xtick=`echo $bds | awk -F "/" '{s=1;print ($2-$1)/s/7}'`
				ytick=`echo $bds | awk -F "/" '{s=1;print ($4-$3)/s/7}'`
			fi
			xtickf=`echo $xtick | awk '{print $1/2}'`
			ytickf=`echo $ytick | awk '{print $1/2}'`

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
		if [ "$U3" != "" ]; then
			TRIANGLES=`grdinfo -C $U3 | awk -v m=$m 'function max(x,y){return (x>y)?x:y};function abs(x){return (0<x)?x:-x}{if (m<max(abs($6),abs($7))){print "-E"}}'`
		fi
	fi

	my_gmt $INDEX
	
	# add trailer information
	echo %% Postscript created with >> $PSFILE
	echo %% $(basename $0) $cmdline >> $PSFILE
	
	# open file for display
	echo $self": Created map "$PSFILE
	
	ps2pdf -sPAPERSIZE="archA" -dPDFSETTINGS=/prepress $PSFILE $PDFFILE
	echo $self": Converted to pdf file "$PDFFILE

	if [ "$xset" != "1" ]; then
		#display -trim $PSFILE &
		#gv -geometry +0+0 -spartan -scale=0.5 $PSFILE &
		xpdf -geometry +0+0 -paper "archA" $PDFFILE -z 100 -g 565x755 -z width >& /dev/null &
	fi
	
	if [ "$#" != "0" ];then
		shift
		if [ "$#" != "0" ];then
			echo ""
		fi
	fi

	# prevent more empty plots (cancel -E option)
	Oset=""
done
