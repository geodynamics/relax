#!/bin/bash

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# script grdview.sh
# display a view of GMT .grd files with
# overlay from user-provided GMT scripts.
#
# Sylvain Barbot (06/13/07) - original form
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

set -e
self=$(basename $0)
selfdir=$(dirname $0)
cmdline=$*
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

usage(){
        echo "usage: $self [-b -5/5/-5/5/-5/5] [-E 180/90] [-p -1.5/1.5/0.1] [...] file1.grd ... fileN.grd"
	echo "or"
        echo "       $self [-b -5/5/-5/5/-5/5] [-E 180/90] [-p -1.5/1.5/0.1] [...] dir1/index1 ... dirN/indexN"
	echo ""
        echo "options:"
        echo "         -b bds sets the map bound to bds"
	echo "         -c palette_name [default my_jet]"
        echo "         -e efile.sh runs efile.sh to add content to map"
        echo "         -g switch to geographic projection (longitude/latitude)"
        echo "         -h display this help message and exit"
	echo "         -i file.grd illuminates the grd image with file.grd"
        echo "         -o output directory"
        echo "         -p palette sets the color scale bounds to palette"
        echo "         -r reverse the y-axis"
        echo "         -s step sets the distance between vectors"
        echo "         -u defines the color scale unit"
        echo "         -v vertical exxageration"
	echo "         -t tick interval"
	echo "         -T title header"
	echo "         -x do not display map (only create .ps file)"
	echo "         -E Sets the viewpoint's azimuth and elevation for perspective view [180/90]"
	echo "         -G drapefile"
	echo "         -N base level"
	echo "         -O file.ps only plot base map with extra scripts"
	echo "         -J overwrites the geographic projections (for -J o the -b option is relative)"
        echo "         -S size of the plot"
        echo "         -X shift the plot horizontally on the page"
        echo "         -Y shift the plot vertically on the page"
	echo ""
	echo "Creates N maps based on grd files file1.grd ... fileN.grd, or based"
	echo "on files dir1/index1-{up,north.east}.grd ... dirN/indexN-{up,north,east}.grd"
	echo ""
        
        exit
}

my_gmt(){

	if [ -e "$U3" ]; then
		echo $self: grdview $U3 -Qi300 -R$bds $PROJ $PROJZ $Eset $Nset \
		    ${AXIS}SWNEZ+ $Eset $Gset -C$colors -P -X${XSHIFT}i -Y${YSHIFT}i $iset
		grdview $U3 -Qi300 -R$bds $PROJ $PROJZ $Eset $Nset \
		    ${AXIS}SWNEZ+ $Eset $Gset \
		    -K -C$colors -P -X1.2i -Y${YSHIFT}i $iset \
		    > $PSFILE
	else
		echo $self: psbasemap -R$bds $PROJ $PROJZ $AXIS $Eset
		psbasemap -R$bds $PROJ $PROJZ \
		    ${AXIS}wsne $Eset \
		     -K -P -X${XSHIFT}i -Y${YSHIFT}i \
		    > $PSFILE
	fi

	# add extra layers
	for subprog in $EXTRA; do
		if [ "" != "$U3" ]; then
			OPTIONP="-p $U3"
		else
			OPTIONP=""
		fi
		if [ -e "$selfdir/$subprog" ]; then
			#echo $self: running $subprog $PSFILE $bds $VECTOR $U3 $HEIGHT
			eval $subprog $gset -b $bds $OPTIONP -H $HEIGHT -E ""$Eset"" -J ""$PROJ"" $PSFILE
		else
			if [ -e "$subprog" ]; then
				EsetP='\"$Eset\"'
				PROJP='\"$PROJ\"'
				PROJZP='\"$PROJZ\"'
				#echo $self: running $subprog $gset -b $bds $OPTIONP -H $HEIGHT -E $EsetP -J $PROJP -Z $PROJZP $PSFILE
				eval $subprog $gset -b $bds $OPTIONP -H $HEIGHT -E $EsetP -J $PROJP -Z $PROJZP $PSFILE
			fi
		fi
	done

	if [ -e "$colors" ]; then
		#-Q0.20c/1.0c/0.4cn1.0c \
		psscale -O -K -B$PSSCALE/:$unit: -D3.5i/-0.8i/7.1i/0.2ih \
			$TRIANGLES -C$colors $REVERT \
			>> $PSFILE 
		
		rm -f $colors
	fi

	if [ ! -e "$U3" ]; then
		echo $self: closing frame
		psbasemap -R$bds $PROJ $PROJZ \
			${AXIS}WSNEZ+ $Eset -O -K -P >> $PSFILE
	fi

	#echo 0 0 | psxy -O $PROJ -R$bds -Sp0.001c >> $PSFILE

}

gmtset LABEL_FONT_SIZE 12p
gmtset HEADER_FONT_SIZE 12p
gmtset ANOT_FONT_SIZE 12p 
gmtset COLOR_BACKGROUND 0/0/255
gmtset COLOR_FOREGROUND 255/0/0
gmtset COLOR_NAN 255/255/255
#gmtset PAPER_MEDIA A0
gmtset PAPER_MEDIA archA

libdir=$(dirname $0)/../share
EXTRA=""

while getopts "b:c:e:ghi:o:p:v:s:t:T:u:xrE:G:HN:O:J:S:X:Y:" flag; do
	case "$flag" in
	b) bset=1;bds=$OPTARG;;
	c) cset="-c";carg=$OPTARG;;
	e) eset=1;EXTRA="$EXTRA $OPTARG";;
	g) gset="-g";;
	h) hset=1;;
	i) iset="-I$OPTARG";;
	o) oset=1;ODIR=$OPTARG;;
	p) pset=1;P=$OPTARG;PSSCALE=`echo $OPTARG | awk -F"/" 'function abs(x){return x<0?-x:x}{print abs($2-$1)/6}'`;;
	r) rset=1;;
	s) sset=1;ADX=$OPTARG;;
	t) tset=1;tickx=$OPTARG;ticky=$OPTARG;tickz=$OPTARG;;
	T) Tset=1;title=$OPTARG;;
	u) uset=1;unit=$OPTARG;;
	v) vset=1;verticalexxageration=$OPTARG;;
	x) xset=1;;
	E) Eset="-E$OPTARG";;
	G) Gset="-G$OPTARG";;
	H) Hset="-H";;
	N) Nset="-N$OPTARG";;
	O) Oset=1;PSFILE=$(dirname $OPTARG)/$(basename $OPTARG .ps).ps;;
	J) Jset=1;PROJ="-J$OPTARG";;
	S) Sset=1;SIZE=$OPTARG;;
	X) Xset=1;Xshift=$OPTARG;;
	Y) Yset=1;Yshift=$OPTARG;;
	esac
done
for item in $bset $cset $iset $oset $pset $sset $tset $vset $Tset $uset $Sset $Xset $Yset $Eset $Gset $Nset $Oset $Jset $EXTRA;do
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
	verticalexxageration=1
fi

# size of plot
if [ "$Sset" != "1" ]; then
	SIZE=5
fi

# horizontal shift of plot
if [ "$Xset" != "1" ]; then
	XSHIFT=1.2
else
	XSHIFT=`echo $Xshift | awk '{print $1+1.2}'`
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

# azimuth and elevation
if [ "$Eset" == "" ]; then
	Eset="-E180/90"
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
		U3=$WDIR/$INDEX-up.grd

		if [ -e "$U3" ]; then
			U3=$U3;
		else
			U3=$WDIR/$INDEX
		fi
		if [ "$bset" != "1" ]; then
			bds=`grdinfo -C $U3 | awk '{s=1;print $2/s"/"$3/s"/"$4/s"/"$5/s"/"$6/s"/"$7/s}'`
		fi

		# tick marks
		if [ "$tset" != "1" ]; then
			tickx=`echo $bds | awk -F "/" '{s=1;print ($2-$1)/s/4}'`
			ticky=`echo $bds | awk -F "/" '{s=1;print ($4-$3)/s/4}'`
			tickz=`echo $bds | awk -F "/" '{s=1;print ($6-$5)/s/4}'`
		fi
		tickxf=`echo $tickx | awk '{print $1/2}'`
		tickyf=`echo $ticky | awk '{print $1/2}'`
		tickzf=`echo $tickz | awk '{print $1/2}'`

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
		else
			TRIANGLES=`grdinfo -C $U3 | \
					awk -v m=$m 'function max(x,y){return (x>y)?x:y}; \
						     function abs(x){return (0<x)?x:-x}{if (m<max(abs($6),abs($7))){print "-E"}}'`
		fi

		if [ "$gset" != "-g" ]; then
			if [ "$Jset" == "" ]; then
				# Cartesian coordinates
				if [ "$Hset" == "" ]; then
					HEIGHT=`echo $bds | awk -F "/" '{printf("%fi\n",($4-$3)/($2-$1)*5)}'`
				else
					HEIGHT="8i"
				fi
				SCALE=`echo $bds | awk -F "/" -v s=$verticalexxageration '{printf("%f\n",5/($2-$1)*s)}'`
				if [ "$rset" != "1" ]; then
					PROJ="-JX5i/"$HEIGHT
					PROJZ=-Jz$SCALE
				else
					PROJ="-JX5i/"-$HEIGHT
					PROJZ=-Jz$SCALE
				fi
			else
				HEIGHT="-"
			fi
			AXIS=-Bf${tickxf}a${tickx}:"x":/f${tickyf}a${ticky}:"y":/f${tickzf}a${tickz}:"z"::."$title":
		else
			# geographic coordinates
			HEIGHT=5i
			PROJ="-JM$HEIGHT"
		        AXIS=-B${tickx}:"":/${ticky}:""::."$title":WSne
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
			echo "$self: -b option should be set with the -O option. exiting."
			echo ""
			usage
		else

			# tick marks
			if [ "$tset" != "1" ]; then
				tickx=`echo $bds | awk -F "/" '{s=1;print ($2-$1)/s/5}'`
				ticky=`echo $bds | awk -F "/" '{s=1;print ($4-$3)/s/5}'`
				tickz=`echo $bds | awk -F "/" '{s=1;print ($6-$5)/s/5}'`
			fi
			tickxf=`echo $tickx | awk '{print $1/2}'`
			tickyf=`echo $ticky | awk '{print $1/2}'`
			tickzf=`echo $tickz | awk '{print $1/2}'`

			# Cartesian vs geographic coordinates
			if [ "$gset" != "-g" ]; then
				if [ "$Hset" == "" ]; then
					HEIGHT=`echo $bds | awk -F "/" -v s=$SIZE '{printf("%fi\n",($4-$3)/($2-$1)*s)}'`
				else
					HEIGHT="7i"
				fi
				SCALE=`echo $bds | awk -F "/" -v v=$verticalexxageration -v s=$SIZE '{printf("%f\n",s/($2-$1)*v)}'`
				if [ "$rset" != "1" ]; then
					PROJ=-JX${SIZE}i/$HEIGHT
					PROJZ=-Jz$SCALE
				else
					PROJ=-JX${SIZE}i/-$HEIGHT
					PROJZ=-Jz$SCALE
				fi
			        AXIS=-Ba${tickx}:"x":/a${ticky}:"y":/f${tickzf}a${tickz}:"z"::."$TITLE":
			else
				HEIGHT=${SIZE}i
				PROJ="-JM0/0/$HEIGHT"
			        AXIS=-B${tickx}:"":/${ticky}:"":/f50:""::."$TITLE":WSne
			fi
		fi
	fi

	# color bounds
	if [ "$pset" == "1" ]; then
		makecpt -C$cptfile -T$P -D > $colors;
		m=`echo $P | awk -F"/" \
			'function max(x,y){return (x>y)?x:y};function abs(x){return (0<x)?x:-x}{print max(abs($1),$2)}'`
	fi

	my_gmt $INDEX
	
	# add trailer information
	echo %% Postscript created with >> $PSFILE
	echo %% $(basename $0) $cmdline >> $PSFILE
	
	# open file for display
	echo $self": Created map "$PSFILE
	
	if [ "$xset" != "1" ]; then
		#display -trim $PSFILE &
		#gv -geometry +0+0 -spartan -scale=0.5 $PSFILE &
		ps2pdf -sPAPERSIZE="archA" -dPDFSETTINGS=/prepress $PSFILE $PDFFILE
		echo $self": Converted to pdf file "$PDFFILE
		xpdf -geometry +0+0 -paper "archA" $PDFFILE -z 100 -g 565x755 -z width >& /dev/null &
	fi
	
	if [ "$#" != "0" ];then
		shift
		if [ "$#" != "0" ];then
			echo ""
		fi
	fi

	# prevent more empty plots (cancel -O option)
	Oset=""
done
