#!/bin/sh

set -e
self=$(basename $0)

my_list_last()
{
  if [ $# -eq 0 ]; then
    echo "";exit 0
fi

while [ $# -ne 1 ]; do
  shift
done
echo $1
return 0
}

my_list_first()
{
if [ $# -eq 0 ]; then
    echo "";exit 0
fi

echo $1
return 0
}

my_list_remove_head()
{
shift
echo $*
return 0
}

my_gmt()
{

psxy $OVERLAY -K -JX5i/1.0i -X$SHIFTX -Y$SHIFTY -R$bds1 -N \
	-Ba$tick:"x":/a$tick1:"u1"::.:WSne \
	-W0.5p/0/0/0 \
	<<EOF >> $PSFILE
`awk '{print $1,$2}' $U1`
EOF

psxy -O -K -JX -R$bds2 -N -Y1.5i \
	-Ba$tick::/a$tick2:"u2"::.:Wsne \
	-W0.5p/0/0/0 \
	<<EOF >> $PSFILE
`awk '{print $1,$2}' $U2`
EOF

psxy -O $MORETHING -JX -R$bds3 -N -Y1.5i \
	-Ba$tick::/a$tick3:"-u3"::."$TITLE":Wsne \
	-W0.5p/0/0/0 \
	<<EOF >> $PSFILE
`awk '{print $1,$2}' $U3`
EOF

}

if (test $# -lt "1"); then
        echo usage: $self [-i tick] [-x -5/5] [-u -0.1/0.1] [...] wdir/index 
	echo or
	echo usage: $self [-o plot.ps] [...] wdir/index1 ... wdir/indexN
	echo ""
        echo options:
	echo          -i tick interval
	echo          -o filename.ps stack all profiles in the same plot
	echo          -x interval sets the range on the x-axis
	echo          -u interval sets the range on the u1-axis
	echo          -v interval sets the range on the u2-axis
	echo          -w interval sets the range on the u3-axis
	echo ""
	echo "plots the profiles index1-up.xy index1-north.xy index1-east.xy"
        
        exit
fi

gmtset LABEL_FONT_SIZE 12p
gmtset HEADER_FONT_SIZE 12p
gmtset ANOT_FONT_SIZE 12p 
gmtset COLOR_BACKGROUND 0/0/255
gmtset COLOR_FOREGROUND 255/0/0
gmtset PAPER_MEDIA a5
gmtset HEADER_OFFSET 0.c

while getopts ":x:i:o:u:v:w:" flag
do
  #echo "$flag" $OPTIND $OPTARG
  case "$flag" in
    i) iset=1;tick=$OPTARG;;
    o) oset=1;TITLE=$(basename $OPTARG);PSFILE=$(dirname $OPTARG)/$(basename $OPTARG .ps).ps;;
    x) xset=1;START=`echo $OPTARG | awk -F"/" '{print $1}'`;END=`echo $OPTARG | awk -F"/" '{print $2}'`;;
    u) uset=1;MIN1=`echo $OPTARG | awk -F"/" '{print $1}'`;MAX1=`echo $OPTARG | awk -F"/" '{print $2}'`;;
    v) vset=1;MIN2=`echo $OPTARG | awk -F"/" '{print $1}'`;MAX2=`echo $OPTARG | awk -F"/" '{print $2}'`;;
    w) wset=1;MIN3=`echo $OPTARG | awk -F"/" '{print $1}'`;MAX3=`echo $OPTARG | awk -F"/" '{print $2}'`;;
  esac
done

for foo in $xset $iset $oset $uset $vset $wset;do
	shift;shift;
done

if [ "$oset" == "1" ]; then
	PDFFILE=$(dirname $PSFILE)/$(basename $PSFILE .ps).pdf
	rm -f $PSFILE $PDFFILE
fi

# initial values
MORETHING="-K"
SHIFTX=2i
SHIFTY=1i

while [ "$#" != "0" ];do

if [ "$oset" == "1" ]; then
	if [ "$#" == "1" ]; then
		MORETHING=""
	fi
fi

WDIR=`dirname $1`
INDEX=`basename $1 .xy`
U1=$WDIR/$INDEX-north.xy
U2=$WDIR/$INDEX-east.xy
U3=$WDIR/$INDEX-up.xy

echo $self": in directory "$WDIR", ploting profile "$INDEX

if [ ! -e $U3 ]; then
	echo $self": file "$U3" does not exist. exiting."
        exit
fi
if [ "$xset" != "1" ]; then
	START=`head -n 1 $U3 | awk '{print $1}'`
	END=`tail -n 1 $U3 | awk '{print $1}'`
fi
if [ "$iset" != "1" ]; then
	tick=`echo $START $END | awk '{printf "%.0e",($2-$1)/5}'`
fi
if [ "$uset" != "1" ]; then
	MAX1=`awk -v x=$END 'BEGIN{m=-1e20}{if($2>m && $1<=x){m=$2}}END{print m}' $U1`
	MIN1=`awk -v x=$END 'BEGIN{m=+1e20}{if($2<m && $1<=x){m=$2}}END{print m}' $U1`
fi
if [ "$vset" != "1" ]; then
	MAX2=`awk -v x=$END 'BEGIN{m=-1e20}{if($2>m && $1<=x){m=$2}}END{print m}' $U2`
	MIN2=`awk -v x=$END 'BEGIN{m=+1e20}{if($2<m && $1<=x){m=$2}}END{print m}' $U2`
fi
if [ "$wset" != "1" ]; then
	MAX3=`awk -v x=$END 'BEGIN{m=-1e20}{if($2>m && $1<=x){m=$2}}END{print m}' $U3`
	MIN3=`awk -v x=$END 'BEGIN{m=+1e20}{if($2<m && $1<=x){m=$2}}END{print m}' $U3`
fi
if [ "$oset" != "1" ]; then
	TITLE=$WDIR/$INDEX
	PSFILE=$WDIR/$INDEX-plot.ps
	PDFFILE=$WDIR/$INDEX-plot.pdf
	rm -f $PSFILE $PDFFILE
fi

bds1=$START/$END/$MIN1/$MAX1
bds2=$START/$END/$MIN2/$MAX2
bds3=$START/$END/$MIN3/$MAX3
echo $self: using u1 range $MIN1/$MAX1, u2 range $MIN2/$MAX2, u3 range $MIN3/$MAX3
tick1=`echo $MIN1 $MAX1 | awk '{printf "%.0e", ($2-$1)/5}'`
tick2=`echo $MIN2 $MAX2 | awk '{printf "%.0e", ($2-$1)/5}'`
tick3=`echo $MIN3 $MAX3 | awk '{printf "%.0e", ($2-$1)/5}'`

echo $self": using x-axis interval ("$START","$END") with tick "$tick

my_gmt

if [ "$oset" == "1" ]; then
	OVERLAY=-O
	SHIFTX=0i
	SHIFTY=-3i
else
	echo $self": Created map "$PSFILE
	#gv $PSFILE &
	ps2pdf -sPAPERSIZE="a4" $PSFILE $PDFFILE
	echo $self": Converted to pdf file "$PDFFILE
	xpdf -paper "A5" $PDFFILE -z -0 >& /dev/null &
fi

shift
echo ""
done

if [ "$oset" == "1" ]; then
	echo $self": Created map "$PSFILE
	#gv $PSFILE &
	ps2pdf -sPAPERSIZE="a4" $PSFILE $PDFFILE
	echo $self": Converted to pdf file "$PDFFILE
	xpdf -paper "A5" $PDFFILE -z -0 >& /dev/null &
fi

