#!/bin/bash

set -e
self=$(basename $0)
trap 'echo $self: Some errors occurred. Exiting.; exit' ERR

my_list_last()
{
if [ $# -eq 0 ]; then
    echo "";exit 0
fi

while [ $# -ne 1 ]; do
  shift
done
echo $1
return 1
}

my_list_first()
{
if [ $# -eq 0 ]; then
    echo "";exit 0
fi

echo $1
return 1
}

my_list_remove_head()
{
shift
echo $*
return 1
}

my_list_count()
{
count=0
while [ $# -ne 0 ]; do
  shift
  count=`echo "$count + 1"|bc`
done
echo $count
return 1
}

if (test $# -lt "2"); then
	echo "usage: $self wdir/file.disp index"
        echo "creates wdir/index-{north,east,up}.xyz from wdir/file.disp"
	exit
fi

WDIR=`dirname $1`
DFILE=$WDIR/`basename $1`
INDEX=$2

echo $self": Using directory "$WDIR
echo $self": Processing file "$DFILE

XYZFILE1=$WDIR/$INDEX-east.xyz
XYZFILE2=$WDIR/$INDEX-north.xyz
XYZFILE3=$WDIR/$INDEX-up.xyz

echo $self": Creating files " $INDEX-east.xyz $INDEX-north.xyz $INDEX-up.xyz

awk '{if (FNR > 3){print $2,$1,$4}}' $DFILE > $XYZFILE1 
awk '{if (FNR > 3){print $2,$1,$3}}' $DFILE > $XYZFILE2
awk '{if (FNR > 3){print $2,$1,-$5}}' $DFILE > $XYZFILE3 

echo $self": Converted all files."




