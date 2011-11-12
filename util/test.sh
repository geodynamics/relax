#!/bin/sh

echo $*
for i in $* ; do

	echo $1, $i
	shift
done

