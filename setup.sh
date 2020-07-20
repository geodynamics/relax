relax=`pwd`

if test ! -f relax; then
    echo
    echo "*** Error! ***"
    echo
    echo "Source this script from the top-level Relax directory:"
    echo
    echo "    cd [directory containing 'setup.sh']"
    echo "    source setup.sh"
    echo 
else
    export PATH="$relax:$relax/util:$PATH"
    export GMT_SHAREDIR="$relax/share"
    export PROJ_LIB="$relax/share"
    echo "Ready to run Relax."
fi
