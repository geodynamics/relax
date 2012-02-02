set relax=`pwd`

if ( ! -f relax) then
    echo
    echo "*** Error! ***"
    echo
    echo "Source this script from the top-level Relax directory:"
    echo
    echo "    cd [directory containing 'setup.sh']"
    echo "    source setup.sh"
    echo 
else
    setenv PATH "{$relax}:{$relax/util}:{$PATH}"
    setenv GMT_SHAREDIR "$relax/share"
    setenv PROJ_LIB "$relax/share"
    echo "Ready to run Relax."
endif
