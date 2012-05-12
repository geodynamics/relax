if not exist co-stress mkdir co-stress


set OMP_NUM_THREADS=4

set GMT_SHAREDIR=%CD%\..\..\share
set PROJ_LIB=%CD%\..\..\share
..\..\relax --no-proj-output < co-stress.input
