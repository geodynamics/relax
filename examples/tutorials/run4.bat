if not exist output4 mkdir output4


set OMP_NUM_THREADS=4

set GMT_SHAREDIR=%CD%\..\..\share
set PROJ_LIB=%CD%\..\..\share
..\..\relax --no-stress-output --no-proj-output < run4.input
