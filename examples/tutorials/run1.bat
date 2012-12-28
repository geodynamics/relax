if not exist output1 mkdir output1


set OMP_NUM_THREADS=4

set GMT_SHAREDIR=%CD%\..\..\share
set PROJ_LIB=%CD%\..\..\share
..\..\relax < run1.input
PAUSE
