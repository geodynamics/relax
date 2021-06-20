if not exist landers mkdir landers


set OMP_NUM_THREADS=4

set GMT_SHAREDIR=%CD%\..\..\share
set PROJ_LIB=%CD%\..\..\share
..\..\relax --no-vtk-output --no-proj-output --no-stress-output < landers.input
PAUSE
