if not exist mogi mkdir mogi


set OMP_NUM_THREADS=4

set GMT_SHAREDIR=%CD%\..\..\share
set PROJ_LIB=%CD%\..\..\share
..\..\relax < mogi.input
