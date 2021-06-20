if not exist relax.exe (goto :error)

SET GMT_SHAREDIR=%CD%\share


PATH=%CD%;%PATH%

echo "Ready to run Relax"



goto :theend

:error
echo

echo *** Error! ***

echo

echo Run this script from the top-level Relax directory:

echo

echo     cd [directory containing 'setup_win.bat']

echo     setup_win.bat

echo


:theend
