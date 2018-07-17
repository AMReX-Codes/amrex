#!/bin/csh
echo       'regression/compile'
make -j10 USE_MPI=FALSE;
echo       'regression/dataArith 1 proc'
./dataArith2d.gnu.DEBUG.ex dataarith.inputs | grep test;
echo       'regression/ebio 1 proc'
ebio2d.gnu.DEBUG.ex ebio.inputs | grep test;
echo       'regression/ebnormalizeTest 1 proc'
./ebnormalizeTest2d.gnu.DEBUG.ex ebnormtest.inputs | grep test;
echo       'regression/fabio 1 proc'
./fabio2d.gnu.DEBUG.ex fabio.inputs | grep test;
echo       'regression/fabfromif 1 proc'
./fabfromif2d.gnu.DEBUG.ex fabfromif.inputs | grep test;
echo       'regression/levelRedist 1 proc'
./levelRedistTest2d.gnu.DEBUG.ex levelredist.inputs | grep test;
echo       'regression/multicelllev 1 proc'
./multicelllev2d.gnu.DEBUG.ex multicelllev.inputs | grep test;
echo       'regression/multicelllev 1 proc'
./multicelllev2d.gnu.DEBUG.ex multicelllev.inputs | grep test;

exit
