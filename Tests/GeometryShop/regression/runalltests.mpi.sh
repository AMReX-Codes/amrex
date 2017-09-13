#!/bin/csh
echo       'regression/compile'
make -j10 USE_MPI=TRUE;
echo       'regression/dataArith 1 proc'
mpirun -np 1 ./dataArith2d.gnu.DEBUG.MPI.ex dataarith.inputs | grep test;
echo       'regression/dataArith 2 proc'
mpirun -np 2 ./dataArith2d.gnu.DEBUG.MPI.ex dataarith.inputs | grep test;
echo       'regression/ebio 1 proc'
mpirun -np 1 ./ebio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/ebio 2 proc'
mpirun -np 2 ./ebio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/ebio 2 proc (bigger prob)'
mpirun -np 2 ./ebio2d.gnu.DEBUG.MPI.ex ebio.2.inputs | grep test;
echo       'regression/ebnormalizeTest 1 proc'
mpirun -np 1 ./ebnormalizeTest2d.gnu.DEBUG.MPI.ex ebnormtest.inputs | grep test;
echo       'regression/ebnormalizeTest 2 proc'
mpirun -np 2 ./ebnormalizeTest2d.gnu.DEBUG.MPI.ex ebnormtest.inputs | grep test;
echo       'regression/fabio 1 proc'
mpirun -np 1 ./fabio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/fabio 2 proc'
mpirun -np 2 ./fabio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/levelRedist 1 proc'
mpirun -np 1 ./levelRedistTest2d.gnu.DEBUG.MPI.ex levelredist.inputs | grep test;
echo       'regression/levelRedist 2 proc'
mpirun -np 2 ./levelRedistTest2d.gnu.DEBUG.MPI.ex levelredist.inputs | grep test;

exit
