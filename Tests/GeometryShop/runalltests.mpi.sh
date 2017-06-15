#!/bin/csh
set thisdir = `pwd`
echo       'ebgraphDistributed 1 proc'
cd $thisdir/ebgraphDistributed; make -j10 all USE_MPI=TRUE;mpirun -np 1 ./ebgraphDist2d.gnu.DEBUG.MPI.ex ./sphere.inputs | grep test;
echo       'ebgraphDistributed 2 proc'
cd $thisdir/ebgraphDistributed; make -j10 all USE_MPI=TRUE;mpirun -np 2 ./ebgraphDist2d.gnu.DEBUG.MPI.ex ./sphere.inputs | grep test;
echo       'regression/compile'
cd $thisdir/regression/; make -j10 all USE_MPI=TRUE;
echo       'regression/dataArith 1 proc'
cd $thisdir/regression/;mpirun -np 1 ./dataArith2d.gnu.DEBUG.MPI.ex dataarith.inputs | grep test;
echo       'regression/dataArith 2 proc'
cd $thisdir/regression/;mpirun -np 2 ./dataArith2d.gnu.DEBUG.MPI.ex dataarith.inputs | grep test;
echo       'regression/ebio 1 proc'
cd $thisdir/regression/;mpirun -np 1 ./ebio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/ebio 2 proc'
cd $thisdir/regression/;mpirun -np 2 ./ebio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/ebnormalizeTest 1 proc'
cd $thisdir/regression/;mpirun -np 1 ./ebnormalizeTest2d.gnu.DEBUG.MPI.ex ebnormtest.inputs | grep test;
echo       'regression/ebnormalizeTest 2 proc'
cd $thisdir/regression/;mpirun -np 2 ./ebnormalizeTest2d.gnu.DEBUG.MPI.ex ebnormtest.inputs | grep test;
echo       'regression/fabio 1 proc'
cd $thisdir/regression/;mpirun -np 1 ./fabio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/fabio 2 proc'
cd $thisdir/regression/;mpirun -np 2 ./fabio2d.gnu.DEBUG.MPI.ex ebio.inputs | grep test;
echo       'regression/levelRedist 1 proc'
cd $thisdir/regression/;mpirun -np 1 ./levelRedistTest2d.gnu.DEBUG.MPI.ex levelredist.inputs | grep test;
echo       'regression/levelRedist 2 proc'
cd $thisdir/regression/;mpirun -np 2 ./levelRedistTest2d.gnu.DEBUG.MPI.ex levelredist.inputs | grep test;
echo       'sparseDataDistributed 1 proc'
cd $thisdir/sparseDataDistributed; make -j10 all USE_MPI=TRUE;mpirun -np 1 ./sparseDataDist2d.gnu.DEBUG.MPI.ex ./sphere.inputs | grep test;
echo       'sparseDataDistributed 2 proc'
cd $thisdir/sparseDataDistributed; make -j10 all USE_MPI=TRUE;mpirun -np 2 ./sparseDataDist2d.gnu.DEBUG.MPI.ex ./sphere.inputs | grep test;
exit
