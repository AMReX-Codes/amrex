#!/bin/csh
set thisdir = `pwd`
echo       'ebgraphDistributed'
cd $thisdir/ebgraphDistributed; make -j10 ; ./ebgraphDist2d.gnu.DEBUG.ex ./sphere.inputs | grep test;
echo       'ebgraphSingleGrid'
cd $thisdir/ebgraphSingleGrid; make -j10 ; ./ebgraphSG2d.gnu.DEBUG.ex ./sphere.inputs | grep test;
echo       'flatPlate'
cd $thisdir/flatPlate; make -j10 ; ./flatPlateTest2d.gnu.DEBUG.ex ./flatplate.inputs | grep test;
echo       'regression/compile'
cd $thisdir/regression/; make -j10 ; 
echo       'regression/dataArith'
cd $thisdir/regression/; ./dataArith2d.gnu.DEBUG.ex dataarith.inputs | grep test;
echo       'regression/ebio'
cd $thisdir/regression/; ./ebio2d.gnu.DEBUG.ex ebio.inputs | grep test;
echo       'regression/ebnormalizeTest'
cd $thisdir/regression/; ./ebnormalizeTest2d.gnu.DEBUG.ex ebnormtest.inputs | grep test;
echo       'regression/fabio'
cd $thisdir/regression/; ./fabio2d.gnu.DEBUG.ex ebio.inputs | grep test;
echo       'regression/levelRedist'
cd $thisdir/regression/; ./levelRedistTest2d.gnu.DEBUG.ex levelredist.inputs | grep test;
echo       'regression/serialization'
cd $thisdir/regression/; ./serialization2d.gnu.DEBUG.ex serialization.inputs | grep test;
echo       'sparseDataDistributed'
cd $thisdir/sparseDataDistributed; make -j10 ; ./sparseDataDist2d.gnu.DEBUG.ex ./sphere.inputs | grep test;
echo       'sparseDataSingleGrid'
cd $thisdir/sparseDataSingleGrid; make -j10 ; ./sparseDataSG2d.gnu.DEBUG.ex ./sphere.inputs | grep test;
exit
