#!/bin/csh
echo       'regression/compile'
make -j10  USE_MPI=FALSE;
echo       'regression/aggpwlfp 1 proc'
 ./aggpwlfpTest2d.gnu.DEBUG.ex aggpwlfp.inputs | grep -i test;
echo       'regression/coarse ave 1 proc'
 ./ebCoarseAveTest2d.gnu.DEBUG.ex ebcoarave.inputs | grep -i test;
echo       'regression/flux reg 1 proc'
./fluxRegTest2d.gnu.DEBUG.ex fluxreg.inputs | grep -i test;
echo       'regression/nwoebquadvfi 1 proc'
./nwoEBQuadCFITest2d.gnu.DEBUG.ex nwoebquadcfi.inputs | grep -i test;
echo       'regression/regfluxregtest 1 proc'
./regFluxRegTest2d.gnu.DEBUG.ex fluxreg.inputs | grep -i test;
echo       'regression/meshrefine 1 proc'
./simpleMeshRefine2d.gnu.DEBUG.ex meshref.inputs | grep -i test;
echo       'regression/fine interp 1 proc'
./ebFineInterpTest2d.gnu.DEBUG.ex ebfineinterp.inputs | grep -i test;
exit
