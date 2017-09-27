#!/bin/csh
echo       'regression/compile'
make -j10  USE_MPI=TRUE;
echo       'regression/aggpwlfp 1 proc'
mpirun -np 1 ./aggpwlfpTest2d.gnu.DEBUG.MPI.ex aggpwlfp.inputs | grep -i test;
echo       'regression/aggpwlfp 2 proc'
mpirun -np 2 ./aggpwlfpTest2d.gnu.DEBUG.MPI.ex aggpwlfp.inputs | grep -i test;
echo       'regression/coarse ave 1 proc'
mpirun -np 1 ebCoarseAveTest2d.gnu.DEBUG.MPI.ex ebcoarave.inputs | grep -i test;
echo       'regression/coarse ave 2 proc'
mpirun -np 2 ebCoarseAveTest2d.gnu.DEBUG.MPI.ex ebcoarave.inputs | grep -i test;
echo       'regression/flux reg 1 proc'
mpirun -np 1 fluxRegTest2d.gnu.DEBUG.MPI.ex fluxreg.inputs | grep -i test;
echo       'regression/flux reg 2 proc'
mpirun -np 2 fluxRegTest2d.gnu.DEBUG.MPI.ex fluxreg.inputs | grep -i test;
echo       'regression/nwoebquadvfi 1 proc'
mpirun -np 1 nwoEBQuadCFITest2d.gnu.DEBUG.MPI.ex nwoebquadcfi.inputs | grep -i test;
echo       'regression/nwoebquadvfi 2 proc'
mpirun -np 2 nwoEBQuadCFITest2d.gnu.DEBUG.MPI.ex nwoebquadcfi.inputs | grep -i test;
echo       'regression/regfluxregtest 1 proc'
mpirun -np 1 regFluxRegTest2d.gnu.DEBUG.MPI.ex fluxreg.inputs | grep -i test;
echo       'regression/regfluxregtest 2 proc'
mpirun -np 2 regFluxRegTest2d.gnu.DEBUG.MPI.ex fluxreg.inputs | grep -i test;
echo       'regression/meshrefine 1 proc'
mpirun -np 1 simpleMeshRefine2d.gnu.DEBUG.MPI.ex meshref.inputs | grep -i test;
echo       'regression/meshrefine 2 proc'
mpirun -np 2 simpleMeshRefine2d.gnu.DEBUG.MPI.ex meshref.inputs | grep -i test;
echo       'regression/fine interp 1 proc'
mpirun -np 1 ebFineInterpTest2d.gnu.DEBUG.MPI.ex ebfineinterp.inputs | grep -i test;
echo       'regression/fine interp 2 proc'
mpirun -np 2 ebFineInterpTest2d.gnu.DEBUG.MPI.ex ebfineinterp.inputs | grep -i test;
exit
