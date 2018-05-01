ls ../Exec/SingleVortex/*.f90 Src_3d/*.f90 > path_names
./mkmf/bin/mkmf -t template.mk -p main3d.MPI.gnu.ex path_names
make
rm -rf plt*
mpirun -np 4 ./main3d.MPI.gnu.ex ../Exec/SingleVortex/inputs

