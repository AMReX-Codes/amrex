./mkmf/bin/mkmf -t template.mk -p main3d.gnu.MPI.ex
make
rm -rf plt*
mpirun -np 8 ./main3d.gnu.MPI.ex inputs
