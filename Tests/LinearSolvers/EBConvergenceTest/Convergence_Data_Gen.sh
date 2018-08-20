# Bash Script for generating scaling data for EB convergence Test
# Written by Steven Reeves, August 6th 2018
# Center for Computational Science and Engineering
# Lawrence Berkeley National Laboratory

# Runs the EB Elliptic Test Solver for n_cell = 32 to 1024 for 2D 
# and n_cell = 16 to 256 for 3D. 
# Then converts the multifabs into a .mat file 
# After all processing it moves the mat files into the Results folder.
# Within the Results folder there is an octave file for plotting the results 

# Make sure you build the Multifab to matlab tools and 
# export PATH="$PATH:/../../../Tools/Postprocessing/C_Src"

if [ $1 -eq 2 ];
then
echo "2D Test Commencing!" 
./main2d.gnu.TEST.MPI.ex inputs n_cell=32 max_grid_size=32
mv phi-0_H phi-0_32_H
MultiFabToMatLab2d.gnu.MPI.ex infile=phi-0_32

./main2d.gnu.TEST.MPI.ex inputs n_cell=64 max_grid_size=64
mv phi-0_H phi-0_64_H
MultiFabToMatLab2d.gnu.MPI.ex infile=phi-0_64

./main2d.gnu.TEST.MPI.ex inputs n_cell=128 max_grid_size=128
mv phi-0_H phi-0_128_H
MultiFabToMatLab2d.gnu.MPI.ex infile=phi-0_128

./main2d.gnu.TEST.MPI.ex inputs n_cell=256 max_grid_size=256
mv phi-0_H phi-0_256_H
MultiFabToMatLab2d.gnu.MPI.ex infile=phi-0_256

./main2d.gnu.TEST.MPI.ex inputs n_cell=512 max_grid_size=512
mv phi-0_H phi-0_512_H
MultiFabToMatLab2d.gnu.MPI.ex infile=phi-0_512

./main2d.gnu.TEST.MPI.ex inputs n_cell=1024 max_grid_size=1024
mv phi-0_H phi-0_1024_H
MultiFabToMatLab2d.gnu.MPI.ex infile=phi-0_1024
elif [ $1 -eq 3 ]; 
then 
echo "Commencing 3D Test!"
./main3d.gnu.TEST.MPI.ex inputs n_cell=16 max_grid_size=16
mv phi-0_H phi-0_16_H
MultiFabToMatLab3d.gnu.MPI.ex infile=phi-0_16

./main3d.gnu.TEST.MPI.ex inputs n_cell=32 max_grid_size=32
mv phi-0_H phi-0_32_H
MultiFabToMatLab3d.gnu.MPI.ex infile=phi-0_32

./main3d.gnu.TEST.MPI.ex inputs n_cell=64 max_grid_size=64
mv phi-0_H phi-0_64_H
MultiFabToMatLab3d.gnu.MPI.ex infile=phi-0_64

./main3d.gnu.TEST.MPI.ex inputs n_cell=128 max_grid_size=128
mv phi-0_H phi-0_128_H
MultiFabToMatLab3d.gnu.MPI.ex infile=phi-0_128

./main3d.gnu.TEST.MPI.ex inputs n_cell=256 max_grid_size=256
mv phi-0_H phi-0_256_H
MultiFabToMatLab3d.gnu.MPI.ex infile=phi-0_256
fi 
mv *.mat Results/
rm phi-0* vfrc-0*
