#!/bin/bash

mpirun -n 4 ./AmrDeriveSpectrum3d.gnu.MPI.ex input_spectrum3d
paste -d ' ' $1/x_vel*_spectrum.dat $1/y_vel*_spectrum.dat $1/z_vel*_spectrum.dat > $1/vel_spectrum.dat
paste -d ' ' $1/x_vort*_spectrum.dat $1/y_vort*_spectrum.dat $1/z_vort*_spectrum.dat > $1/vort_spectrum.dat


