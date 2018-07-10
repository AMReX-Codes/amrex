#!/bin/bash


if [[ -z "$1" ]]; then
    echo "Name of folder containing plotfiles must be supplied."
    exit
fi

for f in `ls -1 $1 | grep plt[0-9]*[0-9]`; do
    echo "Processing $1/${f}..."
    mv $1/$f tmp_plotfile
    mpirun -n 4 ./AmrDeriveSpectrum3d.gnu.MPI.ex tmp_inputs
    mv tmp_plotfile $1/${f}
done

