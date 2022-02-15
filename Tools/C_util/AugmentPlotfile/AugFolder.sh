#!/bin/bash


if [[ -z "$1" ]]; then
    echo "Name of folder containing plotfiles must be supplied."
    exit
fi

echo -e "infile = tmp_infile\noutfile = tmp_outfile\nadd_vorticity = 1\nadd_divergence = 1" > tmp_inputs

for f in `ls -1 $1 | grep plt[0-9]*[0-9]`; do
    echo "Processing $1/${f}..."
    mv $1/$f tmp_infile
    mpirun -n 4 ./AugmentPlotfile3d.gnu.MPI.ex tmp_inputs
    mv tmp_outfile $1/${f}
    rm -R tmp_infile
done

rm tmp_inputs

