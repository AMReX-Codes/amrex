# Least Squares Tests
## Overview
This directory contains tests for:
1. Slope computation on face centroids using the Least Square method,
2. Laplacian operator computation when applying the explicit stencil,

Given a function defined on cell centroids.

# Building & Running
1. Build using `make`
   ```
   $ make DIM=<dim>
   ```
2. Run specifying an inputs file. Eg:
   ```
   $ ./main2d.gnu.DEBUG.MPI.ex inputs.2d.askew-y
   ```
3. View the resulting plot files through `amrvis` visualizaiton tool. Eg:
   ```
   $ amrvis2d plot-askew-y
   $ amrvis2d plot-askew-y-analytic
   ```
4. Compare numerical solution against analytical solution using `fcompare`. Eg:
   ```
   $ fcompare plot-askew-y plot-askew-y-analytic
   ```

## Tests
There are several input files which correspond to different 2D & 3D configurations of `Poisueille Flow through a Channel`. For eg: `inputs.3d.poiseuille.aligned.xy-x` corresponds to the flow through a 3D channel where the walls are along the xy plane and the flow is along the x direction, and the grid is aligned with the walls (does not cut the wall at an angle).
