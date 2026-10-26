# Nodal projection next to a no-slip wall: O(1) error in the viscous term

This test reproduces, with AMReX alone, an error that incflo shows in every viscous run with
no-slip walls: in the cells next to the wall, the wall-normal pressure gradient has an O(1)
error that does not converge. The incflo side, which shows the error with the unchanged
Burggraf benchmark, is in `test_no_eb_2d/wall_gp_error/` on the `wall_gp_error` branch of
https://github.com/cgilet/incflo.

## What it does

`main.cpp` takes the Burggraf velocity on the unit square, which is divergence-free and zero
on the walls except for the lid's tangential velocity. In 3D the velocity does not depend on
z, w = 0, and z is periodic. For n = 32, 64, 128 and 256 cells per side it:

1. projects the velocity once with `MLNodeLaplacian`, as AMReX-Hydro's `NodalProjector` does
   in incflo (walls are Neumann for phi, sigma = 1, zero velocity in the ghost cells), and
2. applies the viscous operators incflo uses, with Dirichlet walls and maxorder 2 or 3, to
   the exact velocity and to the projected one. The operators are `MLABecLaplacian` on each
   component (incflo's scalar solve) and `MLTensorOp` (incflo's default tensor solve).

It prints, in the wall cells and leaving out n/8 cells next to each corner:

- **The projected wall-normal velocity error** in the first four cells from the wall, divided
  by delta = (h^2/8) u_n''. Here u_n'' is the second normal derivative of the normal velocity
  at the wall.
- **The operator's wall error,** as the c in (L_h u - lap u)_n = -c u_n''.

## Build and run

```shell
make -j4 AMREX_HOME=/path/to/amrex       # DIM=3 for 3D
mpiexec -n 2 ./main2d.gnu.MPI.ex inputs
```

`n_cell`, `max_grid_size`, `rtol` and `verbose` can be set in `inputs` or on the command line.
It takes a few seconds. The output does not depend on `max_grid_size` (checked with 16, 32
and 1024), and 3D (8 cells in z) prints the same numbers as 2D.

## Output (AMReX 1aebac90c)

```
Projected Burggraf velocity: max error, and the wall-normal error in the
cells 0-3 from the wall divided by delta = (h^2/8) u_n''
     n   max |u_P - u|   cell 0   cell 1   cell 2   cell 3
    32       2.855e-03    0.903    0.735    0.599    0.489
    64       8.132e-04    0.950    0.858    0.774    0.699
   128       2.196e-04    0.975    0.926    0.880    0.836
   256       5.745e-05    0.987    0.962    0.938    0.914

Viscous operator in the wall cells: (L_h u - lap u)_n = -c u_n''
  operator                         n   c, exact u   c, projected u
  MLABecLaplacian, maxorder 2     32       0.2495           0.4987
                                  64       0.2499           0.4997
                                 128       0.2500           0.4999
                                 256       0.2500           0.5000
  MLABecLaplacian, maxorder 3     32       0.0190           0.3505
                                  64       0.0097           0.3426
                                 128       0.0049           0.3381
                                 256       0.0025           0.3358
  MLTensorOp,      maxorder 2     32       0.3269           0.6573
                                  64       0.3305           0.6626
                                 128       0.3320           0.6648
                                 256       0.3327           0.6658
  MLTensorOp,      maxorder 3     32       0.0340           0.4724
                                  64       0.0176           0.4596
                                 128       0.0090           0.4523
                                 256       0.0045           0.4485
```

## What it shows

1. **The projection.** After the projection the wall-normal velocity is too large by
   delta = (h^2/8) u_n'' in the first few cells: the ratio in the table tends to 1. This error
   is second order, but it is nearly the same in each cell, so a 1/h^2 stencil does not see
   it as smooth. In the wall cell, with the Dirichlet ghost cell, the Laplacian of the offset
   is -2 delta/h^2 = -u_n''/4 with maxorder 2, and -(8/3) delta/h^2 = -u_n''/3 with
   maxorder 3. `MLTensorOp` gives 4/3 of these in the wall-normal component.

   The cause is `mlndlap_divu` (`AMReX_MLNodeLap_2D_K.H`, `AMReX_MLNodeLap_3D_K.H`). A node
   on a wall has a control volume of half a cell, 0 < n < h/2. The stencil takes the
   tangential derivative of the tangential velocity at the first cell centres, n = h/2, for
   the whole half volume, instead of its average over it (at n = h/4 to second order). The
   two differ by (h/4) u_n''. Over the half-cell height h/2, the divergence constraint turns
   this into the offset delta = (h^2/8) u_n''.

2. **The viscous operator.** With maxorder 2 (incflo's default), the Dirichlet stencil is
   O(1) wrong in the wall cell even for the exact velocity: c = 1/4 for `MLABecLaplacian` and
   1/3 for `MLTensorOp`. With maxorder 3 this goes away (O(h)).

In the cells next to the wall, the two errors add up:

| operator | maxorder | c, exact u | c, projected u |
|---|---|---|---|
| `MLABecLaplacian` | 2 | 1/4 | 1/2 |
| `MLABecLaplacian` | 3 | 0 | 1/3 |
| `MLTensorOp` | 2 | 1/3 | 2/3 |
| `MLTensorOp` | 3 | 0 | 4/9 |

## How it shows up in incflo

The projection fixes the wall-normal velocity, so the time step cannot correct an O(1) error
in the wall-normal viscous force with the velocity. The wall-node pressure takes it instead,
which makes the pressure first order there and the pressure gradient O(1) wrong in the wall
cells:

gp_n - dp/dn = -c mu u_n''

with the c of the last column above. In incflo's Burggraf cavity at n = 128 on stock AMReX,
the c fitted on the lid is 0.500, 0.336, 0.666 and 0.450 for the four rows, in that order.
Velocity stays second order.

## A candidate fix, and what is still open

A local AMReX branch makes the nodal divergence second order at walls. At nodes on exactly
one wall, `MLNodeLaplacian::compRHS` adds a correction that takes the tangential velocity
(and the cell-centred source, if any) a quarter cell from the wall, extrapolated linearly from
the first two rows of cells. With it and maxorder 3, c goes to 0 at O(h) along the walls:

```
Projected Burggraf velocity: max error, and the wall-normal error in the
cells 0-3 from the wall divided by delta = (h^2/8) u_n''
     n   max |u_P - u|   cell 0   cell 1   cell 2   cell 3
    32       1.730e-03    0.071   -0.006   -0.068   -0.118
    64       4.130e-04    0.036   -0.007   -0.045   -0.080
   128       9.948e-05    0.018   -0.005   -0.026   -0.046
   256       2.429e-05    0.009   -0.003   -0.014   -0.025

Viscous operator in the wall cells: (L_h u - lap u)_n = -c u_n''
  operator                         n   c, exact u   c, projected u
  MLABecLaplacian, maxorder 2     32       0.2495           0.2770
                                  64       0.2499           0.2642
                                 128       0.2500           0.2573
                                 256       0.2500           0.2537
  MLABecLaplacian, maxorder 3     32       0.0190           0.0556
                                  64       0.0097           0.0288
                                 128       0.0049           0.0147
                                 256       0.0025           0.0074
  MLTensorOp,      maxorder 2     32       0.3269           0.3662
                                  64       0.3305           0.3509
                                 128       0.3320           0.3424
                                 256       0.3327           0.3380
  MLTensorOp,      maxorder 3     32       0.0340           0.0882
                                  64       0.0176           0.0458
                                 128       0.0090           0.0234
                                 256       0.0045           0.0118
```

**Corners remain a problem.** The branch keeps the original stencil at corner nodes. That
leaves a large error in the corner cells wherever u_n'' is not zero at the corner. At
incflo's Burggraf lid corners, the gp error goes from about 0.25 without the fix to only
0.165, 0.143 and 0.121 with it, at n = 32, 64 and 128. Applying the extrapolation at
corners as well makes those converge at first order (0.044, 0.023, 0.011). But in a
manufactured incflo test with div(u) = S != 0, the corner errors then stop converging. No
corner treatment found so far handles both cases, so the branch is not proposed yet. This
test leaves the corners out of its fits.
