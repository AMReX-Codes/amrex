# NodalProjectionWall

Tests `MLNodeLaplacian::setWallDivuCorrection`, the second-order nodal divergence at walls
(Neumann domain faces) and inflow faces.

## The problem

`MLNodeLaplacian::compRHS` computes, at each node, the divergence of the cell-centred velocity
over the node's control volume. At a node on a wall the control volume is cut in half, and the
stencil takes the tangential velocity of the first layer of cells, h/2 from the wall, as constant
over it instead of at its mean height, h/4. That makes the divergence first order there. After a
projection, the wall-normal velocity in the cells next to a no-slip wall is then off by

    delta = (h^2/8) u_n''

where u_n'' is the second normal derivative of the normal velocity at the wall. The offset is
nearly the same in each of the first few cells, a step from 0 at the wall, so a viscous operator
(which divides by h^2) sees an O(1) error in the wall cells, and the wall-normal pressure
gradient takes it: the pressure is only first order at the wall. The cell-centred source of
div(u) = S is integrated over the same half control volume, with the same error.

`setWallDivuCorrection(true)` takes these terms a quarter cell from the wall, extrapolated
linearly from the first two layers of cells, (5/4) layer 1 - (1/4) layer 2, in each wall direction
at domain edges and corners. Inflow faces are treated like walls.

## Build and run

    make -j4 DIM=2          # or DIM=3
    mpiexec -n 4 ./main2d.gnu.MPI.ex inputs-ci              # the checks, correction on
    mpiexec -n 4 ./main2d.gnu.MPI.ex inputs-ci wall_corr=0  # off
    mpiexec -n 4 ./main2d.gnu.MPI.ex inputs                 # the checks and the study

With CMake, the test is built with `-DAMReX_ENABLE_TESTS=ON` and run by CTest
(`ctest -R NodalProjectionWall`), with the correction on and off, in 2D and 3D.

Inputs: `n_cell` (resolutions, cells per side of the unit square, multiples of 8),
`max_grid_size`, `wall_corr` (0 off, 1 on), `checks` (1 to run the checks below), and for the
convergence study `study`, `burggraf`, `channel`, `box3d` (1 to run each part; `box3d` is 3D
only), `box3d_n_cell`, `rtol` and `verbose` (of the projection solves) and `maxorders` (of the
viscous operators).

## The checks

These call `compRHS` without a solve. The ones with a known answer stop the test with
`amrex::Abort` if they fail; it prints `NodalProjectionWall: all checks passed` otherwise. In 3D
the velocity does not depend on z, and z is periodic, except where noted.

- **Wall divergence.** u = (x y, -y^2/2) is divergence-free, with u_n'' = -1 on the bottom wall.
  The stencil is exact for it at interior nodes. On the bottom-wall nodes, rhs / (h/4) is 1
  without the correction and 0 with it. Checked to round-off.
- **Source.** u = (x (y+y^2), y (x+x^2)) with the cell-centred source -S = -div u: rhs is 0 if
  the velocity and the source are integrated over each control volume alike, which they are in
  both settings. Checked to round-off, at interior, wall and corner nodes.
- **Volume conservation.** With a random cell velocity, which changes by O(1) from cell to cell,
  the sum of rhs over the control volumes equals the flux through the domain boundary if the
  fluxes through the faces of the control volumes telescope. Closed box, and walls with inflow
  faces at both x ends; in 3D also with walls at both z ends, so that the box has nodes on three
  walls. Checked to 1e-12, in both settings: the correction redefines each half-cell face flux
  next to a wall once, and both nodes that share the face use it.
- **Two levels.** A fine level (ratio 2) over the coarse cells [3n/8, 5n/8) x [0, n/4), on the
  bottom wall. The correction is added level by level, so the wall nodes of both levels away
  from the coarse/fine boundary behave as in the first check (checked). At the coarse/fine nodes,
  `compRHS` rebuilds the rhs from the coarse cells and the fine level without the correction. The
  test prints them: where the coarse/fine boundary meets the wall, they are corners of the
  coarse/fine boundary, where the stencil is not exact even for u = (y, 0) (rhs -1/8, 1/8 and 1/16
  at every resolution, with and without the correction), so the missing O(h) wall correction
  does not change their order.
- **Burggraf corners** (printed only): rhs / h^2 at the 3x3 nodes at the (0,0) and lid (0,n)
  corners of the Burggraf cavity flow, to show where the corner error sits.

## The source and volume conservation

Without a source, the correction keeps the nodal volume conservation exact. The source part does
not: at a wall node it replaces the source's integral over the half control volume,
(h/2) S_1 in 2D, by (h/2) ((5/4) S_1 - (1/4) S_2), which is not a difference along the wall and
does not telescope. The sum over the control volumes then differs from the sum of the cell
source by O(h^2) for a smooth source. Taking the same amount back from the next node in would
restore volume conservation but brings back the O(1) wall error, since it recreates the step in
the normal velocity one cell further in. In a closed domain the solver removes the difference as
a uniform offset of the rhs; with an outflow, the outflow carries it. Momentum is not affected,
and neither are quantities transported with separately projected face velocities.

## The convergence study

With `study = 1` the test projects known divergence-free velocities with `MLNodeLaplacian`, the
way incflo calls AMReX-Hydro's `NodalProjector` (velocity ghost cells set to 0, sigma = 1), and
applies the viscous operators incflo uses to the projected velocity: `MLABecLaplacian` per
component and `MLTensorOp`, with Dirichlet walls and `maxorder` 2, 3 and 4. It prints the errors
and checks nothing. The numbers below are from `inputs` with 8 MPI ranks, `wall_corr = 0` (off)
and `1` (on); h = 1/n.

### Burggraf cavity (`burggraf`)

The Burggraf flow on [0,1]^2: zero on the walls except the tangential velocity of the lid at
y = 1. In 3D it does not depend on z, and z is periodic.

The wall-normal velocity error in the wall cells, divided by delta = (h^2/8) u_n'' (a fit along
the walls away from the corners), and the largest error in the lid-corner cells, divided by h^2:

| n | wall cell / delta, off | on | lid corner / h^2, off | on |
|---|---|---|---|---|
| 32 | 0.903 | 0.067 | 2.92 | 0.498 |
| 64 | 0.950 | 0.035 | 3.33 | 0.285 |
| 128 | 0.975 | 0.018 | 3.60 | 0.154 |
| 256 | 0.987 | 0.009 | 3.77 | 0.080 |

The viscous operator's error in the wall cells, (L_h u - lap u)_n = -c u_n'' for the normal
component, applied to the projected velocity (applied to the exact velocity, c is 0.0025 and
0.0045 at n = 256 for the two operators with `maxorder = 3`):

| n | `MLABecLaplacian`, maxorder 3, off | on | `MLTensorOp`, maxorder 3, off | on |
|---|---|---|---|---|
| 32 | 0.351 | 0.056 | 0.472 | 0.088 |
| 64 | 0.343 | 0.029 | 0.460 | 0.046 |
| 128 | 0.338 | 0.015 | 0.452 | 0.023 |
| 256 | 0.336 | 0.0074 | 0.449 | 0.012 |

With `maxorder = 2` the operator's own wall stencil keeps c at 0.25 (0.33 for `MLTensorOp`) for
the exact velocity, with or without the correction, so a viscous solve needs `maxorder = 3` as
well to benefit fully.

The correction does not remove the tangential components' O(1) viscous error in the wall cells:
the projection leaves an O(h^2) tangential velocity at the wall, which the Dirichlet stencil sees
over half a cell. Its maximum is 5.3 without the correction and 4.2 with it at n = 256
(`MLABecLaplacian`, maxorder 3), against 0.03 for the exact velocity. Unlike the normal error,
which the projection forces into the pressure, a time step can balance it with an O(h^2) change
of the tangential velocity.

### Inflow and outflow (`channel`)

A channel with no-slip walls at y = 0 and 1, an inflow face at x = 0 (`LinOpBCType::inflow`, the
inflow velocity in the ghost cells) and a pressure outflow at x = 1 (phi = 0), with u_n'' != 0 on
the walls and on the inflow face. The nodal divergence of the exact velocity, max |rhs| / h on
the walls, and the normal velocity error / delta of the projected velocity in the wall cells
(away from the ends), in the wall cell next to the outflow, and in the cells next to the inflow
face:

| n | wall rhs / h, off | on | walls, off | on | next to outflow, off | on | inflow, off | on |
|---|---|---|---|---|---|---|---|---|
| 32 | 0.880 | 0.112 | 0.929 | 0.132 | 0.378 | 0.065 | 0.892 | -0.048 |
| 64 | 0.908 | 0.057 | 0.964 | 0.069 | 0.383 | 0.034 | 0.944 | -0.025 |
| 128 | 0.923 | 0.029 | 0.982 | 0.035 | 0.385 | 0.018 | 0.972 | -0.013 |
| 256 | 0.930 | 0.015 | 0.991 | 0.018 | 0.386 | 0.009 | 0.986 | -0.007 |

An inflow face needs the same correction as a wall, and gets it. The outflow nodes are Dirichlet
and never corrected; the wall nodes next to them are, and converge like the rest of the wall.

Then a uniform inflow, u = 1 in the inflow ghost cells up to the no-slip walls, where the data do
not match at the wall/inflow corners. The flux of the projected velocity through x = 1/4, 1/2 and
3/4 minus the inflow flux: measured with the face fluxes each stencil uses, it is the solver
tolerance, about 1e-12, with and without the correction, so both conserve volume. Measured as
h sum_j u, the corrected stencil differs from it by its O(h^2) wall terms: 2.0e-4, 5.3e-5,
1.3e-5, 3.4e-6 at x = 1/2 for n = 32 to 256.

### Closed box in 3D (`box3d`)

u = grad(Phi) x (1,2,3) with Phi = F(x) F(y) H(z), F = 16 s^2 (1-s)^2, H = 8 s^2 (1-s), in a box
with walls on all six faces, one of them a lid moving tangentially: the case with cells on two
walls (edges) and three (corners). The nodal divergence of the exact velocity, max |rhs| / h, and
the wall-normal error / delta of the projected velocity in the face cells:

| n | faces, off | on | edges, off | on | face cell / delta, off | on |
|---|---|---|---|---|---|---|
| 16 | 91.5 | 23.7 | 113 | 27.6 | 0.763 | 0.204 |
| 32 | 107 | 13.6 | 120 | 15.4 | 0.875 | 0.117 |
| 64 | 117 | 7.35 | 124 | 8.09 | 0.935 | 0.063 |

The divergence is second order on the faces and along the edges with the correction, and first
order without it. The viscous c in the face cells (`MLABecLaplacian`, maxorder 3) is 0.347 off and
0.042 on at n = 64. The tangential components keep their O(1) viscous error, as in 2D.
