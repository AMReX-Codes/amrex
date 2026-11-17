# NodalProjectionWall

Tests `MLNodeLaplacian::setBoundaryDivuCorrection`, the second-order nodal divergence at nodes on
`Neumann` and `inflow` domain faces (walls and inflow faces, both called walls below).

## The problem

`MLNodeLaplacian::compRHS` computes, at each node, the divergence of the cell-centred velocity
over the node's control volume. At a node on a wall the control volume is cut in half, and the
stencil takes the tangential velocity of the first layer of cells, h/2 from the wall, as constant
over it instead of at its mean height, h/4. That makes the divergence first order there, and
after a projection the wall-normal velocity in the cells next to the wall is off by
(h^2/8) u_n'', where u_n'' is the second normal derivative of the normal velocity at the wall.
The cell-centred source of the projection is integrated over the same half control volume, with
the same error.

On its own this is harmless, but the default cell-centred viscous operators (`MLABecLaplacian`,
`MLTensorOp`) with a Dirichlet condition on the wall divide it by h^2: the viscous term then has
an O(1) error in the wall-normal component in the cells next to the wall, which a projection
method passes on to the pressure gradient. With `maxorder = 2` the viscous stencil has an O(1)
wall error of its own, so a viscous solve needs `maxorder = 3` as well.

`setBoundaryDivuCorrection(true)` takes the velocity and the source a quarter cell from the wall,
extrapolated linearly from the first two layers of cells, (5/4) layer 1 - (1/4) layer 2, in each
wall direction at domain edges and corners. Inflow faces are treated like walls: their control
volumes are cut in half in the same way, and the ghost cells hold the inflow velocity on the
face. Outflow faces (`Dirichlet`, phi = 0) have no divergence equation and are not affected.

With EB, the correction is applied at the wall nodes whose nearby cells are all regular (the first
two layers of cells and the ghost cells, on both sides of the node), so the result does not depend
on the box layout. Where it stops along a wall, next to the EB, the cells at the transition have an
O(1) error in the tangential pressure gradient, about as large as the error that the uncorrected
wall has in the normal one.

Without a source, the correction keeps the nodal volume conservation exact. With one, the sum of
the corrected right-hand side over the control volumes differs from the sum of the cell source
by O(h^2): in a closed domain the solver removes the difference as a uniform offset, and with an
outflow, the outflow carries it.

## Build and run

    make -j4 DIM=2          # or DIM=3
    mpiexec -n 4 ./main2d.gnu.MPI.ex inputs-ci              # correction on
    mpiexec -n 4 ./main2d.gnu.MPI.ex inputs-ci wall_corr=0  # off

With CMake, the test is built with `-DAMReX_ENABLE_TESTS=ON` and run by CTest
(`ctest -R NodalProjectionWall`), with the correction on and off, in 2D and 3D. Inputs: `n_cell`
(resolutions, cells per side of the unit square, multiples of 8), `max_grid_size`, and
`wall_corr` (0 off, 1 on).

## The checks

The checks call `compRHS` without a solve, stop the test with `amrex::Abort` if one fails, and
print `NodalProjectionWall: all checks passed` otherwise. In 3D the velocity does not depend on z,
and z is periodic, except where noted.

- **Wall divergence.** u = (x y, -y^2/2) is divergence-free, with u_n'' = -1 on the bottom wall.
  The stencil is exact for it at interior nodes. On the bottom-wall nodes, rhs / (h/4) is 1
  without the correction and 0 with it. Checked to round-off.
- **Source.** u = (x (y+y^2), y (x+x^2)) with the cell-centred source -S = -div u: rhs is 0 if
  the velocity and the source are integrated over each control volume alike, which they are in
  both settings. Checked to round-off, at interior, wall and corner nodes.
- **Volume conservation.** With a random cell velocity, the sum of rhs over the control volumes
  equals the flux through the domain boundary if the fluxes through the faces of the control
  volumes telescope. Closed box, and walls with inflow faces at both x ends; in 3D also with walls
  at both z ends, so that the box has nodes on three walls. Checked to 1e-12 in both settings.
- **Two levels.** A fine level (ratio 2) on the bottom wall: the wall nodes of both levels away
  from the coarse/fine boundary behave as in the first check, and give 0 for the linear shear
  u = (y, 0). The coarse/fine nodes are rebuilt by `compRHS` without the correction.
