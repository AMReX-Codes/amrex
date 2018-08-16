.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


MLMG and Linear Operator Classes
================================

``MLMG`` is a class for solving the linear system using the geometric
multigrid method.  The constructor of ``MLMG`` takes the reference to
:cpp:`MLLinOp`, an abstract base class of various linear operator
classes, :cpp:`MLABecLaplacian`, :cpp:`MLPoisson`,
:cpp:`MLNodeLaplacian`, etc.  We choose the type of linear operator
class according to the type the linear system to solve.

- :cpp:`MLABecLaplacian` for cell-centered canonical form (equation :eq:`eqn::abeclap`).

- :cpp:`MLPoisson` for cell-centered constant coefficient Poisson's
  equation :math:`\nabla^2 \phi = f`.

- :cpp:`MLNodeLaplacian` for nodal variable coefficient Poisson's
  equation :math:`\nabla \cdot (\sigma \nabla \phi) = f`.

The constructors of these linear operator classes are in the form like
below

.. highlight:: c++

::

    MLABecLaplacian (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo(),
                     const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
     
It takes :cpp:`Vectors` of :cpp:`Geometry`, :cpp:`BoxArray` and
:cpp:`DistributionMapping`.  They are :cpp:`Vectors` because MLMG can
do multi-level composite solve.  If you are using it for single-level,
you can do

.. highlight:: c++

::

    // Given Geometry geom, BoxArray grids, and DistributionMapping dmap on single level
    MLABecLaplacian mlabeclap({geom}, {grids}, {dmap});

to let the compiler construct :cpp:`Vectors` for you.  We have seen
classes :cpp:`Vector`, :cpp:`Geometry`, :cpp:`BoxArray`, and
:cpp:`DistributionMapping` in chapter :ref:`Chap:Basics`.  There are
two new classes that are optional parameters.  :cpp:`LPInfo` is a
class for passing parameters.  :cpp:`FabFactory` is used in problems
with embedded boundaries (chapter :ref:`Chap:EB`).

After the linear operator is built, we need to set up boundary
conditions.  This will be discussed later in section
:ref:`sec:linearsolver:bc`.

For :cpp:`MLABecLaplacian`, we then need to call member functions

.. highlight:: c++

::

    void setScalars (Real alpha, Real beta);
    void setACoeffs (int amrlev, const MultiFab& A);
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& B);

to set up the coefficients.  For :cpp:`MLPoisson`, there are no
coefficients to set.  For :cpp:`MLNodeLaplacian`, one needs to call
member function

.. highlight:: c++

::

    void setSigma (int amrlev, const MultiFab& a_sigma);

The :cpp:`int amrlev` parameter should be zero for single-level
solves.  For multi-level solves, each level needs to be provided with
``A`` and ``B``, or ``Sigma``.  For composite solves, level 0 here
means the lowest level for the solver, not necessarily the lowest
level in the entire AMR hierarchy, because it might be doing a solve
on AMR levels 3 to 5.

After boundary conditions and coefficients are set up, the linear
operator is ready for an MLMG object like below.

.. highlight:: C++

::

    MLMG mlmg(mlabeclaplacian);

Various parameters (section :ref:`sec:linearsolver:pars`) can be
optionally set, and then we can ``MLMG`` member function

.. highlight:: C++

::

    Real solve (const Vector<MultiFab*>& a_sol,
                const Vector<MultiFab const*>& a_rhs,
                Real a_tol_rel, Real a_tol_abs);

to solve the problem given an initial guess and a right-hand side.
Zero is a perfectly fine initial guess.  The two :cpp:`Reals` are the
targeted relative and absolute tolerances.  The solver will stop when
one of the targets is met.  Set the absolute tolerance to zero if one
does not have a good value for it.  The return value of :cpp:`solve`
is the max-norm error.

After the solver returns successfully, if needed, we can call 

.. highlight:: c++

::

    void compResidual (const Vector<MultiFab*>& a_res,
                       const Vector<MultiFab*>& a_sol,
                       const Vector<MultiFab const*>& a_rhs);

to compute residual (i.e., :math:`f - L(\phi)`) given the solution and
the righ-hand side.  For cell-centered solvers, we can also call the
following functions to compute gradient :math:`\nabla \phi` and fluxes
:math:`-B \nabla \phi`.

.. highlight:: c++

::

    void getGradSolution (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_grad_sol);
    void getFluxes       (const Vector<Array<MultiFab*,AMREX_SPACEDIM> >& a_fluxes);


.. _sec:linearsolver:bc:

Boundary Conditions
===================

We now discuss how to set up boundary conditions for linear operators.
In the following, physical domain boundaries refer to the boundaries
of the physical domain, whereas coarse/fine boundaries refer to the
level boundaries of fine AMR levels.  The following steps must be
followed in the exact order.

First we need to set physical domain boundary types via :cpp:`MLLinOp` member
function

.. highlight:: c++

::

    void setDomainBC (const Array<BCType,AMREX_SPACEDIM>& lobc,  // for lower ends
                      const Array<BCType,AMREX_SPACEDIM>& hibc); // for higher ends

The supported BC types at
the physical domain boundaries are

- :cpp:`LinOpBCType::Periodic` for periodic boundary.

- :cpp:`LinOpBCType::Dirichlet` for Dirichlet boundary condition.

- :cpp:`LinOpBCType::Neumann` for homogeneous Neumann boundary condition.

- :cpp:`LinOpBCType::reflect_odd` for reflection with sign changed.

The second step is to set up coarse/fine Dirichlet boundary
conditions.  This step is not always needed.  But when it's needed, it
must be called before step 3.  This step is not needed when the
coarsest level in the solver covers the whole computational domain
(e.g., the coarsest level is AMR level 0).  Note that this step is
still needed in the case that the solver is used to do a single-level
solve on a fine AMR level not covering the whole domain.  The
:cpp:`MLLinOp` member function for this step is

.. highlight:: c++

::

    void setCoarseFineBC (const MultiFab* crse, int crse_ratio);

Here :cpp:`const MultiFab* crse` contains the Dirichlet boundary
values at the coarse resolution, and :cpp:`int crse_ratio` (e.g., 2)
is the refinement ratio between the coarsest solver level and the AMR
level below it.

In the third step, for each level, we call :cpp:`MLLinOp` member
function

.. highlight:: c++

::

    virtual void setLevelBC (int amrlev, const MultiFab* levelbcdata) = 0;

to set up Dirichlet boundary values.  Here the ``MultiFab`` must have
one ghost cell.  Although ``MultiFab`` is used to pass the data, only
the values in the ghost cells at Dirichlet boundaries are used.  If
there are no Dirichlet boundaries, we still need to make this call,
but we could call it with :cpp:`nullptr`.  It should be emphasized
that the data in ``levelbcdata`` for Dirichlet boundaries are assumed
to be exactly on the face of the physical domain even though for
cell-centered solvers the sought unknowns are at cell centers.  And
for cell-centered solvers, the ``MultiFab`` argument must have the
cell-centered type.

.. _sec:linearsolver:pars:

Parameters
==========

There are many parameters that can be set.  Here we discuss some
commonly used ones.

:cpp:`MLLinOp::setVerbose(int)`, :cpp:`MLMG::setVerbose(int)` and
:cpp:`MLMG:setBottomVerbose(int)` can be control the verbosity of the
linear operator, multigrid solver and the bottom solver, respectively.

The multigrid solver is an iterative solver.  The maximal number of
iterations can be changed with :cpp:`MLMG::setMaxIter(int)`.  We can
also do a fixed number of iterations with
:cpp:`MLMG::setFixedIter(int)`.  By default, V-cycle is used.  We can
use :cpp:`MLMG::setMaxFmgIter(int)` to control how many full multigrid
cycles can be done before switching to V-cycle.

:cpp:`LPInfo::setMaxCoarseningLevel(int)` can be used to control the
maximal number of multigrid levels.  We usually should not call this
function.  However, we sometimes build the solver to simply apply the
operator (e.g., :math:`L(\phi)`) without needing to solve the system.
We can do something as follows to avoid the cost of building coarsened
operators for the multigrid.

.. highlight:: c++

::

    MLABecLaplacian mlabeclap({geom}, {grids}, {dmap}, LPInfo().setMaxCoarseningLevel(0));
    // set up BC
    // set up coefficients
    MLMG mlmg(mlabeclap);
    // out = L(in)
    mlmg.apply(out, in);  // here both in and out are const Vector<MultiFab*>&

At the bottom of the multigrid cycles, we use the biconjugate gradient
stabilized method as the bottom solver.  :cpp:`MLMG` member method

.. highlight:: c++

::

    void setBottomSolver (BottomSolver s);

can be used to change the bottom solver.  Available choices are

- :cpp:`MLMG::BottomSolver::smoother`: Smoother such as Gauss-Seidel.

- :cpp:`MLMG::BottomSolver::bicgstab`: The default.

- :cpp:`MLMG::BottomSolver::cg`: The conjugate gradient method.  The
  matrix must be symmetric.

- :cpp:`MLMG::BottomSolver::Hypre`: BoomerAMG in HYPRE.  Currently for
  cell-centered only.

Curvilinear Coordinates
=======================

The linear solvers support curvilinear coordinates including 1D
spherical and 2d cylindrical :math:`(r,z)`.  In those cases, the
divergence operator has extra metric terms.  If one does not want the
solver to include the metric terms because they have been handled in
other ways, one can call :cpp:`setMetricTerm(bool)` with :cpp:`false`
on the :cpp:`LPInfo` object passed to the constructor of linear
operators. 

HYPRE
=====

AMReX can use HYPRE BoomerAMG as a bottom solver (currently for
cell-centered problems only), as we have mentioned.  For challenging
problems, our geometric multigrid solver may have difficulty solving,
whereas an algebraic multigrid method might be more robust.  We note
that by default our solver always tries to geometrically coarsen the
problem as much as possible.  However, as we have mentioned, we can
call :cpp:`setMaxCoarseningLevel(0)` on the :cpp:`LPInfo` object
passed to the constructor of a linear operator to disable the
coarsening completely.  In that case the bottom solver is solving the
residual correction form of the original problem.  To use HYPRE, one
must include ``amrex/Src/Extern/HYPRE`` in the build system.

Embedded Boundaries
===================

AMReX support solving linear systems with embedded boundaries.  See
chapter :ref:`Chap:EB` for more details.

.. solver reuse

