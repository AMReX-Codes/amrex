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
:cpp:`DistributionMapping`.  The arguments are :cpp:`Vectors` because MLMG can
do multi-level composite solve.  If you are using it for single-level,
you can do

.. highlight:: c++

::

    // Given Geometry geom, BoxArray grids, and DistributionMapping dmap on single level
    MLABecLaplacian mlabeclaplacian({geom}, {grids}, {dmap});

to let the compiler construct :cpp:`Vectors` for you.  Recall that the
classes :cpp:`Vector`, :cpp:`Geometry`, :cpp:`BoxArray`, and
:cpp:`DistributionMapping` are defined in chapter :ref:`Chap:Basics`.  There are
two new classes that are optional parameters.  :cpp:`LPInfo` is a
class for passing parameters.  :cpp:`FabFactory` is used in problems
with embedded boundaries (chapter :ref:`Chap:EB`).

After the linear operator is built, we need to set up boundary
conditions.  This will be discussed later in section
:ref:`sec:linearsolver:bc`.

For :cpp:`MLABecLaplacian`, we next need to call member functions

.. highlight:: c++

::

    void setScalars (Real A, Real B);
    void setACoeffs (int amrlev, const MultiFab& alpha);
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta);

to set up the coefficients for equation :eq:`eqn::abeclap`. This is unnecessary for
:cpp:`MLPoisson`, as there are no coefficients to set.  For :cpp:`MLNodeLaplacian`,
one needs to call the member function

.. highlight:: c++

::

    void setSigma (int amrlev, const MultiFab& a_sigma);

The :cpp:`int amrlev` parameter should be zero for single-level
solves.  For multi-level solves, each level needs to be provided with
``alpha`` and ``beta``, or ``Sigma``.  For composite solves, :cpp:`amrlev` 0 will
mean the lowest level for the solver, which is not necessarily the lowest
level in the AMR hierarchy. This is so solves can be done on different sections
of the AMR hierarchy, e.g. on AMR levels 3 to 5.

After boundary conditions and coefficients are prescribed, the linear
operator is ready for an MLMG object like below.

.. highlight:: C++

::

    MLMG mlmg(mlabeclaplacian);

Optional parameters can be set (see section :ref:`sec:linearsolver:pars`),
and then we can use the ``MLMG`` member function

.. highlight:: C++

::

    Real solve (const Vector<MultiFab*>& a_sol,
                const Vector<MultiFab const*>& a_rhs,
                Real a_tol_rel, Real a_tol_abs);

to solve the problem given an initial guess and a right-hand side.
Zero is a perfectly fine initial guess.  The two :cpp:`Reals` in the argument
list are the targeted relative and absolute error tolerances.
The solver will terminate when one of these targets is met.
Set the absolute tolerance to zero if one
does not have a good value for it.  The return value of :cpp:`solve`
is the max-norm error.

After the solver returns successfully, if needed, we can call

.. highlight:: c++

::

    void compResidual (const Vector<MultiFab*>& a_res,
                       const Vector<MultiFab*>& a_sol,
                       const Vector<MultiFab const*>& a_rhs);

to compute residual (i.e., :math:`f - L(\phi)`) given the solution and
the right-hand side.  For cell-centered solvers, we can also call the
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
boundaries between AMR levels. The following steps must be
followed in the exact order.

1) For any type of solver, we first need to set physical domain boundary types via the :cpp:`MLLinOp` member
function

.. highlight:: c++

::

    void setDomainBC (const Array<BCType,AMREX_SPACEDIM>& lobc,  // for lower ends
                      const Array<BCType,AMREX_SPACEDIM>& hibc); // for higher ends

The supported BC types at the physical domain boundaries are

- :cpp:`LinOpBCType::Periodic` for periodic boundary.

- :cpp:`LinOpBCType::Dirichlet` for Dirichlet boundary condition.

- :cpp:`LinOpBCType::Neumann` for homogeneous Neumann boundary condition.

- :cpp:`LinOpBCType::inhomogNeumann` for inhomogeneous Neumann boundary condition.

- :cpp:`LinOpBCType::reflect_odd` for reflection with sign changed.

2) Cell-centered solvers only: 
if we want to do a linear solve where the boundary conditions on the 
coarsest AMR level of the solve come from a coarser level (e.g. the
base AMR level of the solve is > 0 and does not cover the entire domain), 
we must explicitly provide the coarser data.  Boundary conditions from a 
coarser level are always Dirichlet.  

Note that this step, if needed, must be performed before the step below.  
The :cpp:`MLLinOp` member function for this step is

.. highlight:: c++

::

    void setCoarseFineBC (const MultiFab* crse, int crse_ratio);

Here :cpp:`const MultiFab* crse` contains the Dirichlet boundary
values at the coarse resolution, and :cpp:`int crse_ratio` (e.g., 2)
is the refinement ratio between the coarsest solver level and the AMR
level below it.  The MultiFb crse does not need to have ghost cells itself. 
If the coarse grid bc's for the solve are identically zero, :cpp:`nullptr` 
can be passed instead of :cpp:`crse`.

3) Cell-centered solvers only: 
before the solve one must always call the :cpp:`MLLinOp` member function 

.. highlight:: c++

::

    virtual void setLevelBC (int amrlev, const MultiFab* levelbcdata) = 0;

If we want to supply any inhomogeneous Dirichlet or Neumann boundary 
conditions at the domain boundaries, we must supply those values 
in ``MultiFab* levelbcdata``, which must have at least one ghost cell. 
Note that the argument :cpp:`amrlev` is relative to the solve, not
necessarily the full AMR hierarchy; amrlev = 0 refers to the coarsest
level of the solve.

If the boundary condition is Dirichlet the ghost cells outside the
domain boundary of ``levelbcdata`` must hold the value of the solution
at the domain boundary; 
if the boundary condition is Neumann those ghost cells must hold
the value of the gradient of the solution normal to the boundary
(e.g. it would hold dphi/dx on both the low and high facees in the x-direction).

If the boundary conditions contain no inhomogeneous Dirichlet or Neumann boundaries,
we can pass :cpp:`nullptr` instead of a MultiFab.

We can use the solution array itself to hold these values;
the values are copied to internal arrays and will not be over-written
when the solution array itself is being updated by the solver. 
Note, however, that this call does not provide an initial guess for the solve.

It should be emphasized that the data in ``levelbcdata`` for 
Dirichlet or Neumann boundaries are assumed to be exactly on the face 
of the physical domain; storing these values in the ghost cell of
a cell-centered array is a convenience of implementation.

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

- :cpp:`MLMG::BottomSolver::bicgstab`: The default.

- :cpp:`MLMG::BottomSolver::cg`: The conjugate gradient method.  The
  matrix must be symmetric.

- :cpp:`MLMG::BottomSolver::smoother`: Smoother such as Gauss-Seidel.

- :cpp:`MLMG::BottomSolver::bicgcg`: Start with bicgstab. Switch to cg
  if bicgstab fails.  The matrix must be symmetric.

- :cpp:`MLMG::BottomSolver::cgbicg`: Start with cg. Switch to bicgstab
  if cg fails.  The matrix must be symmetric.

- :cpp:`MLMG::BottomSolver::hypre`: BoomerAMG in hypre.

- :cpp:`MLMG::BottomSolver::petsc`: Currently for cell-centered only.

Curvilinear Coordinates
=======================

The linear solvers support curvilinear coordinates including 1D
spherical and 2d cylindrical :math:`(r,z)`.  In those cases, the
divergence operator has extra metric terms.  If one does not want the
solver to include the metric terms because they have been handled in
other ways, one can call :cpp:`setMetricTerm(bool)` with :cpp:`false`
on the :cpp:`LPInfo` object passed to the constructor of linear
operators.

External Solvers
================

AMReX can use the `hypre <https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`_ algebraic multigrid solver, BoomerAMG, 
as a bottom solver for both cell-centered and node-based problems.
For challenging problems, our geometric multigrid solver may have difficulty solving,
whereas an algebraic multigrid method might be more robust.  
We note that by default our solver always tries to geometrically coarsen the
problem as much as possible.  However, as we have mentioned, we can
call :cpp:`setMaxCoarseningLevel(0)` on the :cpp:`LPInfo` object
passed to the constructor of a linear operator to disable the
coarsening completely.  In that case the bottom solver is solving the
residual correction form of the original problem.  

To use hypre, one must include ``amrex/Src/Extern/HYPRE`` in the build system. 
For an example of using hypre, we refer the reader to
``Tutorials/LinearSolvers/ABecLaplacian_C``.

AMReX can also use `PETSc <https://www.mcs.anl.gov/petsc/>`_ as a bottom solver for cell-centered
problems.  To use PETSc, one must include ``amrex/Src/Extern/PETSc``
in the build system.  For an example of using PETSc, we refer the
reader to ``Tutorials/LinearSolvers/ABecLaplacian_C``.

MAC Projection
=========================

Some codes define a velocity field :math:`U = (u,v,w)` on faces, i.e. 
:math:`u` is defined on x-faces, :math:`v` is defined on y-faces,
and :math:`w` is defined on z-faces.   We refer to the exact projection 
of this velocity field as a MAC projection, in which we solve 

.. math::

   D( \beta \nabla \phi) = D(U^*) - S

for :math:`\phi` and then set 

.. math::

   U = U^* - \beta \nabla \phi


where :math:`U^*` is a vector field (typically velocity) that we want to satisfy 
:math:`D(U) = S`.  For incompressible flow,  :math:`S = 0`.

The MacProjection class can be defined and used to perform the MAC projection without explicitly
calling the solver directly.  In addition to solving the variable coefficient Poisson equation,
the MacProjector internally computes the divergence of the vector field, :math:`D(U^*)`,
to compute the right-hand-side, and after the solve, subtracts the weighted gradient term to
make the vector field result satisfy the divergence constraint.  

In the simplest form of the call, :math:`S` is assumed to be zero and does not need to be specified.
Typically, the user does not allocate the solution array, but it is also possible to create and pass
in the solution array and have :math:`\phi` returned as well as :math:`U`.  

Caveat:  Currently the MAC projection only works when the base level covers the full domain; it does
not yet have the interface to pass boundary conditions for a fine level that come from coarser data.

Also note that any Dirichlet or Neumann boundary conditions at domain boundaries
are assumed to be homogeneous.  The call to the :cpp:`MLLinOp` member function 
:cpp:`setLevelBC` occurs inside the MacProjection class; one does not need to call that
explicitly when using the MacProjection class.

The code below is taken from 
``Tutorials/LinearSolvers/MAC_Projection_EB/main.cpp`` and demonstrates how to set up 
the MACProjector object and use it to perform a MAC projection.

.. highlight:: c++

::

    EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

    // allocate face-centered velocities and face-centered beta coefficient
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        vel[idim].define (amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 1,
                          MFInfo(), factory);
        beta[idim].define(amrex::convert(grids,IntVect::TheDimensionVector(idim)), dmap, 1, 0,
	                  MFInfo(), factory);
        beta[idim].setVal(1.0);  // set beta to 1
    }

    // If we want to use phi elsewhere, we must create an array in which to return the solution 
    // MultiFab phi_inout(grids, dmap, 1, 1, MFInfo(), factory);

    // If we want to supply a non-zero S we must allocate and fill it outside the solver
    // MultiFab S(grids, dmap, 1, 0, MFInfo(), factory);
    // Set S here ... 

    // set initial velocity to U=(1,0,0)
    AMREX_D_TERM(vel[0].setVal(1.0);,
                 vel[1].setVal(0.0);,
                 vel[2].setVal(0.0););

    LPInfo lp_info;

    // If we want to use hypre to solve the full problem we do not need to coarsen the GMG stencils
    if (use_hypre_as_full_solver)
        lp_info.setMaxCoarseningLevel(0);

    MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // face-based velocity
                         {amrex::GetArrOfConstPtrs(beta)}, // beta
                         {geom},                           // the geometry object
                         lp_info);                         // structure for passing info to the operator

    // Here we specifiy the desired divergence S
    // MacProjector macproj({amrex::GetArrOfPtrs(vel)},       // face-based velocity
    //                      {amrex::GetArrOfConstPtrs(beta)}, // beta
    //                      {geom},                           // the geometry object
    //                      lp_info,                          // structure for passing info to the operator
    //                      {&S});                            // defines the specified RHS divergence

    // Set bottom-solver to use hypre instead of native BiCGStab 
    if (use_hypre_as_full_solver || use_hypre_as_bottom_solver) 
       macproj.setBottomSolver(MLMG::BottomSolver::hypre);

    // Set boundary conditions.
    //  Here we use Neumann on the low x-face, Dirichlet on the high x-face,
    //  and periodic in the other two directions  
    //  (the first argument is for the low end, the second is for the high end)
    // Note that Dirichlet and Neumann boundary conditions are assumed to be homogeneous.
    macproj.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                      LinOpBCType::Periodic,
                                      LinOpBCType::Periodic)},
                        {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                      LinOpBCType::Periodic,
                                      LinOpBCType::Periodic)});

    macproj.setVerbose(mg_verbose);
    macproj.setCGVerbose(cg_verbose);

    // Define the relative tolerance
    Real reltol = 1.e-8;

    // Define the absolute tolerance; note that this argument is optional
    Real abstol = 1.e-15;

    // Solve for phi and subtract from the velocity to make it divergence-free
    macproj.project(reltol,abstol);

    // If we want to use phi elsewhere, we can pass in an array in which to return the solution 
    // macproj.project({&phi_inout},reltol,abstol);


See ``Tutorials/LinearSolvers/MAC_Projection_EB`` for the complete working example.

Nodal Projection
================

Some codes define a velocity field :math:`U = (u,v,w)` with all
components co-located on cell centers.  The nodal solver in AMReX 
can be used to compute an approximate projection of the cell-centered
velocity field, with pressure and velocity divergence defined on nodes.
When we use the nodal solver this way, and subtract only the cell average
of the gradient from the velocity, it is effectively an approximate projection.

As with the MAC projection, consider that we want to solve 

.. math::

   D( \beta \nabla \phi) = D(U^*) - S

for :math:`\phi` and then set 

.. math::

   U = U^* - \beta \nabla \phi

where :math:`U^*` is a vector field defined on cell centers and we want to satisfy
:math:`D(U) = S`.  For incompressible flow,  :math:`S = 0`.

Currently this nodal approximate projection does not exist in a separate
operator like the MAC projection; instead we demonstrate below the steps needed
to compute the approximate projection.  This means we must compute explicitly the
right-hand-side , including the the divergence of the vector field, :math:`D(U^*)`,
solve the variable coefficient Poisson equation, then subtract the weighted
gradient term to make the vector field result satisfy the divergence constraint.

.. highlight:: c++

::
                  
   //
   // Given a cell-centered velocity (vel) field, a cell-centered
   // scalar field (sigma) field, and a source term S (either node-
   // or cell-centered )solve:
   //
   //   div( sigma * grad(phi) ) = div(vel) - S
   //
   // and then perform the projection:
   //
   //     vel = vel - sigma * grad(phi)
   // 

   //
   // Create the EB factory
   // 
   EBFArrayBoxFactory factory(eb_level, geom, grids, dmap, ng_ebs, ebs);

   //
   //  Create the cell-centered velocity field we want to project  
   //
   MultiFab vel(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);

   // Set velocity field to (1,0,0) including ghost cells for this example
   vel.setVal(1.0, 0, 1, 1);
   vel.setVal(0.0, 1, AMREX_SPACEDIM-1, 1);

   //
   // Setup linear operator, AKA the nodal Laplacian
   // 
   LPInfo lp_info;

   // If we want to use hypre to solve the full problem we do not need to coarsen the GMG stencils
   // if (use_hypre_as_full_solver)
   //     lp_info.setMaxCoarseningLevel(0);

   MLNodeLaplacian matrix({geom}, {grids}, {dmap}, lp_info,
                          Vector<EBFArrayBoxFactory const*>{&factory});

   // Set boundary conditions.
   // Here we use Neumann on the low x-face, Dirichlet on the high x-face,
   // and periodic in the other two directions
   // (the first argument is for the low end, the second is for the high end)
   // Note that Dirichlet boundary conditions are assumed to be homogeneous (i.e. phi = 0)
   matrix.setDomainBC({AMREX_D_DECL(LinOpBCType::Neumann,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)},
                      {AMREX_D_DECL(LinOpBCType::Dirichlet,
                                    LinOpBCType::Periodic,
                                    LinOpBCType::Periodic)});

   // Set matrix attributes to be used by MLMG solver
   matrix.setGaussSeidel(true);
   matrix.setHarmonicAverage(false);

   //
   // Compute RHS 
   //
   // NOTE: it's up to the user to compute the RHS. as opposed
   //       to the MAC projection case !!!
   //
   // NOTE: do this operation AFTER setting up the linear operator so
   //       that compRHS method can be used
   // 

   // RHS is nodal
   const BoxArray & nd_grids = amrex::convert(grids, IntVect{1,1,1}); // nodal grids

   // Multifab to host RHS
   MultiFab rhs(nd_grids, dmap, 1, 1, MFInfo(), factory);

   // Cell-centered contributions to RHS
   MultiFab S_cc(grids, dmap, 1, 1, MFInfo(), factory);
   S_cc.setVal(0.0); // Set it to zero for this example

   // Node-centered contributions to RHS
   MultiFab S_nd(nd_grids, dmap, 1, 1, MFInfo(), factory);
   S_nd.setVal(0.0); // Set it to zero for this example

   // Compute RHS -- vel must be cell-centered
   matrix.compRHS({&rhs}, {&vel}, {&S_nd}, {&S_cc});

   //
   // Create the cell-centered sigma field and set it to 1 for this example
   //
   MultiFab sigma(grids, dmap, 1, 1, MFInfo(), factory);
   sigma.setVal(1.0);

   // Set sigma
   matrix.setSigma(0, sigma);

   //
   // Create node-centered phi
   //
   MultiFab phi(nd_grids, dmap, 1, 1, MFInfo(), factory);
   phi.setVal(0.0);

   //
   // Setup MLMG solver
   //
   MLMG nodal_solver(matrix);

   // We can specify the maximum number of iterations
   nodal_solver.setMaxIter(mg_maxiter);
   nodal_solver.setCGMaxIter(mg_cg_maxiter);

   nodal_solver.setVerbose(mg_verbose);
   nodal_solver.setCGVerbose(mg_cg_verbose);

   // Set bottom-solver to use hypre instead of native BiCGStab 
   //   ( we could also have set this to cg, bicgcg, cgbicg)
   // if (use_hypre_as_full_solver || use_hypre_as_bottom_solver) 
   //     nodal_solver.setBottomSolver(MLMG::BottomSolver::hypre);

   // Define the relative tolerance
   Real reltol = 1.e-8;

   // Define the absolute tolerance; note that this argument is optional
   Real abstol = 1.e-15;

   //
   // Solve div( sigma * grad(phi) ) = RHS
   //
   nodal_solver.solve( {&phi}, {&rhs}, reltol, abstol);

   //
   // Create cell-centered multifab to hold value of -sigma*grad(phi) at cell-centers
   // 
   //
   MultiFab fluxes(grids, dmap, AMREX_SPACEDIM, 1, MFInfo(), factory);
   fluxes.setVal(0.0);

   // Get fluxes from solver
   nodal_solver.getFluxes( {&fluxes} );

   //
   // Apply projection explicitly --  vel = vel - sigma * grad(phi)  
   // 
   MultiFab::Add( *vel, *fluxes, 0, 0, AMREX_SPACEDIM, 0);

See ``Tutorials/LinearSolvers/Nodal_Projection_EB`` for the complete working example.

Tensor Solve
============

Application codes that solve the Navier-Stokes equations need to evaluate
the viscous term;  solving for this term implicitly requires a multi-component
solve with cross terms.  Because this is a commonly used motif, we provide
a tensor solve for cell-centered velocity components.

Consider a velocity field :math:`U = (u,v,w)` with all
components co-located on cell centers.  The viscous term can be written in vector form as

.. math::

 \nabla \cdot (\eta \nabla U) + \nabla \cdot (\eta (\nabla U)^T ) + \nabla \cdot ( (\kappa - \frac{2}{3} \eta) (\nabla \cdot U) )

and in 3-d Cartesian component form as

.. math::

 ( (\eta u_x)_x + (\eta u_y)_y + (\eta u_z)_z ) + ( (\eta u_x)_x + (\eta v_x)_y + (\eta w_x)_z ) +  ( (\kappa - \frac{2}{3} \eta) (u_x+v_y+w_z) )_x

 ( (\eta v_x)_x + (\eta v_y)_y + (\eta v_z)_z ) + ( (\eta u_y)_x + (\eta v_y)_y + (\eta w_y)_z ) +  ( (\kappa - \frac{2}{3} \eta) (u_x+v_y+w_z) )_y

 ( (\eta w_x)_x + (\eta w_y)_y + (\eta w_z)_z ) + ( (\eta u_z)_x + (\eta v_z)_y + (\eta w_z)_z ) +  ( (\kappa - \frac{2}{3} \eta) (u_x+v_y+w_z) )_z


Here :math:`\eta` is the dynamic viscosity and :math:`\kappa` is the bulk viscosity.  

We evaluate the following terms from the above using the ``MLABecLaplacian`` and ``MLEBABecLaplacian`` operators;

.. math::

   ( (\frac{4}{3} \eta + \kappa) u_x)_x + (              \eta           u_y)_y + (\eta u_z)_z 

                 (\eta           v_x)_x + ( (\frac{4}{3} \eta + \kappa) v_y)_y + (\eta v_z)_z 

    (\eta w_x)_x                        + (              \eta           w_y)_y + ( (\frac{4}{3} \eta + \kappa) w_z)_z 

the following cross-terms are evaluted separately using the ``MLTensorOp`` and ``MLEBTensorOp`` operators.

.. math::

    ( (\kappa - \frac{2}{3} \eta) (v_y + w_z) )_x + (\eta v_x)_y  + (\eta w_x)_z

    (\eta u_y)_x + ( (\kappa - \frac{2}{3} \eta) (u_x + w_z) )_y  + (\eta w_y)_z

    (\eta u_z)_x + (\eta v_z)_y - ( (\kappa - \frac{2}{3} \eta) (u_x + v_y) )_z

The code below is an example of how to set up the solver to compute the
viscous term `divtau` explicitly:

.. highlight:: c++

::

   Box domain(geom[0].Domain());

   // Set BCs for Poisson solver in bc_lo, bc_hi
   ...

   //
   // First define the operator "ebtensorop"
   // Note we call LPInfo().setMaxCoarseningLevel(0) because we are only applying the operator,
   //      not doing an implicit solve
   //
   //       (alpha * a - beta * (del dot b grad)) sol
   //
   // LPInfo                       info;
   MLEBTensorOp ebtensorop(geom, grids, dmap, LPInfo().setMaxCoarseningLevel(0),
                           amrex::GetVecOfConstPtrs(ebfactory));

   // It is essential that we set MaxOrder of the solver to 2
   // if we want to use the standard sol(i)-sol(i-1) approximation
   // for the gradient at Dirichlet boundaries.
   // The solver's default order is 3 and this uses three points for the
   // gradient at a Dirichlet boundary.
   ebtensorop.setMaxOrder(2);

   // LinOpBCType Definitions are in amrex/Src/Boundary/AMReX_LO_BCTYPES.H
   ebtensorop.setDomainBC ( {(LinOpBCType)bc_lo[0], (LinOpBCType)bc_lo[1], (LinOpBCType)bc_lo[2]},
                            {(LinOpBCType)bc_hi[0], (LinOpBCType)bc_hi[1], (LinOpBCType)bc_hi[2]} );

   // Return div (eta grad)) phi
   ebtensorop.setScalars(0.0, -1.0);

   amrex::Vector<amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>> b;
   b.resize(max_level + 1);

   // Compute the coefficients
   for (int lev = 0; lev < nlev; lev++)
   {
       // We average eta onto faces
       for(int dir = 0; dir < AMREX_SPACEDIM; dir++)
       {
           BoxArray edge_ba = grids[lev];
           edge_ba.surroundingNodes(dir);
           b[lev][dir].reset(new MultiFab(edge_ba, dmap[lev], 1, nghost, MFInfo(), *ebfactory[lev]));
       }

       average_cellcenter_to_face( GetArrOfPtrs(b[lev]), *etan[lev], geom[lev] );

       b[lev][0] -> FillBoundary(geom[lev].periodicity());
       b[lev][1] -> FillBoundary(geom[lev].periodicity());
       b[lev][2] -> FillBoundary(geom[lev].periodicity());

       ebtensorop.setShearViscosity  (lev, GetArrOfConstPtrs(b[lev]));
       ebtensorop.setEBShearViscosity(lev, (*eta[lev]));

       ebtensorop.setLevelBC ( lev, GetVecOfConstPtrs(vel)[lev] );
   }

   MLMG solver(ebtensorop);

   solver.apply(GetVecOfPtrs(divtau), GetVecOfPtrs(vel));


Multi-Component Operators
=========================

This section discusses solving linear systems in which the solution variable :math:`\mathbf{\phi}` has multiple components.
An example (implemented in the ``MultiComponent`` tutorial) might be:

.. math::

   D(\mathbf{\phi})_i = \sum_{i=1}^N \alpha_{ij} \nabla^2 \phi_j

(Note: only operators of the form :math:`D:\mathbb{R}^n\to\mathbb{R}^n` are currently allowed.)

- To implement a multi-component *cell-based* operator, inherit from the ``MLCellLinOp`` class.
  Override the ``getNComp`` function to return the number of components (``N``)that the operator will use.
  The solution and rhs fabs must also have at least one ghost node.
  ``Fapply``, ``Fsmooth``, ``Fflux`` must be implemented such that the solution and rhs fabs all have ``N`` components.

- Implementing a multi-component *node-based* operator is slightly different.
  A MC nodal operator must specify that the reflux-free coarse/fine strategy is being used by the solver.

  .. code::

     solver.setCFStrategy(MLMG::CFStrategy::ghostnodes);

  The reflux-free method circumvents the need to implement a special ``reflux`` at the coarse-fine boundary.
  This is accomplished by using ghost nodes.
  Each AMR level must have 2 layers of ghost nodes.
  The second (outermost) layer of nodes is treated as constant by the relaxation, essentially acting as a Dirichlet boundary.
  The first layer of nodes is evolved using the relaxation, in the same manner as the rest of the solution.
  When the residual is restricted onto the coarse level (in ``reflux``) this allows the residual at the coarse-fine boundary to be interpolated using the first layer of ghost nodes.
  :numref:`fig::refluxfreecoarsefine` illustrates the how the coarse-fine update takes place.

  .. _fig::refluxfreecoarsefine:

  .. figure:: ./LinearSolvers/refluxfreecoarsefine.png
	      :height: 2cm
	      :align: center

	      : Reflux-free coarse-fine boundary update.
	      Level 2 ghost nodes (small dark blue) are interpolated from coarse boundary.
	      Level 1 ghost nodes are updated during the relaxation along with all the other interior fine nodes.
	      Coarse nodes (large blue) on the coarse/fine boundary are updated by restricting with interior nodes
	      and the first level of ghost nodes.
	      Coarse nodes underneath level 2 ghost nodes are not updated.
	      The remaining coarse nodes are updates by restriction.
	      
  The MC nodal operator can inherit from the ``MCNodeLinOp`` class.
  ``Fapply``, ``Fsmooth``, and ``Fflux`` must update level 1 ghost nodes that are inside the domain.
  `interpolation` and `restriction` can be implemented as usual.
  `reflux` is a straightforward restriction from fine to coarse, using level 1 ghost nodes for restriction as described above.
  
  See ``Tutorials/LinearSolvers/MultiComponent`` for a complete working example.

   

.. solver reuse

