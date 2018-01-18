.. _Chap:Boundary:

Boundary Conditions
===================

This chapter describes how to implement domain boundary conditions in .
A ghost cell that is outside of the valid region can be thought of as either
“interior” (for periodic and coarse-fine ghost cells), or “physical”.
Physical boundary conditions can include inflow, outflow, slip/no-slip walls,
but are ultimately linked to mathematical Dirichlet or Neumann conditions.

The basic idea behind physical boundary conditions is as follows:

-  Create a BCRec object, which is essentially a multidimensional integer array of
   2*DIM components. Each component defines a boundary condition type for
   the lo/hi side of the domain, for each direction.
   See Src/Base/AMReX_BC_TYPES.H for common physical and mathematical types.
   If there is more than one variable, we can create an array of BCRec objects,
   and pass in a pointer to the 0-index component since the arrays for all the
   components are contiguous in memory.
   Here we need to provide boundary types to each component of the
   MultiFab. Below is an example of setting up Vector<BCRec>
   before the call to ghost cell routines.

   ::

         // Set up BC; see Src/Base/AMReX_BC_TYPES.H for supported types
         Vector<BCRec> bc(phi.nComp());
         for (int n = 0; n < phi.nComp(); ++n)
         {
             for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
             {
                 if (Geometry::isPeriodic(idim))
                 {
                     bc[n].setLo(idim, BCType::int_dir); // interior
                     bc[n].setHi(idim, BCType::int_dir);
                 }
                 else
                 {
                     bc[n].setLo(idim, BCType::foextrap); // first-order extrapolation
                     bc[n].setHi(idim, BCType::foextrap);
                 }
             }
         }

   amrex::BCType has the following types,

       int_dir
           Interior, including periodic boundary

       ext_dir
           “External Dirichlet”. It is the user’s responsibility to write a routine
           to fill ghost cells (more details below).

       foextrap
           “First Order Extrapolation”
           First order extrapolation from last cell in interior.

       reflect_even
           Reflection from interior cells with sign
           unchanged, :math:`q(-i) = q(i)`.

       reflect_odd
           Reflection from interior cells with sign
           unchanged, :math:`q(-i) = -q(i)`.

-  We have interfaces to a fortran routine that fills ghost cells at domain
   boundaries based on the boundary condition type defined in the BCRec object.
   It is the user’s responsibility to have a consisent definition of what the ghost cells
   represent. A common option used in  codes is to fill the domain ghost cells
   with the value that lies on the boundary (as opposed to another common option where
   the value in the ghost cell represents an extrapolated value based on the boundary
   condition type). Then in our stencil based “work” codes, we also pass in the
   BCRec object and use modified stencils near the domain boundary that know the value
   in the first ghost cell represents the value on the boundary.

Depending on the level of complexity of your code, there are various options
for filling domain boundary ghost cells.

For single-level codes built from Src/Base (excluding the
Src/AmrCore and Src/Amr source code directories), you will have
single-level MultiFabs filled with data in the valid region where you need
to fill the ghost cells on each grid. There are essentially three ways to fill the ghost
cells. (refer to Tutorials/Basic/HeatEquation_EX2_C for an example).

::

    MultiFab mf;
    Geometry geom;
    Vector<BCRec> bc;

    // ...

    // fills interior and periodic domain boundary ghost cells
    mf.FillBoundary(geom.periodicity());

    // fills interior (but not periodic domain boundary) ghost cells
    mf.FillBoundary();

    // fills physical domain boundary ghost cells
    FillDomainBoundary(mf, geom, bc);

FillDomainBoundary() is a function is in Src/Base/AMReX_BCUtil.cpp,
and is essentially an interface to fortran subroutine amrex_fab_filcc()
in Src/Base/AMReX_filcc_mod.F90, which ultimately calls fortran
subroutine filcc() in Src/Base/AMReX_FILCC_XD.F. To create more
custom boundary conditions, create a local modified copy of
Src/Base/AMReX_FILCC_XD.F and put it your local source code.

For multi-level codes using the Src_AmrCore source code, the
functions described above still work, however additional classes need to
be set up since the FillPatch routines call them.
In fact it is possible to avoid using the single-level calls directly if
you fill all your grids and ghost cells using the FillPatch routines.
Refer to Tutorials/Amr/Advection_AmrCore/ for an example.
The class PhysBCFunct in Src/Base/AMReX_PhysBCFunct.cpp
is derived from PhysBCFunctBase and contains a BCRec, Geometry,
and a pointer to a BndryFunctBase function.

Note that PhyBCFunct is an example of how to derive from PhysBCFunctBase and is
not meant to be a base class. PhysBCFunctBase is the base class.
PhysBCFunctBase is designed for users to derive and extend.
You could/should write your own class derived from PhysBCFuncBase.
There you can make modifications such as storing a vector of BCRecs for, e.g.,
multiple component MultiFabs.

The function FillBoundary fills physical ghost cells and has a similar functionality
to the single-level case described above, where FillDomainBoundary
fills the physical ghost cells. In fact you can have your BndryFunctBase
point to the same filcc routines called by the single-level routines.
