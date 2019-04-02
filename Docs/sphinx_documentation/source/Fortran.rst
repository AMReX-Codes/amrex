.. role:: cpp(code)
   :language: c++

.. role:: fortran(code)
   :language: fortran


Getting Started
===============

We have discussed AMReX’s build systems in the chapter on
:ref:`Chap:BuildingAMReX`.  To build with GNU Make, we need to include the
Fortran interface source tree into the make system. The source codes for the
Fortran interface are in ``amrex/Src/F_Interfaces`` and there are several
sub-directories. The "Base" directory includes sources for the basic
functionality, the "AmrCore" directory wraps around the :cpp:`AmrCore` class
(see the chapter on :ref:`Chap:AmrCore`), and the "Octree" directory adds
support for octree type of AMR grids. Each directory has a "Make.package" file
that can be included in make files (see ``amrex/Tutorials/Basic/HelloWorld_F`` and
``amrex/Tutorials/Amr/Advection_F`` for examples). The libamrex approach includes the
Fortran interface by default. 

A simple example can be found at ``amrex/Tutorials/Basic/HelloWorld_F/``. The source code
is shown below in its entirety.

.. highlight:: fortran

::

    program main
      use amrex_base_module
      implicit none
      call amrex_init()
      if (amrex_parallel_ioprocessor()) then
         print *, "Hello world!"
      end if
      call amrex_finalize()
    end program main

To access the AMReX Fortran interfaces, we can use these three
modules, amrex_base_module for the basics functionalities
(Section `2 <#sec:fi:basics>`__), amrex_amrcore_module for AMR
support (Section `3 <#sec:fi:amrcore>`__) and amrex_octree_module
for octree style AMR (Section `4 <#sec:fi:octree>`__).

.. _sec:fi:basics:

The Basics
==========

Module :fortran:`amrex_base_module` is a collection of various Fortran modules
providing interfaces to most of the basics of AMReX C++ library (see the
chapter on :ref:`Chap:Basics`). These modules shown in this section can be used
without being explicitly included because they are included by
:fortran:`amrex_base_module`.

The spatial dimension is an integer parameter :fortran:`amrex_spacedim`.  We
can also use the :fortran:`AMREX_SPACEDIM` macro in preprocessed Fortran codes
(e.g., .F90 files) just like in the C++ codes. Unlike in C++, the convention
for AMReX Fortran interface is that coordinate direction index starts at 1.

There is an integer parameter :fortran:`amrex_real`, a Fortran kind parameter
for :fortran:`real`. Fortran :fortran:`real(amrex_real)` corresponds to
:cpp:`amrex::Real` in C++, which is either double or single precision depending
the setting of precision.

The module :fortran:`amrex_parallel_module` (
``amrex/Src/F_Interfaces/Base/AMReX_parallel_mod.F90``) includes wrappers to the
:cpp:`ParallelDescriptor` namespace, which is in turn a wrapper to the parallel
communication library used by AMReX (e.g. MPI).

The module :cpp:`amrex_parmparse_module` (
``amrex/Src/Base/AMReX_parmparse_mod.F90``) provides interface to
:cpp:`ParmParse` (see the section on :ref:`sec:basics:parmparse`). Here are some
examples.

.. highlight:: fortran

::

      type(amrex_parmparse) :: pp
      integer :: n_cell, max_grid_size
      call amrex_parmparse_build(pp)
      call pp%get("n_cell", n_cell)
      max_grid_size = 32 ! default size
      call pp%query("max_grid_size", max_grid_size)
      call amrex_parmpase_destroy(pp) ! optional if compiler supports finalization

Finalization is a Fortran 2003 feature that some compilers may not support. For
those compilers, we must explicitly destroy the objects, otherwise there will
be memory leaks. This applies to many other derived types.

:fortran:`amrex_box` is a derived type in :fortran:`amrex_box_module`
``amrex/Src/F_Interfaces/Base/AMReX_box_mod.F90``. It has three members, :fortran:`lo`
(lower corner), :fortran:`hi` (upper corner) and :fortran:`nodal` (logical flag
for index type).

:fortran:`amrex_geometry` is a wrapper for the :cpp:`Geometry` class
containing information for the physical domain. Below is an example
of building it.

.. highlight:: fortran

::

      integer :: n_cell
      type(amrex_box) :: domain
      type(amrex_geometry) : geom
      ! n_cell = ...
      ! Define a single box covering the domain
      domain = amrex_box((/0,0,0/), (/n_cell-1, n_cell-1, n_cell-1/))
      ! This defines a amrex_geometry object.
      call amrex_geometry_build(geom, domain)
      !
      ! ...
      !
      call amrex_geometry_destroy(geom)

:fortran:`amrex_boxarray` ( ``amrex/Src/F_Interfaces/Base/AMReX_boxarray_mod.F90``) is a
wrapper for the :cpp:`BoxArray` class, and :fortran:`amrex_distromap` (
``amrex/Src/F_Interfaces/Base/AMReX_distromap_mod.F90``) is a wrapper for the
:cpp:`DistributionMapping` class. Here is an example of building a
:cpp:`BoxArray` and a :cpp:`DistributionMapping`.

.. highlight:: fortran

::

      integer :: n_cell
      type(amrex_box) :: domain
      type(amrex_boxarray) : ba
      type(amrex_distromap) :: dm
      ! n_cell = ...
      ! Define a single box covering the domain
      domain = amrex_box((/0,0,0/), (/n_cell-1, n_cell-1, n_cell-1/))
      ! Initialize the boxarray "ba" from the single box "bx"
      call amrex_boxarray_build(ba, domain)
      ! Break up boxarray "ba" into chunks no larger than "max_grid_size"
      call ba%maxSize(max_grid_size)
      ! Build a DistributionMapping for the boxarray
      call amrex_distromap_build(dm, ba)
      !
      ! ...
      !
      call amrex_distromap_distromap(dm)
      call amrex_boxarray_destroy(ba)

Given :fortran:`amrex_boxarray` and :fortran:`amrex_distromap`, we can build
:cpp:`amrex_multifab`, a wrapper for the :cpp:`MultiFab` class, as follows.

.. highlight:: fortran

::

      integer :: ncomp, nghost
      type(amrex_boxarray) : ba
      type(amrex_distromap) :: dm
      type(amrex_multifab) :: mf, ndmf
      ! Build amrex_boxarray and amrex_distromap
      ! ncomp = ...
      ! nghost = ...
      ! ...
      ! Build amrex_multifab with ncomp component and nghost ghost cells
      call amrex_multifab_build(mf, ba, dm, ncomp, nghost)
      ! Build a nodal multifab
      call amrex_multifab_build(ndmf,ba,dm,ncomp,nghost,(/.true.,.true.,.true./))
      !
      ! ...
      !
      call amrex_multifab_destroy(mf)
      call amrex_multifab_destroy(ndmf)

There are many type-bound procedures for :fortran:`amrex_multifab`. For example

::

      ncomp   ! Return the number of components
      nghost  ! Return the number of ghost cells
      setval  ! Set the data to the given value 
      copy    ! Copy data from given amrex_multifab to this amrex_multifab

Note that the copy function here only works on copying data from another
:fortran:`amrex_multifab` built with the same :fortran:`amrex_distromap`, like
the :cpp:`MultiFab::Copy` function in C++.  :fortran:`amrex_multifab` also has
two parallel communication procedures, :fortran:`fill_boundary` and
:fortran:`parallel_copy`. Their and interface and usage are very similar to
functions :cpp:`FillBoundary` and :cpp:`ParallelCopy` for :cpp:`MultiFab` in
C++.

.. highlight:: fortran

::

      type(amrex_geometry) :: geom
      type(amrex_multifab) :: mf, mfsrc
      ! ...
      call mf%fill_boundary(geom)       ! Fill all components
      call mf%fill_boundary(geom, 1, 3) ! Fill 3 components starting with component 1

      call mf%parallel_copy(mfsrc, geom) ! Parallel copy from another multifab

It should be emphasized that the component index for :fortran:`amrex_multifab`
starts with 1 following Fortran convention. This is different from the C++ part
of AMReX.

AMReX provides a Fortran interface to :fortran:`MFIter` for iterating over the
data in :fortran:`amrex_multifab`. The Fortran type for this is
:fortran:`amrex_mfiter`. Here is an example of using :fortran:`amrex_mfiter` to
loop over :fortran:`amrex_multifab` with tiling and launch a kernel function.

.. highlight:: fortran

::

      integer :: plo(4), phi(4)
      type(amrex_box) :: bx
      real(amrex_real), contiguous, dimension(:,:,:,:), pointer :: po, pn
      type(amrex_multifab) :: old_phi, new_phi
      type(amrex_mfiter) :: mfi
      ! Define old_phi and new_phi ...
      ! In this example they are built with the same boxarray and distromap.
      ! And they have the same number of ghost cells and 1 component.
      call amrex_mfiter_build(mfi, old_phi, tiling=.true.)
      do while (mfi%next())
        bx = mfi%tilebox()
        po => old_phi%dataptr(mfi)
        pn => new_phi%dataptr(mfi)
        plo = lbound(po)
        phi = ubound(po)
        call update_phi(bx%lo, bx&hi, po, pn, plo, phi)
      end do
      call amrex_mfiter_destroy(mfi)

Here procedure :fortran:`update_phi` is

::

     subroutine update_phi (lo, hi, pold, pnew, plo, phi)
      integer, intent(in) :: lo(3), hi(3), plo(3), phi(3)
       real(amrex_real),intent(in   ) pold(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
       real(amrex_real),intent(inout) pnew(plo(1):phi(1),plo(2):phi(2),plo(3):phi(3))
       ! ...
     end subroutine update_phi

Note that amrex_multifab’s procedure :fortran:`dataptr` takes
:fortran:`amrex_mfiter` and returns a 4-dimensional Fortran pointer. For
performance, we should declare the pointer as :fortran:`contiguous`. In C++,
the similar operation returns a reference to :cpp:`FArrayBox`.  However,
:cpp:`FArrayBox` and Fortran pointer have a similar capability of containing
array bound information. We can call :fortran:`lbound` and :fortran:`ubound` on
the pointer to return its lower and upper bounds. The first three dimensions of
the bounds are spatial and the fourth is for the number of component.

Many of the derived Fortran types in  (e.g., :fortran:`amrex_multifab`,
:fortran:`amrex_boxarray`, :fortran:`amrex_distromap`, :fortran:`amrex_mfiter`,
and :fortran:`amrex_geometry`) contain a :fortran:`type(c_ptr)` that points a
C++ object. They also contain a :fortran:`logical` type indicating whether or
not this object owns the underlying object (i.e., responsible for deleting the
object). Due to the semantics of Fortran, one should not return these types
with functions. Instead we should pass them as arguments to procedures
(preferably with :fortran:`intent` specified). These five types all have
assignment(=) operator that performs a shallow copy. After the assignment, the
original objects still owns the data and the copy is just an alias. For
example,

.. highlight:: fortran

::

      type(amrex_multifab) :: mf1, mf2
      call amrex_multifab_build(mf1, ...)
      call amrex_multifab_build(mf2, ...)
      ! At this point, both mf1 and mf2 are data owners
      mf2 = mf1   ! This will destroy the original data in mf2.
                  ! Then mf2 becomes a shallow copy of mf1.
                  ! mf1 is still the owner of the data.
      call amrex_multifab_destroy(mf1)
      ! mf2 no longer contains a valid pointer because mf1 has been destroyed. 
      call amrex_multifab_destroyed(mf2)  ! But we still need to destroy it.

If we need to transfer the ownership, :fortran:`amrex_multifab`,
:fortran:`amrex_boxarray` and :fortran:`amrex_distromap` provide type-bound
:fortran:`move` procedure. We can use it as follows

.. highlight:: fortran

::

      type(amrex_multifab) :: mf1, mf2
      call amrex_multifab_build(mf1, ...)
      call mf2%move(mf1)   ! mf2 is now the data owner and mf1 is not.
      call amrex_multifab_destroy(mf1)
      call amrex_multifab_destroyed(mf2)

:fortran:`amrex_multifab` also has a type-bound :fortran:`swap` procedure for
exchanging the data.

AMReX also provides :fortran:`amrex_plotfile_module` for writing plotfiles. The
interface is similar to the C++ versions.


.. _sec:fi:amrcore:

Amr Core Infrastructure
=======================

The module :fortran:`amrex_amr_module` provides interfaces to AMR core
infrastructure. With AMR, the main program might look like below,

.. highlight:: fortran

::

      program main
        use amrex_amr_module
        implicit none  
        call amrex_init()
        call amrex_amrcore_init()
        call my_amr_init()       ! user's own code, not part of AMReX
        ! ...
        call my_amr_finalize()   ! user's own code, not part of AMReX
        call amrex_amrcore_finalize()
        call amrex_finalize()
      end program main

Here we need to call :fortran:`amrex_amrcore_init` and
:fortran:`amrex_amrcore_finalize`. And usually we need to call application code
specific procedures to provide some “hooks” needed by AMReX.  In C++, this is
achieved by using virtual functions. In Fortran, we need to call

.. highlight:: fortran

::

      subroutine amrex_init_virtual_functions (mk_lev_scrtch, mk_lev_crse, &
                                               mk_lev_re, clr_lev, err_est)

        ! Make a new level from scratch using provided boxarray and distromap
        ! Only used during initialization.
        procedure(amrex_make_level_proc)  :: mk_lev_scrtch
        ! Make a new level using provided boxarray and distromap, and fill
        ! with interpolated coarse level data.
        procedure(amrex_make_level_proc)  :: mk_lev_crse
        ! Remake an existing level using provided boxarray and distromap,
        ! and fill with existing fine and coarse data.
        procedure(amrex_make_level_proc)  :: mk_lev_re
        ! Delete level data
        procedure(amrex_clear_level_proc) :: clr_lev
        ! Tag cells for refinement
        procedure(amrex_error_est_proc)   :: err_est
      end subroutine amrex_init_virtual_functions

We need to provide five functions and these functions have three types of
interfaces:

.. highlight:: fortran

::

      subroutine amrex_make_level_proc (lev, time, ba, dm) bind(c)
        import
        implicit none
        integer, intent(in), value :: lev
        real(amrex_real), intent(in), value :: time
        type(c_ptr), intent(in), value :: ba, dm
      end subroutine amrex_make_level_proc
      
      subroutine amrex_clear_level_proc (lev) bind(c)
        import
        implicit none
        integer, intent(in) , value :: lev
      end subroutine amrex_clear_level_proc
      
      subroutine amrex_error_est_proc (lev, tags, time, tagval, clearval) bind(c)
        import
        implicit none
        integer, intent(in), value :: lev
        type(c_ptr), intent(in), value :: tags
        real(amrex_real), intent(in), value :: time
        character(c_char), intent(in), value :: tagval, clearval
      end subroutine amrex_error_est_proc

Tutorials/Amr/Advection_F/Source/my_amr_mod.F90 shows an
example of the setup process. The user provided
:fortran:`procedure(amrex_error_est_proc)` has a tags argument that
is of type :fortran:`c_ptr` and its value is a pointer to a  
:fortran:`TagBoxArray` object. We need to convert this into a Fortran
:fortran:`amrex_tagboxarray` object.

::

      type(amrex_tagboxarray) :: tag
      tag = tags

The module :fortran:`amrex_fillpatch_module` provides interface to
C++ functions :cpp:`FillPatchSinglelevel` and :cpp:`FillPatchTwoLevels`. To use
it, the application code needs to provide procedures for interpolation and
filling physical boundaries.  See
Tutorials/Amr/Advection_F/Source/fillpatch_mod.F90 for an example.

Module :fortran:`amrex_fluxregister_module` provides interface to
:cpp:`FluxRegister` (see the section on :ref:`sec:amrcore:fluxreg`). Its usage
is demonstrated in the tutorial at Tutorials/Amr/Advection_F/.


.. _sec:fi:octree:

Octree
======

In AMReX, the union of fine level grids is properly contained within the union
of coarse level grids. There are no required direct parent-child connections
between levels. Therefore, grids in AMReX in general cannot be represented by
trees. Nevertheless, octree type grids are supported via Fortran interface,
because  grids are more general than octree grids. A tutorial example using
amrex_octree_module ( ``amrex/Src/F_Interfaces/Octree/AMReX_octree_mod.f90``) is
available at ``amrex/Tutorials/Amr/Advection_octree_F/``. Procedures
:fortran:`amrex_octree_init` and :fortran:`amrex_octree_finalize` must be
called as follows,

.. highlight:: fortran

::

      program main
        use amrex_amrcore_module
        use amrex_octree_module
        implicit none
        call amrex_init()
        call amrex_octree_int()  ! This should be called before amrex_amrcore_init.
        call amrex_amrcore_init()
        call my_amr_init()       ! user's own code, not part of AMReX
        ! ...
        call my_amr_finalize()   ! user's own code, not part of AMReX
        call amrex_amrcore_finalize()
        call amrex_octree_finalize()
        call amrex_finalize()
      end program main

By default, the grid size is :math:`8^3`, and this can be changed via
:cpp:`ParmParse` parameter ``amr.max_grid_size``. The module
:fortran:`amrex_octree_module` provides :fortran:`amrex_octree_iter` that can
be used to iterate over leaves of octree. For example,

.. highlight:: fortran

::

      type(amrex_octree_iter) :: oti
      type(multifab) :: phi_new(*)   ! one multifab for each level
      integer :: ilev, igrd
      type(amrex_box) :: bx
      real(amrex_real), contiguous, pointer, dimension(:,:,:,:) :: pout
      call amrex_octree_iter_build(oti)
      do while(oti%next())
         ilev = oti%level()
         igrd = oti%grid_index()
         bx   = oti%box()
         pout => phi_new(ilev)%dataptr(igrd)
         ! ...
      end do
      call amrex_octree_iter_destroy(oti)
