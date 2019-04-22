
module amrex_multigrid_module

  use iso_c_binding
  use amrex_base_module
  use amrex_linop_module
  implicit none

  integer, parameter, public :: amrex_bottom_smoother = 0
  integer, parameter, public :: amrex_bottom_bicgstab = 1
  integer, parameter, public :: amrex_bottom_cg       = 2
  integer, parameter, public :: amrex_bottom_hypre    = 3
  integer, parameter, public :: amrex_bottom_petsc    = 4
  integer, parameter, public :: amrex_bottom_default  = 1

  private
  public :: amrex_multigrid_build, amrex_multigrid_destroy

  type, public :: amrex_multigrid
     logical     :: owner = .false.
     type(c_ptr) :: p = c_null_ptr
   contains
     generic :: assignment(=) => amrex_multigrid_assign  ! shallow copy
     procedure :: solve                 => amrex_multigrid_solve
     procedure :: get_grad_solution     => amrex_multigrid_get_grad_solution
     procedure :: get_fluxes            => amrex_multigrid_get_fluxes
     procedure :: comp_residual         => amrex_multigrid_comp_residual
     procedure :: set_verbose           => amrex_multigrid_set_verbose
     procedure :: set_max_iter          => amrex_multigrid_set_max_iter
     procedure :: set_max_fmg_iter      => amrex_multigrid_set_max_fmg_iter
     procedure :: set_bottom_solver     => amrex_multigrid_set_bottom_solver
     procedure :: set_cg_verbose        => amrex_multigrid_set_cg_verbose
     procedure :: set_always_use_bnorm  => amrex_multigrid_set_always_use_bnorm
     procedure :: set_final_fill_bc     => amrex_multigrid_set_final_fill_bc
     procedure, private :: amrex_multigrid_assign
     procedure, private :: amrex_multigrid_solve
     procedure, private :: amrex_multigrid_get_grad_solution
     procedure, private :: amrex_multigrid_get_fluxes
     procedure, private :: amrex_multigrid_comp_residual
     procedure, private :: amrex_multigrid_set_verbose
     procedure, private :: amrex_multigrid_set_max_iter
     procedure, private :: amrex_multigrid_set_max_fmg_iter
     procedure, private :: amrex_multigrid_set_cg_verbose
     procedure, private :: amrex_multigrid_set_always_use_bnorm
     procedure, private :: amrex_multigrid_set_final_fill_bc

#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_multigrid_destroy
#endif
  end type amrex_multigrid

  ! interfaces to C++ functions

  interface
     subroutine amrex_fi_new_multigrid (mg, linop) bind(c)
       import
       implicit none
       type(c_ptr), intent(inout) :: mg
       type(c_ptr), value :: linop
     end subroutine amrex_fi_new_multigrid
  
     subroutine amrex_fi_delete_multigrid (mg) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
     end subroutine amrex_fi_delete_multigrid

     function amrex_fi_multigrid_solve (mg, sol, rhs, tol_rel, tol_abs) result(r) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       type(c_ptr), intent(in) :: sol(*), rhs(*)
       real(amrex_real), value :: tol_rel, tol_abs
       real(amrex_real) :: r
     end function amrex_fi_multigrid_solve

     subroutine amrex_fi_multigrid_get_grad_solution (mg, gradsol) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       type(c_ptr), intent(in) :: gradsol(*)
     end subroutine amrex_fi_multigrid_get_grad_solution

     subroutine amrex_fi_multigrid_get_fluxes (mg, flx) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       type(c_ptr), intent(in) :: flx(*)
     end subroutine amrex_fi_multigrid_get_fluxes

     subroutine amrex_fi_multigrid_comp_residual (mg, res, sol, rhs) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       type(c_ptr), intent(in) :: res(*), sol(*), rhs(*)
     end subroutine amrex_fi_multigrid_comp_residual

     subroutine amrex_fi_multigrid_set_verbose (mg, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       integer(c_int), value :: v
     end subroutine amrex_fi_multigrid_set_verbose

     subroutine amrex_fi_multigrid_set_max_iter (mg, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       integer(c_int), value :: n
     end subroutine amrex_fi_multigrid_set_max_iter

     subroutine amrex_fi_multigrid_set_max_fmg_iter (mg, n) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       integer(c_int), value :: n
     end subroutine amrex_fi_multigrid_set_max_fmg_iter

     subroutine amrex_fi_multigrid_set_bottom_solver (mg, s) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       integer(c_int), value :: s
     end subroutine amrex_fi_multigrid_set_bottom_solver

     subroutine amrex_fi_multigrid_set_cg_verbose (mg, v) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       integer(c_int), value :: v
     end subroutine amrex_fi_multigrid_set_cg_verbose

     subroutine amrex_fi_multigrid_set_always_use_bnorm (mg, f) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       integer(c_int), value :: f
     end subroutine amrex_fi_multigrid_set_always_use_bnorm

     subroutine amrex_fi_multigrid_set_final_fill_bc (mg, f) bind(c)
       import
       implicit none
       type(c_ptr), value :: mg
       integer(c_int), value :: f
     end subroutine amrex_fi_multigrid_set_final_fill_bc
  end interface

contains

  subroutine amrex_multigrid_assign (dst, src)
    class(amrex_multigrid), intent(inout) :: dst
    type(amrex_multigrid), intent(in) :: src
    call amrex_multigrid_destroy(dst)
    dst%owner = .false.
    dst%p = src%p
  end subroutine amrex_multigrid_assign


  subroutine amrex_multigrid_build (mg, linop)
    type(amrex_multigrid), intent(inout) :: mg
    class(amrex_linop), intent(in) :: linop
    mg%owner = .true.
    call amrex_fi_new_multigrid(mg%p, linop%p)
  end subroutine amrex_multigrid_build


  subroutine amrex_multigrid_destroy (this)
    type(amrex_multigrid), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_multigrid(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_multigrid_destroy


  function amrex_multigrid_solve (mg, sol, rhs, tol_rel, tol_abs) result(r)
    class(amrex_multigrid), intent(inout) :: mg
    type(amrex_multifab) :: sol(0:)
    type(amrex_multifab), intent(in) :: rhs(0:)
    real(amrex_real), intent(in) :: tol_rel, tol_abs
    real(amrex_real) :: r
    type(c_ptr) :: psol(0:size(rhs)-1), prhs(0:size(rhs)-1)
    integer :: i
    do i = 0, size(rhs)-1
       psol(i) = sol(i)%p
       prhs(i) = rhs(i)%p
    end do
    r = amrex_fi_multigrid_solve(mg%p, psol, prhs, tol_rel, tol_abs)
  end function amrex_multigrid_solve


  subroutine amrex_multigrid_get_grad_solution (mg, gradsol)
    class(amrex_multigrid), intent(inout) :: mg
    type(amrex_multifab), intent(inout) :: gradsol(1:,0:)
    type(c_ptr) :: pg(0:size(gradsol)-1)
    integer :: i, idim, ilev
    i = 0
    do ilev = 0, size(gradsol,2)-1
       do idim = 1, amrex_spacedim
          pg(i) = gradsol(idim,ilev)%p
          i = i+1
       end do
    end do
    call amrex_fi_multigrid_get_grad_solution(mg%p, pg)
  end subroutine amrex_multigrid_get_grad_solution


  subroutine amrex_multigrid_get_fluxes (mg, flx)
    class(amrex_multigrid), intent(inout) :: mg
    type(amrex_multifab), intent(inout) :: flx(1:,0:)
    type(c_ptr) :: pf(0:size(flx)-1)
    integer :: i, idim, ilev
    i = 0
    do ilev = 0, size(flx,2)-1
       do idim = 1, amrex_spacedim
          pf(i) = flx(idim,ilev)%p
          i = i+1
       end do
    end do
    call amrex_fi_multigrid_get_fluxes(mg%p, pf)
  end subroutine amrex_multigrid_get_fluxes


  subroutine amrex_multigrid_comp_residual (mg, res, sol, rhs)
    class(amrex_multigrid), intent(inout) :: mg
    type(amrex_multifab), intent(inout) :: res(0:)
    type(amrex_multifab), intent(in) :: sol(0:), rhs(0:)
    type(c_ptr) :: pres(0:size(res)-1), psol(0:size(res)-1), prhs(0:size(res)-1)
    integer :: i
    do i = 0, size(res)-1
       pres(i) = res(i)%p
       psol(i) = sol(i)%p
       prhs(i) = rhs(i)%p
    end do
    call amrex_fi_multigrid_comp_residual(mg%p, pres, psol, prhs)
  end subroutine amrex_multigrid_comp_residual

  
  subroutine amrex_multigrid_set_verbose (mg, v)
    class(amrex_multigrid), intent(inout) :: mg
    integer, intent(in) :: v
    call amrex_fi_multigrid_set_verbose(mg%p, v)
  end subroutine amrex_multigrid_set_verbose


  subroutine amrex_multigrid_set_max_iter (mg, n)
    class(amrex_multigrid), intent(inout) :: mg
    integer, intent(in) :: n
    call amrex_fi_multigrid_set_max_iter(mg%p, n)
  end subroutine amrex_multigrid_set_max_iter


  subroutine amrex_multigrid_set_max_fmg_iter (mg, n)
    class(amrex_multigrid), intent(inout) :: mg
    integer, intent(in) :: n
    call amrex_fi_multigrid_set_max_fmg_iter(mg%p, n)
  end subroutine amrex_multigrid_set_max_fmg_iter


  subroutine amrex_multigrid_set_bottom_solver (mg, s)
    class(amrex_multigrid), intent(inout) :: mg
    integer, intent(in) :: s
    call amrex_fi_multigrid_set_bottom_solver(mg%p, s)
  end subroutine amrex_multigrid_set_bottom_solver


  subroutine amrex_multigrid_set_cg_verbose (mg, v)
    class(amrex_multigrid), intent(inout) :: mg
    integer, intent(in) :: v
    call amrex_fi_multigrid_set_cg_verbose(mg%p, v)
  end subroutine amrex_multigrid_set_cg_verbose


  subroutine amrex_multigrid_set_always_use_bnorm (mg, f)
    class(amrex_multigrid), intent(inout) :: mg
    logical, intent(in) :: f
    integer :: iflag
    if (f) then
       iflag = 1
    else
       iflag = 0
    end if
    call amrex_fi_multigrid_set_always_use_bnorm(mg%p, iflag)
  end subroutine amrex_multigrid_set_always_use_bnorm


  subroutine amrex_multigrid_set_final_fill_bc (mg, f)
    class(amrex_multigrid), intent(inout) :: mg
    logical, intent(in) :: f
    integer :: iflag
    if (f) then
       iflag = 1
    else
       iflag = 0
    end if
    call amrex_fi_multigrid_set_final_fill_bc(mg%p, iflag)
  end subroutine amrex_multigrid_set_final_fill_bc

end module amrex_multigrid_module
