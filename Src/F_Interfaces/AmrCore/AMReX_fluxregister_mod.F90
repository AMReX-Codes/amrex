
module amrex_fluxregister_module

  use iso_c_binding
  use amrex_base_module
  
  implicit none

  private

  public :: amrex_fluxregister_destroy ! List first to avoid XL compiler bug
  public :: amrex_fluxregister_build

  type, public :: amrex_fluxregister
     logical     :: owner = .false.
     type(c_ptr) :: p     = c_null_ptr
     integer     :: flev  = -1
   contains
     generic   :: assignment(=) => amrex_fluxregister_assign   ! shallow copy
     generic   :: fineadd       => amrex_fluxregister_fineadd, amrex_fluxregister_fineadd_1fab
     procedure :: crseinit      => amrex_fluxregister_crseinit
     procedure :: crseadd       => amrex_fluxregister_crseadd
     procedure :: setval        => amrex_fluxregister_setval
     procedure :: reflux        => amrex_fluxregister_reflux
     procedure :: overwrite     => amrex_fluxregister_overwrite
     procedure, private :: amrex_fluxregister_assign
     procedure, private :: amrex_fluxregister_fineadd
     procedure, private :: amrex_fluxregister_fineadd_1fab
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_fluxregister_destroy
#endif
  end type amrex_fluxregister

  interface
     subroutine amrex_fi_new_fluxregister (fr, ba, dm, rr, flev, nc) bind(c)
       import
       implicit none
       type(c_ptr) :: fr
       type(c_ptr), value :: ba, dm
       integer, value :: rr, flev, nc
     end subroutine amrex_fi_new_fluxregister

     subroutine amrex_fi_delete_fluxregister (fr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
     end subroutine amrex_fi_delete_fluxregister

     subroutine amrex_fi_fluxregister_fineadd (fr, flxs, scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       type(c_ptr), intent(in) :: flxs(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_fineadd

     subroutine amrex_fi_fluxregister_fineadd_1fab_1dir (fr, fabdata, flo,fhi, dir, boxno, zeroFirst, nfluxes, scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       type(c_ptr), value,intent(in) :: fabdata
       integer(c_int),intent(in) :: flo(*), fhi(*)
       integer(c_int),intent(IN),value :: dir, boxno, zeroFirst, nfluxes
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_fineadd_1fab_1dir

     subroutine amrex_fi_fluxregister_crseinit (fr, flxs, scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       type(c_ptr), intent(in) :: flxs(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_crseinit

     subroutine amrex_fi_fluxregister_crseadd (fr, flxs, scale, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr, geom
       type(c_ptr), intent(in) :: flxs(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_crseadd

     subroutine amrex_fi_fluxregister_setval (fr, val) bind(c)
       import 
       implicit none
       type(c_ptr), value :: fr
       real(amrex_real), value :: val
     end subroutine amrex_fi_fluxregister_setval

     subroutine amrex_fi_fluxregister_reflux (fr, mf, scale, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr, mf, geom
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_reflux

     subroutine amrex_fi_fluxregister_overwrite (fr, flxs, scale, geom) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr, geom
       type(c_ptr), intent(in) :: flxs(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_fluxregister_overwrite
  end interface

contains

  subroutine amrex_fluxregister_build (fr, ba, dm, ref_ratio, fine_lev, ncomp)
    type(amrex_fluxregister) :: fr
    type(amrex_boxarray), intent(in) :: ba
    type(amrex_distromap), intent(in) :: dm
    integer, intent(in) :: ref_ratio, fine_lev, ncomp
    fr%owner = .true.
    fr%flev  = fine_lev
    call amrex_fi_new_fluxregister(fr%p, ba%p, dm%p, ref_ratio, fine_lev, ncomp)
  end subroutine amrex_fluxregister_build
  
  impure elemental subroutine amrex_fluxregister_destroy (this)
    type(amrex_fluxregister), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_fluxregister(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr
  end subroutine amrex_fluxregister_destroy
    
  subroutine amrex_fluxregister_assign (dst, src)
    class(amrex_fluxregister), intent(inout) :: dst
    type (amrex_fluxregister), intent(in   ) :: src
    dst%owner = .false.
    dst%flev  = src%flev
    dst%p     = src%p
  end subroutine amrex_fluxregister_assign

  subroutine amrex_fluxregister_fineadd (this, fluxes, scale)
    class(amrex_fluxregister), intent(inout) :: this
    type(amrex_multifab), intent(in) :: fluxes(amrex_spacedim)
    real(amrex_real), intent(in) :: scale
    integer :: dim
    type(c_ptr) :: mf(amrex_spacedim)
    do dim = 1, amrex_spacedim
       mf(dim) = fluxes(dim)%p
    end do
    call amrex_fi_fluxregister_fineadd(this%p, mf, scale)
  end subroutine amrex_fluxregister_fineadd

  subroutine amrex_fluxregister_fineadd_1fab (this, fluxfabs, gridIdx, scale, zeroFirst)
    class(amrex_fluxregister), intent(inout) :: this
    type(amrex_fab),  intent(in) :: fluxfabs(amrex_spacedim)
    integer(c_int),   intent(in) :: gridIdx
    real(amrex_real), intent(in) :: scale
    logical,optional, intent(in) :: zeroFirst
    integer :: dir, nc, myZeroFirst
    type(c_ptr) :: cp
    integer,dimension(3) :: flo=(/1,1,1/),fhi=(/1,1,1/)
    real(amrex_real), pointer, dimension(:,:,:,:) :: fp
    do dir = 1, amrex_spacedim
       associate(fab => fluxfabs(dir) )
         flo(1:amrex_spacedim) = fab%bx%lo(1:amrex_spacedim)
         fhi(1:amrex_spacedim) = fab%bx%hi(1:amrex_spacedim)  
         fp => fab%dataPtr()
         cp = c_loc(fp(flo(1),flo(2),flo(3),1))
         nc =  fab%nc
       end associate
       myZeroFirst = 0
       if (present(zeroFirst)) then
          if (zeroFirst) myZeroFirst = 1
       end if
       call amrex_fi_fluxregister_fineadd_1fab_1dir(this%p, cp, flo,fhi, dir-1, gridIdx, myZeroFirst, nc, scale)
    end do
  end subroutine amrex_fluxregister_fineadd_1fab

  subroutine amrex_fluxregister_crseinit (this, fluxes, scale)
    class(amrex_fluxregister), intent(inout) :: this
    type(amrex_multifab), intent(in) :: fluxes(amrex_spacedim)
    real(amrex_real), intent(in) :: scale
    integer :: dim
    type(c_ptr) :: mf(amrex_spacedim)
    do dim = 1, amrex_spacedim
       mf(dim) = fluxes(dim)%p
    end do
    call amrex_fi_fluxregister_crseinit(this%p, mf, scale)
  end subroutine amrex_fluxregister_crseinit

  subroutine amrex_fluxregister_crseadd (this, fluxes, scale)
    use amrex_amrcore_module, only : amrex_geom
    class(amrex_fluxregister), intent(inout) :: this
    type(amrex_multifab), intent(in) :: fluxes(amrex_spacedim)
    real(amrex_real), intent(in) :: scale
    integer :: dim
    type(c_ptr) :: mf(amrex_spacedim)
    do dim = 1, amrex_spacedim
       mf(dim) = fluxes(dim)%p
    end do
    call amrex_fi_fluxregister_crseadd(this%p, mf, scale, amrex_geom(this%flev-1)%p)
  end subroutine amrex_fluxregister_crseadd

  subroutine amrex_fluxregister_setval (this, val)
    class(amrex_fluxregister), intent(inout) :: this
    real(amrex_real), intent(in) :: val
    call amrex_fi_fluxregister_setval(this%p, val)
  end subroutine amrex_fluxregister_setval

  subroutine amrex_fluxregister_reflux (this, mf, scale)
    use amrex_amrcore_module, only : amrex_geom
    class(amrex_fluxregister), intent(inout) :: this
    type(amrex_multifab), intent(in) :: mf
    real(amrex_real), intent(in) :: scale
    call amrex_fi_fluxregister_reflux(this%p, mf%p, scale, amrex_geom(this%flev-1)%p)
  end subroutine amrex_fluxregister_reflux

  subroutine amrex_fluxregister_overwrite (this, fluxes, scale)
    use amrex_amrcore_module, only : amrex_geom
    class(amrex_fluxregister), intent(inout) :: this
    type(amrex_multifab), intent(inout) :: fluxes(amrex_spacedim)
    real(amrex_real), intent(in) :: scale
    integer :: dim
    type(c_ptr) :: mf(amrex_spacedim)
    do dim = 1, amrex_spacedim
       mf(dim) = fluxes(dim)%p
    end do
    call amrex_fi_fluxregister_overwrite(this%p, mf, scale, amrex_geom(this%flev-1)%p)
  end subroutine amrex_fluxregister_overwrite

end module amrex_fluxregister_module

