
module amrex_flash_fluxregister_module

  use iso_c_binding
  use amrex_base_module
  use amrex_amrcore_module

  implicit none

  private

  public :: amrex_flash_fluxregister_build, amrex_flash_fluxregister_destroy

  type, public :: amrex_flash_fluxregister
     logical     :: owner = .false.
     type(c_ptr) :: p     = c_null_ptr
     integer     :: flev  = -1
   contains
     generic   :: assignment(=) => amrex_flash_fluxregister_assign ! shallow copy
     procedure :: load          => amrex_flash_fluxregister_load
     procedure :: store         => amrex_flash_fluxregister_store
     procedure :: communicate   => amrex_flash_fluxregister_communicate
     procedure, private :: amrex_flash_fluxregister_assign
     procedure, private :: amrex_flash_fluxregister_load
     procedure, private :: amrex_flash_fluxregister_store
#if !defined(__GFORTRAN__) || (__GNUC__ > 4)
     final :: amrex_flash_fluxregister_destroy
#endif
  end type amrex_flash_fluxregister

  interface
     subroutine amrex_fi_new_flash_fluxregister (fr,fba,cba,fdm,cdm,fgm,cgm,rr,nc) bind(c)
       import
       implicit none
       type(c_ptr) :: fr
       type(c_ptr), value :: fba,cba,fdm,cdm,fgm,cgm
       integer(c_int), value :: rr, nc
     end subroutine amrex_fi_new_flash_fluxregister

     subroutine amrex_fi_delete_flash_fluxregister (fr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
     end subroutine amrex_fi_delete_flash_fluxregister

     subroutine amrex_fi_flash_fluxregister_load (fr,cgid,dir,fabdata,flo,fhi,nc,scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       real(amrex_real), intent(inout) :: fabdata(*)
       integer(c_int), value, intent(in) :: cgid, dir, nc
       integer(c_int), intent(in) :: flo(*), fhi(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_flash_fluxregister_load

     subroutine amrex_fi_flash_fluxregister_store (fr,cgid,dir,fabdata,flo,fhi,nc,scale) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
       real(amrex_real), intent(in) :: fabdata(*)
       integer(c_int), value, intent(in) :: cgid, dir, nc
       integer(c_int), intent(in) :: flo(*), fhi(*)
       real(amrex_real), value :: scale
     end subroutine amrex_fi_flash_fluxregister_store

     subroutine amrex_fi_flash_fluxregister_communicate (fr) bind(c)
       import
       implicit none
       type(c_ptr), value :: fr
     end subroutine amrex_fi_flash_fluxregister_communicate
  end interface

contains

  subroutine amrex_flash_fluxregister_build (fr, fba, cba, fdm, cdm, fine_lev, ncomp)
    type(amrex_flash_fluxregister) :: fr
    integer, intent(in) :: fine_lev, ncomp
    type(amrex_boxarray), intent(in) :: fba, cba
    type(amrex_distromap), intent(in) :: fdm, cdm
    !
    type(amrex_geometry) :: fgm, cgm
    fr%owner = .true.
    fr%flev  = fine_lev
    fgm = amrex_get_geometry(fine_lev)
    cgm = amrex_get_geometry(fine_lev-1)
    call amrex_fi_new_flash_fluxregister(fr%p,fba%p,cba%p,fdm%p,cdm%p,fgm%p,cgm%p, &
         amrex_ref_ratio(fine_lev-1), ncomp)
  end subroutine amrex_flash_fluxregister_build

  subroutine amrex_flash_fluxregister_destroy (this)
    type(amrex_flash_fluxregister), intent(inout) :: this
    if (this%owner) then
       if (c_associated(this%p)) then
          call amrex_fi_delete_flash_fluxregister(this%p)
       end if
    end if
    this%owner = .false.
    this%p = c_null_ptr;
  end subroutine amrex_flash_fluxregister_destroy

  subroutine amrex_flash_fluxregister_assign (dst, src)
    class(amrex_flash_fluxregister), intent(inout) :: dst
    type (amrex_flash_fluxregister), intent(in   ) :: src
    dst%owner = .false.
    dst%flev  = src%flev
    dst%p     = src%p
  end subroutine amrex_flash_fluxregister_assign

  subroutine amrex_flash_fluxregister_communicate (this)
    class(amrex_flash_fluxregister), intent(inout) :: this
    call amrex_fi_flash_fluxregister_communicate(this%p)
  end subroutine amrex_flash_fluxregister_communicate

  subroutine amrex_flash_fluxregister_store (this, flux, flo, fhi, grid_idx, dir, scale)
    class(amrex_flash_fluxregister), intent(inout) :: this
    integer, intent(in) :: flo(*), fhi(*), grid_idx, dir
    real(amrex_real), intent(in) :: flux(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),flo(4):fhi(4))
    real(amrex_real), optional, intent(in) :: scale
    !
    real(amrex_real) :: my_scale
    if (present(scale)) then
       my_scale = scale
    else
       my_scale = 1._amrex_real
    end if
    call amrex_fi_flash_fluxregister_store(this%p, grid_idx, dir, &
         flux, flo, fhi, fhi(4)-flo(4)+1, my_scale)
  end subroutine amrex_flash_fluxregister_store

  subroutine amrex_flash_fluxregister_load (this, flux, flo, fhi, grid_idx, dir, scale)
    class(amrex_flash_fluxregister), intent(inout) :: this
    integer, intent(in) :: flo(*), fhi(*), grid_idx, dir
    real(amrex_real), intent(inout) :: flux(flo(1):fhi(1),flo(2):fhi(2),flo(3):fhi(3),flo(4):fhi(4))
    real(amrex_real), optional, intent(in) :: scale
    !
    real(amrex_real) :: my_scale
    if (present(scale)) then
       my_scale = scale
    else
       my_scale = 1._amrex_real
    end if
    call amrex_fi_flash_fluxregister_load(this%p, grid_idx, dir, &
         flux, flo, fhi, fhi(4)-flo(4)+1, my_scale)
  end subroutine amrex_flash_fluxregister_load

end module amrex_flash_fluxregister_module
