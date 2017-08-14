
module amrex_plotfile_module

  use iso_c_binding
  use amrex_fort_module
  use amrex_multifab_module
  use amrex_geometry_module
  use amrex_string_module

  implicit none

  private

  public :: amrex_write_plotfile

  interface
     subroutine amrex_fi_write_plotfile (name, nlevs, mf, varname, geom, t, steps, rr) bind(c)
       import
       implicit none
       character(kind=c_char), intent(in) :: name(*)
       integer(c_int), value :: nlevs
       type(c_ptr), intent(in) :: mf(*), varname(*), geom(*)
       real(amrex_real), value :: t
       integer(c_int), intent(in) :: steps(*), rr(*)
     end subroutine amrex_fi_write_plotfile
  end interface

contains

  subroutine amrex_write_plotfile (name, nlevs, mf, varname, geom, t, steps, rr)
    character(len=*), intent(in) :: name
    integer, intent(in) :: nlevs, steps(*), rr(*)
    type(amrex_multifab), intent(in) :: mf(*)
    type(amrex_string), intent(in), target :: varname(*)
    type(amrex_geometry), intent(in) :: geom(*)
    real(amrex_real), intent(in) :: t

    integer :: lev, icomp
    type(c_ptr) :: mf_c(nlevs), geom_c(nlevs)
    type(c_ptr) :: varname_c(mf(1)%ncomp())

    do lev = 1, nlevs
       mf_c(lev) = mf(lev)%p
       geom_c(lev) = geom(lev)%p
    end do

    do icomp = 1, size(varname_c)
       varname_c(icomp) = c_loc(varname(icomp)%data(1))
    end do

    call amrex_fi_write_plotfile(amrex_string_f_to_c(name), nlevs, mf_c, varname_c, geom_c, &
         t, steps, rr)
  end subroutine amrex_write_plotfile

end module amrex_plotfile_module
