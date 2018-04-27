
module amrex_paralleldescriptor_module
  use iso_c_binding
  implicit none

  interface
     integer(c_int) function amrex_fi_pd_myproc () bind(c)
       import
       implicit none
     end function amrex_fi_pd_myproc

     integer(c_int) function amrex_fi_pd_nprocs () bind(c)
       import
       implicit none
     end function amrex_fi_pd_nprocs
     
     integer(c_int) function amrex_fi_pd_ioprocessor () bind(c)
       import
       implicit none
     end function amrex_fi_pd_ioprocessor
  end interface

  private
  public :: amrex_pd_myproc, amrex_pd_nprocs, amrex_pd_ioprocessor

contains

  integer function amrex_pd_myproc ()
    amrex_pd_myproc = amrex_fi_pd_myproc()
  end function amrex_pd_myproc

  integer function amrex_pd_nprocs ()
    amrex_pd_nprocs = amrex_fi_pd_nprocs()
  end function amrex_pd_nprocs

  logical function amrex_pd_ioprocessor ()
    integer(c_int) :: i
    i = amrex_fi_pd_ioprocessor()
    amrex_pd_ioprocessor = (i.ne.0)
  end function amrex_pd_ioprocessor

end module amrex_paralleldescriptor_module
