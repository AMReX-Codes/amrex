
module amrex_paralleldescriptor_module
  use iso_c_binding
  use amrex_fort_module
  implicit none

  interface
     function amrex_fi_pd_myproc () bind(c)
       import
       implicit none
       integer(c_int) :: amrex_fi_pd_myproc
     end function amrex_fi_pd_myproc

     function amrex_fi_pd_nprocs () bind(c)
       import
       implicit none
       integer(c_int) :: amrex_fi_pd_nprocs
     end function amrex_fi_pd_nprocs
     
     function amrex_fi_pd_ioprocessor () bind(c)
       import
       implicit none
       integer(c_int) :: amrex_fi_pd_ioprocessor
     end function amrex_fi_pd_ioprocessor

     function amrex_fi_pd_ioprocessor_number () bind(c)
       import
       implicit none
       integer(c_int) amrex_fi_pd_ioprocessor_number
     end function amrex_fi_pd_ioprocessor_number

     subroutine amrex_fi_pd_bcast_r (x, n, root) bind(c)
       import
       implicit none
       real(amrex_real) :: x(*)
       integer(c_int), intent(in), value :: n, root
     end subroutine amrex_fi_pd_bcast_r

     function amrex_fi_pd_wtime () bind(c)
       import
       implicit none
       real(amrex_real) :: amrex_fi_pd_wtime
     end function amrex_fi_pd_wtime
  end interface

  interface amrex_pd_bcast
     module procedure amrex_pd_bcast_r
     module procedure amrex_pd_bcast_rv
     module procedure amrex_pd_bcast_r2v
     module procedure amrex_pd_bcast_r3v
  end interface amrex_pd_bcast

  private
  public :: amrex_pd_myproc, amrex_pd_nprocs, amrex_pd_ioprocessor, amrex_pd_ioprocessor_number, &
       amrex_pd_bcast, amrex_pd_wtime

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

  integer function amrex_pd_ioprocessor_number ()
    amrex_pd_ioprocessor_number = amrex_fi_pd_ioprocessor_number()
  end function amrex_pd_ioprocessor_number

  subroutine amrex_pd_bcast_r (x, a_root)
    real(amrex_real), target :: x
    integer, intent(in), optional :: a_root
    integer :: root
    real(amrex_real) :: r(1)
    if (present(a_root)) then
       root = a_root
    else
       root = amrex_pd_ioprocessor_number()
    end if
    if (root .eq. amrex_pd_myproc()) then
       r(1) = x
    end if
    call amrex_fi_pd_bcast_r(r, 1, root)
    if (root .ne. amrex_pd_myproc()) then
       x = r(1)
    end if
  end subroutine amrex_pd_bcast_r

  subroutine amrex_pd_bcast_rv (x, a_root)
    real(amrex_real) :: x(:)
    integer, intent(in), optional :: a_root
    integer :: root
    if (present(a_root)) then
       root = a_root
    else
       root = amrex_pd_ioprocessor_number()
    end if
    call amrex_fi_pd_bcast_r(x, size(x), root)
  end subroutine amrex_pd_bcast_rv

  subroutine amrex_pd_bcast_r2v (x, a_root)
    real(amrex_real) :: x(:,:)
    integer, intent(in), optional :: a_root
    integer :: root
    if (present(a_root)) then
       root = a_root
    else
       root = amrex_pd_ioprocessor_number()
    end if
    call amrex_fi_pd_bcast_r(x, size(x), root)
  end subroutine amrex_pd_bcast_r2v

  subroutine amrex_pd_bcast_r3v (x, a_root)
    real(amrex_real) :: x(:,:,:)
    integer, intent(in), optional :: a_root
    integer :: root
    if (present(a_root)) then
       root = a_root
    else
       root = amrex_pd_ioprocessor_number()
    end if
    call amrex_fi_pd_bcast_r(x, size(x), root)
  end subroutine amrex_pd_bcast_r3v

  function amrex_pd_wtime () result(r)
    real(amrex_real) :: r
    r = amrex_fi_pd_wtime()
  end function amrex_pd_wtime

end module amrex_paralleldescriptor_module
