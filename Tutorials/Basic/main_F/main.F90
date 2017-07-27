!
! This example demonstrates:
!
! (1) By default, AMReX initializes MPI and uses MPI_COMM_WORLD as its communicator.
!     However, applications could choose to initialize MPI themselves and pass in an
!     existing communicator.
!
! (2) By default, AMReX treats command line arguments as inputs parameters.  The expected
!     format of argv is
!         executable inputs_file parm=value
!     Here, `excutable` is the filename of the executable, `inputs_file` is the file containing
!     runtime parameters used to build AMReX ParmParse database, and `parm=value` is an input
!     parameter that will override its value in `inputs_file`.  Both `inputs_file` and
!     `parm=value` are optional.  At most one `inputs_file` is allowed. Howeer, there can be
!     multiple `parm=value`s.
!
!     The parsing of the command line arguments is performed in amrex_init.  Applications can
!     choose to skip command line parsing.  Applications can also provide a procedure that
!     adds parameters to AMReX ParmParse database.  The procedure must have the signature of
!     subroutine () bind(c) (i.e., void arguments and c binding).
!

program main

  use mpi
  use amrex_base_module

  implicit none

  interface
     subroutine add_parameters () bind(c)
     end subroutine add_parameters
  end interface

  integer :: ierr

  call mpi_init(ierr)

  ! We pass MPI_COMM_WORLD.  Its duplicate will be used by AMReX.
  !
  ! arg_parmparse=.false. : No command line arguments are used to build ParmParse database.
  !
  ! proc_parmparse : We pass a procedure that adds parameters to ParmParse database..
  !
  ! All three are optional arguments.

  call amrex_init(comm=MPI_COMM_WORLD, arg_parmparse=.false., proc_parmparse=add_parameters)

  ! Testing
  call test_parameters()

  ! ...
  
  call amrex_finalize()

  call mpi_finalize(ierr)  ! We have to call this because we called MPI_Init.

end program main


subroutine add_parameters () bind(c)
  use amrex_base_module
  implicit none

  type(amrex_parmparse) :: pp

  ! prefix "amrex"
  call amrex_parmparse_build(pp,"amrex")
  call pp%add("fpe_trap_invalid", 1)   ! turn on NaN trapping, which is off by default.
  call amrex_parmparse_destroy(pp)

  ! anonymous prefix
  call amrex_parmparse_build(pp)
  call pp%add("an_int_scalar", 2)      ! integer scalar: an_int_scalar
  call pp%add("a_bool_scalar", .true.) ! logical scalar: a_bool_scalar
  call pp%addarr("a_real_array", [1._amrex_real, 2._amrex_real, 3._amrex_real]) ! real array: a_real_array
  call amrex_parmparse_destroy(pp)

  ! prefix "a_prefix"
  call amrex_parmparse_build(pp, "a_prefix")
  call pp%addarr("an_int_array", [2, 3, 4])      ! integer array: a_prefix.an_int_array
  call pp%add("a_real_scalar", 3.14_amrex_real)  ! real scalar  : a_prefix.a_real_scalar
  call pp%add("a_string", "vonNeumann")          ! character    : a_prefix.a_string
  call amrex_parmparse_destroy(pp)
end subroutine add_parameters


subroutine test_parameters ()
  use amrex_base_module
  implicit none
  
  type(amrex_parmparse) :: pp
  integer :: i
  integer, allocatable :: ia(:)
  logical :: b
  real(amrex_real) :: r
  real(amrex_real), allocatable :: ra(:)
  character(len=:), allocatable :: s

  ! anonymous prefix
  call amrex_parmparse_build(pp)
  call pp%get("an_int_scalar", i)
  call pp%get("a_bool_scalar",b)
  call pp%getarr("a_real_array", ra)
  call amrex_parmparse_destroy(pp)

  ! prefix "a_prefix"
  call amrex_parmparse_build(pp, "a_prefix")
  call pp%getarr("an_int_array", ia)
  call pp%get("a_real_scalar", r)
  call pp%get("a_string", s)
  call amrex_parmparse_destroy(pp)

  if (amrex_parallel_ioprocessor()) then
     print *, "an_int_scalar = ", i
     print *, "a_bool_scalar = ", b
     print *, "a_real_array = ", ra
     print *, "a_prefix.an_int_array = ", ia
     print *, "a_prefix.a_real_scalar = ", r
     print *, "a_prefix.a_string = ", s
  end if
end subroutine test_parameters
