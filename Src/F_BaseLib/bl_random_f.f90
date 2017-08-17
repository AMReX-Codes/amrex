
! This module uses C++11 random number generator, Mersenne Twister
module bl_random_module

  use iso_c_binding
  use bl_error_module
  use parallel
  use fabio_module
  use bl_types

  implicit none

  private 
  public :: bl_rng_get, &
       bl_rng_build_engine, bl_rng_destroy_engine, bl_rng_save_engine, &
       bl_rng_copy_engine, bl_rng_restore_engine, bl_rng_build_distro, &
       bl_rng_destroy_distro, bl_rng_save_distro, bl_rng_restore_distro, &
       bl_rng_engine, bl_rng_uniform_real, bl_rng_normal, bl_rng_poisson, &
       bl_rng_binomial

  type bl_rng_engine
     type(c_ptr) :: p = c_null_ptr
  end type bl_rng_engine

  type bl_rng_uniform_real
     type(c_ptr) :: p = c_null_ptr
  end type bl_rng_uniform_real

  type bl_rng_normal
     type(c_ptr) :: p = c_null_ptr
  end type bl_rng_normal

  type bl_rng_poisson
     type(c_ptr) :: p = c_null_ptr
  end type bl_rng_poisson

  type bl_rng_binomial
     type(c_ptr) :: p = c_null_ptr
  end type bl_rng_binomial

  interface bl_rng_build_distro
     module procedure bl_rng_build_uniform_real
     module procedure bl_rng_build_normal
     module procedure bl_rng_build_poisson
     module procedure bl_rng_build_binomial
  end interface bl_rng_build_distro

  interface bl_rng_destroy_distro
     module procedure bl_rng_destroy_uniform_real
     module procedure bl_rng_destroy_normal
     module procedure bl_rng_destroy_poisson
     module procedure bl_rng_destroy_binomial
  end interface bl_rng_destroy_distro

  interface bl_rng_get
     module procedure bl_rng_get_uniform_real
     module procedure bl_rng_get_normal
     module procedure bl_rng_get_poisson
     module procedure bl_rng_get_binomial
  end interface bl_rng_get

  interface bl_rng_save_distro
     module procedure bl_rng_save_uniform_real
     module procedure bl_rng_save_normal
     module procedure bl_rng_save_poisson
     module procedure bl_rng_save_binomial
  end interface bl_rng_save_distro

  interface bl_rng_restore_distro
     module procedure bl_rng_restore_uniform_real
     module procedure bl_rng_restore_normal
     module procedure bl_rng_restore_poisson
     module procedure bl_rng_restore_binomial
  end interface bl_rng_restore_distro

  interface
     function bl_rng_random_uint_c() result (r) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       integer(c_int) :: r
     end function bl_rng_random_uint_c
  end interface

  ! engine
  interface
     subroutine bl_rng_new_engine_c(eng, s, rank, nprocs) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: eng
       integer(c_int), intent(in), value :: s, rank, nprocs
     end subroutine bl_rng_new_engine_c
     
     subroutine bl_rng_delete_engine_c(eng) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: eng
     end subroutine bl_rng_delete_engine_c

     subroutine bl_rng_save_engine_c(eng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: eng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_save_engine_c

     subroutine bl_rng_copy_engine_c(eng_dst, eng_src) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: eng_dst
       type(c_ptr), value :: eng_src
     end subroutine bl_rng_copy_engine_c

     subroutine bl_rng_restore_engine_c(eng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: eng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_restore_engine_c
  end interface

  ! uniform real distribution
  interface
     subroutine bl_rng_new_uniform_real_c(rng, a, b) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       real(c_double), intent(in), value :: a, b
     end subroutine bl_rng_new_uniform_real_c

     subroutine bl_rng_delete_uniform_real_c(rng) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
     end subroutine bl_rng_delete_uniform_real_c

     function bl_rng_get_uniform_real_c(rng,eng) result(r) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng, eng
       real(c_double) :: r
     end function bl_rng_get_uniform_real_c

     subroutine bl_rng_save_uniform_real_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_save_uniform_real_c

     subroutine bl_rng_restore_uniform_real_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_restore_uniform_real_c
  end interface

  ! normal distribution
  interface
     subroutine bl_rng_new_normal_c(rng, mean, stddev) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       real(c_double), intent(in), value :: mean, stddev
     end subroutine bl_rng_new_normal_c

     subroutine bl_rng_delete_normal_c(rng) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
     end subroutine bl_rng_delete_normal_c

     function bl_rng_get_normal_c(rng, eng) result(r) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng, eng
       real(c_double) :: r
     end function bl_rng_get_normal_c

     subroutine bl_rng_save_normal_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_save_normal_c

     subroutine bl_rng_restore_normal_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_restore_normal_c
  end interface

  ! poisson distribution
  interface
     subroutine bl_rng_new_poisson_c(rng, mean) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       real(c_double), intent(in), value :: mean
     end subroutine bl_rng_new_poisson_c

     subroutine bl_rng_delete_poisson_c(rng) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
     end subroutine bl_rng_delete_poisson_c

     function bl_rng_get_poisson_c(rng, eng) result(r) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng, eng
       integer(c_int) :: r
     end function bl_rng_get_poisson_c

     subroutine bl_rng_save_poisson_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_save_poisson_c

     subroutine bl_rng_restore_poisson_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_restore_poisson_c
  end interface

  ! binomial distribution
  interface
     subroutine bl_rng_new_binomial_c(rng, t, p) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       integer(c_int), intent(in), value :: t
       real(c_double), intent(in), value :: p
     end subroutine bl_rng_new_binomial_c

     subroutine bl_rng_delete_binomial_c(rng) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
     end subroutine bl_rng_delete_binomial_c

     function bl_rng_get_binomial_c(rng, eng) result(r) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng, eng
       integer(c_int) :: r
     end function bl_rng_get_binomial_c

     subroutine bl_rng_save_binomial_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr), value :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_save_binomial_c

     subroutine bl_rng_restore_binomial_c(rng, name) bind(c)
       use, intrinsic :: iso_c_binding
       implicit none
       type(c_ptr) :: rng
       character(kind=c_char), intent(in) :: name(*)
     end subroutine bl_rng_restore_binomial_c
  end interface

contains

  subroutine bl_rng_filename(filename, dirname)
    character(kind=c_char), pointer, intent(inout) :: filename(:)
    character(len=*), intent(in) :: dirname
    integer :: i, n
    character(len=16) :: procname
    character(len=128) :: fullname
    write(procname, *) parallel_myproc()
    fullname = trim(dirname) // "/s" // trim(adjustl(procname))
    n = len_trim(fullname)
    allocate(filename(n+1))
    do i = 1, n
       filename(i) = fullname(i:i)
    end do
    filename(n+1) = c_null_char
  end subroutine bl_rng_filename

  function bl_rng_init_seed(s) result(r)
    integer(c_int), intent(in) :: s
    integer :: r
    integer, save :: rstatic
    if (s .eq. 0) then
       !$omp master
       rstatic = bl_rng_random_uint_c()
       call parallel_bcast(rstatic)
       if (parallel_IOProcessor()) then
          print*,'seed = 0 --> picking a random root seed',rstatic
       end if
       !$omp end master
       !$omp barrier
       r = rstatic
    else if (s .gt. 0) then
       r = s
    else
       call bl_error("bl_rng: seed must >= 0")
    end if
  end function bl_rng_init_seed

  !
  ! engine
  !
  subroutine bl_rng_build_engine(eng, seed)
    type(bl_rng_engine), intent(inout) :: eng
    integer(c_int), intent(in) :: seed
    integer(c_int) :: lseed
    lseed = bl_rng_init_seed(seed)
    call bl_rng_new_engine_c(eng%p, lseed, parallel_myproc(), parallel_nprocs())
  end subroutine bl_rng_build_engine
  !
  subroutine bl_rng_destroy_engine(eng)
    type(bl_rng_engine), intent(inout) :: eng
    if (c_associated(eng%p)) call bl_rng_delete_engine_c(eng%p)
    eng%p = c_null_ptr
  end subroutine bl_rng_destroy_engine
  !
  subroutine bl_rng_save_engine(eng, dirname)
    type(bl_rng_engine), intent(in) :: eng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (parallel_IOProcessor()) then
       call fabio_mkdir(dirname)
    end if
    call parallel_barrier()
    call bl_rng_filename(filename, dirname)
    call bl_rng_save_engine_c(eng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_save_engine
  !
  subroutine bl_rng_copy_engine(eng_dst, eng_src)
    type(bl_rng_engine), intent(inout) :: eng_dst
    type(bl_rng_engine), intent(in   ) :: eng_src
    call parallel_barrier()
    call bl_rng_copy_engine_c(eng_dst%p,eng_src%p)
  end subroutine bl_rng_copy_engine
  !
  subroutine bl_rng_restore_engine(eng, dirname)
    type(bl_rng_engine), intent(inout) :: eng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (c_associated(eng%p)) call bl_rng_destroy_engine(eng) 
    call bl_rng_filename(filename, dirname)
    call bl_rng_restore_engine_c(eng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_restore_engine

  ! 
  ! uniform real distribution
  !
  subroutine bl_rng_build_uniform_real(rng, a, b)
    type(bl_rng_uniform_real), intent(inout) :: rng
    real(c_double), intent(in) :: a, b
    call bl_rng_new_uniform_real_c(rng%p, a, b)
  end subroutine bl_rng_build_uniform_real
  !
  subroutine bl_rng_destroy_uniform_real(rng)
    type(bl_rng_uniform_real), intent(inout) :: rng
    if (c_associated(rng%p)) call bl_rng_delete_uniform_real_c(rng%p)
    rng%p = c_null_ptr
  end subroutine bl_rng_destroy_uniform_real
  !
  function bl_rng_get_uniform_real(rng, eng) result(r)
    type(bl_rng_uniform_real), intent(inout) :: rng
    type(bl_rng_engine), intent(inout) :: eng
    real(c_double) :: r
    r = bl_rng_get_uniform_real_c(rng%p, eng%p)    
  end function bl_rng_get_uniform_real
  !
  subroutine bl_rng_save_uniform_real(rng, dirname)
    type(bl_rng_uniform_real), intent(in) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (parallel_IOProcessor()) then
       call fabio_mkdir(dirname)
    end if
    call parallel_barrier()
    call bl_rng_filename(filename, dirname)
    call bl_rng_save_uniform_real_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_save_uniform_real
  !
  subroutine bl_rng_restore_uniform_real(rng, dirname)
    type(bl_rng_uniform_real), intent(inout) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (c_associated(rng%p)) call bl_rng_destroy_distro(rng) 
    call bl_rng_filename(filename, dirname)
    call bl_rng_restore_uniform_real_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_restore_uniform_real

  ! 
  ! normal distribution
  !
  subroutine bl_rng_build_normal(rng,  mean, stddev)
    type(bl_rng_normal), intent(inout) :: rng
    real(c_double), intent(in) :: mean, stddev
    call bl_rng_new_normal_c(rng%p, mean, stddev)
  end subroutine bl_rng_build_normal
  !
  subroutine bl_rng_destroy_normal(rng)
    type(bl_rng_normal), intent(inout) :: rng
    if (c_associated(rng%p)) call bl_rng_delete_normal_c(rng%p)
    rng%p = c_null_ptr
  end subroutine bl_rng_destroy_normal
  !
  function bl_rng_get_normal(rng, eng) result(r)
    type(bl_rng_normal), intent(inout) :: rng
    type(bl_rng_engine), intent(inout) :: eng
    real(c_double) :: r
    r = bl_rng_get_normal_c(rng%p, eng%p)    
  end function bl_rng_get_normal
  !
  subroutine bl_rng_save_normal(rng, dirname)
    type(bl_rng_normal), intent(in) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (parallel_IOProcessor()) then
       call fabio_mkdir(dirname)
    end if
    call parallel_barrier()
    call bl_rng_filename(filename, dirname)
    call bl_rng_save_normal_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_save_normal
  !
  subroutine bl_rng_restore_normal(rng, dirname)
    type(bl_rng_normal), intent(inout) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (c_associated(rng%p)) call bl_rng_destroy_distro(rng) 
    call bl_rng_filename(filename, dirname)
    call bl_rng_restore_normal_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_restore_normal

  ! 
  ! poisson distribution
  !
  subroutine bl_rng_build_poisson(rng, mean)
    type(bl_rng_poisson), intent(inout) :: rng
    real(c_double), intent(in) :: mean
    call bl_rng_new_poisson_c(rng%p, mean)
  end subroutine bl_rng_build_poisson
  !
  subroutine bl_rng_destroy_poisson(rng)
    type(bl_rng_poisson), intent(inout) :: rng
    if (c_associated(rng%p)) call bl_rng_delete_poisson_c(rng%p)
    rng%p = c_null_ptr
  end subroutine bl_rng_destroy_poisson
  !
  function bl_rng_get_poisson(rng, eng) result(r)
    type(bl_rng_poisson), intent(inout) :: rng
    type(bl_rng_engine), intent(inout) :: eng
    integer(c_int) :: r
    r = bl_rng_get_poisson_c(rng%p, eng%p)    
  end function bl_rng_get_poisson
  !
  subroutine bl_rng_save_poisson(rng, dirname)
    type(bl_rng_poisson), intent(in) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (parallel_IOProcessor()) then
       call fabio_mkdir(dirname)
    end if
    call parallel_barrier()
    call bl_rng_filename(filename, dirname)
    call bl_rng_save_poisson_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_save_poisson
  !
  subroutine bl_rng_restore_poisson(rng, dirname)
    type(bl_rng_poisson), intent(inout) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (c_associated(rng%p)) call bl_rng_destroy_distro(rng) 
    call bl_rng_filename(filename, dirname)
    call bl_rng_restore_poisson_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_restore_poisson

  ! 
  ! binomial distribution
  !
  subroutine bl_rng_build_binomial(rng, t, p)
    type(bl_rng_binomial), intent(inout) :: rng
    integer(c_int), intent(in) :: t
    real(c_double), intent(in) :: p
    call bl_rng_new_binomial_c(rng%p, t, p)
  end subroutine bl_rng_build_binomial
  !
  subroutine bl_rng_destroy_binomial(rng)
    type(bl_rng_binomial), intent(inout) :: rng
    if (c_associated(rng%p)) call bl_rng_delete_binomial_c(rng%p)
    rng%p = c_null_ptr
  end subroutine bl_rng_destroy_binomial
  !
  function bl_rng_get_binomial(rng, eng) result(r)
    type(bl_rng_binomial), intent(inout) :: rng
    type(bl_rng_engine), intent(inout) :: eng
    integer(c_int) :: r
    r = bl_rng_get_binomial_c(rng%p, eng%p)    
  end function bl_rng_get_binomial
  !
  subroutine bl_rng_save_binomial(rng, dirname)
    type(bl_rng_binomial), intent(in) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (parallel_IOProcessor()) then
       call fabio_mkdir(dirname)
    end if
    call parallel_barrier()
    call bl_rng_filename(filename, dirname)
    call bl_rng_save_binomial_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_save_binomial
  !
  subroutine bl_rng_restore_binomial(rng, dirname)
    type(bl_rng_binomial), intent(inout) :: rng
    character(len=*), intent(in) :: dirname
    character(kind=c_char), pointer :: filename(:)
    if (c_associated(rng%p)) call bl_rng_destroy_distro(rng) 
    call bl_rng_filename(filename, dirname)
    call bl_rng_restore_binomial_c(rng%p,filename)
    deallocate(filename)
  end subroutine bl_rng_restore_binomial

end module bl_random_module

