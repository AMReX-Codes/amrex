module amrex_fort_module

  use iso_c_binding, only : c_char, c_short, c_int, c_long, c_float, c_double, c_size_t, c_ptr

  implicit none

  integer, parameter ::    bl_spacedim = AMREX_SPACEDIM
  integer, parameter :: amrex_spacedim = AMREX_SPACEDIM

#ifdef BL_USE_FLOAT
  integer, parameter :: amrex_real = c_float
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 4_c_size_t
#else
  integer, parameter :: amrex_real = c_double
  ! We could/should use Fortran 2008 c_sizeof here.
  integer (kind=c_size_t), parameter :: amrex_real_size = 8_c_size_t
#endif

#ifdef BL_SINGLE_PRECISION_PARTICLES
  integer, parameter :: amrex_particle_real = c_float
#else
  integer, parameter :: amrex_particle_real = c_double
#endif

  interface
     function amrex_malloc (s) bind(c,name='amrex_malloc')
       import
       integer(c_size_t), intent(in), value :: s
       type(c_ptr) :: amrex_malloc
     end function amrex_malloc

     subroutine amrex_free (p) bind(c,name='amrex_free')
       import
       type(c_ptr), value :: p
     end subroutine amrex_free

     function amrex_random () bind(c,name='amrex_random')
       import
       real(c_double) :: amrex_random
     end function amrex_random

     function amrex_random_int (n) bind(c,name='amrex_random_int')
       import
       integer(c_long), intent(in), value :: n
       integer(c_long) :: amrex_random_int
     end function amrex_random_int
  end interface

contains

  function amrex_coarsen_intvect (n, iv, rr) result(civ)
    integer, intent(in) :: n, rr
    integer, intent(in) :: iv(n)
    integer :: civ(n)
    integer :: i
    do i = 1, n
       if (iv(i) .lt. 0) then
          civ(i) = -abs(iv(i)+1)/rr - 1
       else
          civ(i) = iv(i)/rr
       end if
    end do
  end function amrex_coarsen_intvect


  subroutine amrex_get_loop_bounds(blo, bhi, lo, hi) bind(c, name='amrex_get_loop_bounds')

    implicit none

    integer, intent(in   ) :: lo(3), hi(3)
    integer, intent(inout) :: blo(3), bhi(3)

    !$gpu

#ifdef AMREX_USE_GPU_PRAGMA
    ! Get our spatial index based on the CUDA thread index

    blo(1) = lo(1) + (threadIdx%x - 1) + blockDim%x * (blockIdx%x - 1)
    blo(2) = lo(2) + (threadIdx%y - 1) + blockDim%y * (blockIdx%y - 1)
    blo(3) = lo(3) + (threadIdx%z - 1) + blockDim%z * (blockIdx%z - 1)

    ! If we have more threads than zones, set hi < lo so that the
    ! loop iteration gets skipped.

    if (blo(1) .gt. hi(1) .or. blo(2) .gt. hi(2) .or. blo(3) .gt. hi(3)) then
       bhi = blo - 1
    else
       bhi = blo
    endif
#else
    blo = lo
    bhi = hi
#endif

  end subroutine amrex_get_loop_bounds

  subroutine amrex_add(x, y)

    implicit none

    ! Add y to x.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = x + y

  end subroutine amrex_add



#ifdef AMREX_USE_CUDA_FORTRAN
  ! Note that the device versions of these
  ! functions are intentionally constructed
  ! by hand rather than scripted.

  attributes(device) subroutine amrex_add_device(x, y)

    implicit none

    ! Add y to x atomically on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

    t = atomicAdd(x, y)

  end subroutine amrex_add_device
#endif



#ifdef AMREX_USE_CUDA_FORTRAN
  attributes(device) function warpReduceSum(x) result(y)
    ! Reduce within a warp.
    ! https://devblogs.nvidia.com/faster-parallel-reductions-kepler/

    implicit none

    real(amrex_real), intent(in) :: x

    real(amrex_real) :: y

    integer :: offset

    offset = warpsize / 2

    y = x

    do while (offset > 0)

       y = y + __shfl_down(y, offset)

       offset = offset / 2

    end do

  end function warpReduceSum



  attributes(device) function blockReduceSum(x) result(y)
    ! Reduce within a threadblock.
    ! https://devblogs.nvidia.com/faster-parallel-reductions-kepler/

    implicit none

    real(amrex_real), intent(in) :: x

    real(amrex_real) :: y

    real(amrex_real), shared :: s(0:(AMREX_CUDA_MAX_THREADS/warpsize) - 1)

    integer :: lane, wid

    lane = mod(threadIdx%x-1, warpsize)
    wid = (threadIdx%x-1) / warpsize

    y = x
    y = warpReduceSum(y)

    if (lane == 0) then
       s(wid) = y
    end if

    call syncthreads()

    if ((threadIdx%x-1) < blockDim%x / warpsize) then
       y = s(lane)
    else
       y = 0
    end if

    if (wid == 0) then
       y = warpReduceSum(y)
    end if

  end function blockReduceSum
#endif



  subroutine amrex_reduce_add(x, y)

    implicit none

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = x + y

  end subroutine amrex_reduce_add



#ifdef AMREX_USE_CUDA_FORTRAN
  attributes(device) subroutine amrex_reduce_add_device(x, y)
    ! Do a shared memory reduction within a threadblock,
    ! then do an atomic add with a single thread in the block.

    implicit none

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

    t = y

    t = blockReduceSum(t)

    if (threadIdx%x == 1) then

       t = atomicAdd(x, t)

    end if

  end subroutine amrex_reduce_add_device
#endif



  subroutine amrex_subtract(x, y)

    implicit none

    ! Subtract y from x.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = x - y

  end subroutine amrex_subtract



#ifdef AMREX_USE_CUDA_FORTRAN
  attributes(device) subroutine amrex_subtract_device(x, y)

    implicit none

    ! Subtract y from x atomically on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

    t = atomicSub(x, y)

  end subroutine amrex_subtract_device
#endif



  subroutine amrex_max(x, y)

    implicit none

    ! Set in x the maximum of x and y.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = max(x, y)

  end subroutine amrex_max



#ifdef AMREX_USE_CUDA_FORTRAN
  attributes(device) subroutine amrex_max_device(x, y)

    implicit none

    ! Set in x the maximum of x and y atomically on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

    t = atomicMax(x, y)

  end subroutine amrex_max_device
#endif



  subroutine amrex_min(x, y)

    implicit none

    ! Set in x the minimum of x and y.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    x = min(x, y)

  end subroutine amrex_min



#ifdef AMREX_USE_CUDA_FORTRAN
  attributes(device) subroutine amrex_min_device(x, y)

    implicit none

    ! Set in x the minimum of x and y atomically on the GPU.

    real(amrex_real), intent(in   ) :: y
    real(amrex_real), intent(inout) :: x

    real(amrex_real) :: t

    t = atomicMin(x, y)

  end subroutine amrex_min_device
#endif

end module amrex_fort_module
