module cuda_module

  use cudafor, only: cuda_stream_kind, cudaSuccess, cudaGetErrorString, cudaDeviceProp

  implicit none

  public

  integer, parameter :: max_cuda_streams = 100
  integer(kind=cuda_stream_kind) :: cuda_streams(0:max_cuda_streams) ! Note: zero will be the default stream.

  integer, private :: cuda_device_id

  type(cudaDeviceProp), private :: prop

  ! Some C++ static class members will initialize at program load,
  ! before we have the ability to initialize this module. Guard against
  ! this by only testing on the device properties if we have actually
  ! loaded them.

  logical, private :: have_prop = .false.

  integer :: stream_index

  !$omp threadprivate(stream_index)

contains

  subroutine initialize_cuda() bind(c, name='initialize_cuda')

    use cudafor, only: cudaStreamCreate, cudaGetDeviceProperties, cudaforSetDefaultStream

    implicit none

    integer :: i, cudaResult, ilen

    cudaResult = cudaStreamCreate(cuda_streams(0))
    call gpu_error_test(cudaResult)

    ! Setting stream 0 as the default stream is currently causing
    ! bugs in the PGI compiler. Disable it for now.
!    cudaResult = cudaforSetDefaultStream(cuda_streams(0))
!    call gpu_error_test(cudaResult)

    stream_index = -1

    do i = 1, max_cuda_streams
       cudaResult = cudaStreamCreate(cuda_streams(i))
       call gpu_error_test(cudaResult)
    enddo

    cuda_device_id = 0

    cudaResult = cudaGetDeviceProperties(prop, cuda_device_id)
    call gpu_error_test(cudaResult)

    have_prop = .true.

    if (prop%major < 3) then
       call bl_abort("CUDA functionality unsupported on GPUs with compute capability earlier than 3.0")
    end if

    ilen = verify(prop%name, ' ', .true.)

    print *, ""
    print *, "Using GPU: ", prop%name(1:ilen)
    print *, ""

  end subroutine initialize_cuda



  subroutine finalize_cuda() bind(c, name='finalize_cuda')

    use cudafor, only: cudaStreamDestroy, cudaProfilerStop, cudaDeviceReset

    implicit none

    integer :: i, cudaResult

    do i = 1, max_cuda_streams
       cudaResult = cudaStreamDestroy(cuda_streams(i))
       call gpu_error_test(cudaResult, abort=.false.)
    end do

    call cudaProfilerStop()

    cudaResult = cudaDeviceReset()
    call gpu_error_test(cudaResult, abort=.false.)

  end subroutine finalize_cuda



  subroutine set_gpu_stack_limit(size) bind(c, name='set_gpu_stack_limit')

    use cudafor, only: cudaDeviceSetLimit, cudaLimitStackSize, cuda_count_kind

    implicit none

    integer(cuda_count_kind), intent(in) :: size

    integer :: cudaResult

    ! Set the GPU stack size, in bytes.

    cudaResult = cudaDeviceSetLimit(cudaLimitStackSize, size)
    call gpu_error_test(cudaResult)

  end subroutine set_gpu_stack_limit



  subroutine get_cuda_device_id(id) bind(c, name='get_cuda_device_id')

    implicit none

    integer :: id

    id = cuda_device_id

  end subroutine get_cuda_device_id



  subroutine set_stream_idx(index_in) bind(c, name='set_stream_idx')

    implicit none

    integer, intent(in), value :: index_in

    stream_index = index_in

  end subroutine set_stream_idx



  subroutine get_stream_idx(index) bind(c, name='get_stream_idx')

    implicit none

    integer, intent(out) :: index

    index = stream_index

  end subroutine get_stream_idx



  integer function stream_from_index(idx)

    implicit none

    integer :: idx

    stream_from_index = MOD(idx, max_cuda_streams) + 1

  end function stream_from_index



#ifdef BL_SPACEDIM
  subroutine threads_and_blocks(lo, hi, numBlocks, numThreads)

    use cudafor, only: dim3

    implicit none

    integer, intent(in)       :: lo(BL_SPACEDIM), hi(BL_SPACEDIM)
    type(dim3), intent(inout) :: numBlocks, numThreads

    integer :: tile_size(BL_SPACEDIM)

    tile_size = hi - lo + 1

    if (BL_SPACEDIM .eq. 1) then

       numThreads % x = 256
       numThreads % y = 1
       numThreads % z = 1

       numBlocks % x = (tile_size(1) + numThreads % x - 1) / numThreads % x
       numBlocks % y = 1
       numBlocks % z = 1

    else if (BL_SPACEDIM .eq. 2) then

       numThreads % x = 16
       numThreads % y = 16
       numThreads % z = 1

       numBlocks % x = (tile_size(1) + numThreads % x - 1) / numThreads % x
       numBlocks % y = (tile_size(2) + numThreads % y - 1) / numThreads % y
       numBlocks % z = 1

    else

       numThreads % x = 8
       numThreads % y = 8
       numThreads % z = 8

       numBlocks % x = (tile_size(1) + numThreads % x - 1) / numThreads % x
       numBlocks % y = (tile_size(2) + numThreads % y - 1) / numThreads % y
       numBlocks % z = (tile_size(3) + numThreads % z - 1) / numThreads % z

    endif

  end subroutine threads_and_blocks
#endif



  subroutine gpu_malloc(x, sz) bind(c, name='gpu_malloc')

    use cudafor, only: cudaMalloc, c_devptr
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: x
    integer(c_size_t) :: sz

    integer :: cudaResult

    cudaResult = cudaMalloc(x, sz)
    call gpu_error_test(cudaResult)

  end subroutine gpu_malloc



  subroutine gpu_hostalloc(x, sz) bind(c, name='gpu_hostalloc')

    use cudafor, only: cudaHostAlloc, cudaHostAllocMapped, cudaHostAllocWriteCombined
    use iso_c_binding, only: c_size_t, c_ptr

    implicit none

    type(c_ptr) :: x
    integer :: sz

    integer :: cudaResult

    cudaResult = cudaHostAlloc(x, sz, ior(cudaHostAllocMapped, cudaHostAllocWriteCombined))
    call gpu_error_test(cudaResult)

  end subroutine gpu_hostalloc



  subroutine gpu_malloc_managed(x, sz) bind(c, name='gpu_malloc_managed')

    use cudafor, only: cudaMallocManaged, cudaMemAttachGlobal, c_devptr
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: x
    integer(c_size_t) :: sz

    integer :: cudaResult

    if ((.not. have_prop) .or. (have_prop .and. prop%managedMemory == 1)) then

       cudaResult = cudaMallocManaged(x, sz, cudaMemAttachGlobal)
       call gpu_error_test(cudaResult)

    else

       call bl_abort("The GPU does not support managed memory allocations")

    end if

  end subroutine gpu_malloc_managed



  subroutine gpu_free(x) bind(c, name='gpu_free')

    use cudafor, only: cudaFree, c_devptr, cudaErrorCudartUnloading

    implicit none

    type(c_devptr), value :: x

    integer :: cudaResult

    cudaResult = cudaFree(x)
    call gpu_error_test(cudaResult, abort=.false.)

  end subroutine gpu_free



  subroutine gpu_freehost(x) bind(c, name='gpu_freehost')

    use cudafor, only: cudaFreeHost, cudaErrorCudartUnloading
    use iso_c_binding, only: c_ptr

    implicit none

    type(c_ptr), value :: x

    integer :: cudaResult

    cudaResult = cudaFreeHost(x)
    call gpu_error_test(cudaResult, abort=.false.)

  end subroutine gpu_freehost



  subroutine gpu_host_device_ptr(x, y) bind(c, name='gpu_host_device_ptr')

    use cudafor, only: cudaHostGetDevicePointer, c_devptr
    use iso_c_binding, only: c_ptr

    type(c_devptr) :: x
    type(c_ptr), value :: y
    integer :: flags

    integer :: cudaResult

    flags = 0 ! This argument does nothing at present, but is part of the API.

    cudaResult = cudaHostGetDevicePointer(x, y, flags)
    call gpu_error_test(cudaResult)

  end subroutine gpu_host_device_ptr



  subroutine gpu_htod_memcpy_async(p_d, p_h, sz, idx) bind(c, name='gpu_htod_memcpy_async')

    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, c_devptr, cuda_stream_kind
    use iso_c_binding, only: c_ptr, c_size_t

    implicit none

    type(c_devptr), value :: p_d
    type(c_ptr), value :: p_h
    integer(c_size_t) :: sz
    integer :: idx

    integer :: cudaResult

    cudaResult = cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, cuda_streams(stream_from_index(idx)))
    call gpu_error_test(cudaResult)

  end subroutine gpu_htod_memcpy_async



  subroutine gpu_dtoh_memcpy_async(p_h, p_d, sz, idx) bind(c, name='gpu_dtoh_memcpy_async')

    use cudafor, only: cudaMemcpyAsync, cudaMemcpyDeviceToHost, c_devptr
    use iso_c_binding, only: c_ptr, c_size_t

    implicit none

    type(c_ptr), value :: p_h
    type(c_devptr), value :: p_d
    integer(c_size_t) :: sz
    integer :: idx

    integer :: cudaResult

    cudaResult = cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, cuda_streams(stream_from_index(idx)))
    call gpu_error_test(cudaResult)

  end subroutine gpu_dtoh_memcpy_async



  subroutine gpu_htod_memprefetch_async(p, sz, idx) bind(c, name='gpu_htod_memprefetch_async')

    use cudafor, only: cudaMemPrefetchAsync, c_devptr
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: p
    integer(c_size_t) :: sz
    integer :: idx

    integer :: cudaResult

    if ((.not. have_prop) .or. (have_prop .and. prop%managedMemory == 1 .and. prop%concurrentManagedAccess == 1)) then

       cudaResult = cudaMemPrefetchAsync(p, sz, cuda_device_id, cuda_streams(stream_from_index(idx)))
       call gpu_error_test(cudaResult)

    end if

  end subroutine gpu_htod_memprefetch_async



  subroutine gpu_dtoh_memprefetch_async(p, sz, idx) bind(c, name='gpu_dtoh_memprefetch_async')

    use cudafor, only: cudaMemPrefetchAsync, c_devptr, cudaCpuDeviceId
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: p
    integer(c_size_t) :: sz
    integer :: idx

    integer :: cudaResult

    if ((.not. have_prop) .or. (have_prop .and. prop%managedMemory == 1)) then

       cudaResult = cudaMemPrefetchAsync(p, sz, cudaCpuDeviceId, cuda_streams(stream_from_index(idx)))
       call gpu_error_test(cudaResult)

    end if

  end subroutine gpu_dtoh_memprefetch_async



  subroutine gpu_synchronize() bind(c, name='gpu_synchronize')

    use cudafor, only: cudaDeviceSynchronize

    implicit none

    integer :: cudaResult

    cudaResult = cudaDeviceSynchronize()
    call gpu_error_test(cudaResult)

  end subroutine gpu_synchronize



  subroutine gpu_stream_synchronize(index) bind(c, name='gpu_stream_synchronize')

    use cudafor, only: cudaStreamSynchronize

    implicit none

    integer, intent(in), value :: index

    integer :: cudaResult

    cudaResult = cudaStreamSynchronize(cuda_streams(stream_from_index(index)))
    call gpu_error_test(cudaResult)

  end subroutine gpu_stream_synchronize



  subroutine mem_advise_set_preferred(p, sz, device) bind(c, name='mem_advise_set_preferred')

    use cudafor, only: c_devptr, cudaMemAdvise, cudaMemAdviseSetPreferredLocation
    use iso_c_binding, only: c_size_t, c_int

    type(c_devptr), value :: p
    integer(c_size_t), value :: sz
    integer(c_int) :: device

    integer :: cudaResult

    ! Note: we do not error test in this subroutine because the error
    ! code seems to be broken in PGI.

    if ((.not. have_prop) .or. (have_prop .and. prop%concurrentManagedAccess == 1)) then
       cudaResult = cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device)
    end if

  end subroutine mem_advise_set_preferred



  subroutine mem_advise_set_readonly(p, sz) bind(c, name='mem_advise_set_readonly')

    use cudafor, only: c_devptr, cudaMemAdvise, cudaMemAdviseSetReadMostly, cudaCpuDeviceId
    use iso_c_binding, only: c_size_t

    type(c_devptr), value :: p
    integer(c_size_t), value :: sz

    integer :: cudaResult

    ! Note: we do not error test in this subroutine because the error
    ! code seems to be broken in PGI.

    ! Note: the device argument in this call is ignored, so we arbitrarily pick the CPU.
    if ((.not. have_prop) .or. (have_prop .and. prop%concurrentManagedAccess == 1)) then
       cudaResult = cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId)
    end if

  end subroutine mem_advise_set_readonly



  subroutine gpu_error(cudaResult, abort) bind(c, name='gpu_error')

    integer, intent(in) :: cudaResult
    logical, intent(in), optional :: abort

    character(len=32) :: cudaResultStr
    character(len=128) :: error_string
    logical :: do_abort

    do_abort = .true.

    if (present(abort)) then
       do_abort = abort
    end if

    write(cudaResultStr, "(I32)") cudaResult

    error_string = "CUDA Error " // trim(adjustl(cudaResultStr)) // ": " // cudaGetErrorString(cudaResult)

    if (abort) then
       call bl_abort(error_string)
    else
       call bl_warning(error_string)
    end if

  end subroutine gpu_error



  subroutine gpu_error_test(cudaResult, abort) bind(c, name='gpu_error_test')

    implicit none

    integer, intent(in) :: cudaResult
    logical, intent(in), optional :: abort

    logical :: do_abort

    do_abort = .true.

    if (present(abort)) then
       do_abort = abort
    end if

    if (cudaResult /= cudaSuccess) then
       call gpu_error(cudaResult, do_abort)
    end if

  end subroutine gpu_error_test



  ! Check if any kernels or previous API calls have returned an error code.
  ! Abort if so.

  subroutine check_for_gpu_errors() bind(c, name='check_for_gpu_errors')

    use cudafor, only: cudaGetLastError

    implicit none

    integer :: cudaResult

    cudaResult = cudaGetLastError()

    call gpu_error_test(cudaResult)

  end subroutine check_for_gpu_errors

end module cuda_module
