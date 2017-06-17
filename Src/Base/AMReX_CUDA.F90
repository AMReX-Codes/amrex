module cuda_module

  use cudafor, only: cuda_stream_kind

  implicit none

  integer, parameter :: max_cuda_streams = 100
  integer(kind=cuda_stream_kind) :: cuda_streams(max_cuda_streams)

  integer, save :: cuda_device_id

contains

  subroutine initialize_cuda() bind(c, name='initialize_cuda')

    use cudafor, only: cudaStreamCreate
    
    implicit none

    integer :: i, cudaResult

    do i = 1, max_cuda_streams
       cudaResult = cudaStreamCreate(cuda_streams(i))
    enddo

    cuda_device_id = 0

  end subroutine initialize_cuda



  subroutine get_cuda_device_id(id) bind(c, name='get_cuda_device_id')

    implicit none

    integer :: id

    id = cuda_device_id

  end subroutine get_cuda_device_id



  integer function stream_from_index(idx)

    implicit none

    integer :: idx

    if (idx < 0) then
       stream_from_index = 0
    else
       stream_from_index = MOD(idx, max_cuda_streams) + 1
    endif

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

  end subroutine gpu_malloc



  subroutine gpu_malloc_managed(x, sz) bind(c, name='gpu_malloc_managed')

    use cudafor, only: cudaMallocManaged, cudaMemAttachGlobal, c_devptr
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: x
    integer(c_size_t) :: sz

    integer :: cudaResult

    cudaResult = cudaMallocManaged(x, sz, cudaMemAttachGlobal)

  end subroutine gpu_malloc_managed



  subroutine gpu_free(x) bind(c, name='gpu_free')

    use cudafor, only: cudaFree, c_devptr

    implicit none

    type(c_devptr), value :: x

    integer :: cudaResult

    cudaResult = cudaFree(x)

  end subroutine gpu_free



  subroutine gpu_htod_memcpy_async(p_d, p_h, sz, idx) bind(c, name='gpu_htod_memcpy_async')

    use cudafor, only: cudaMemcpyAsync, cudaMemcpyHostToDevice, c_devptr, cuda_stream_kind
    use iso_c_binding, only: c_ptr, c_size_t

    implicit none

    type(c_devptr), value :: p_d
    type(c_ptr), value :: p_h
    integer(c_size_t) :: sz
    integer :: idx

    integer :: s
    integer :: cudaResult

    if (idx < 0) then

       cudaResult = cudaMemcpyAsync(p_d, p_h, sz)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, cuda_streams(s))

    endif

  end subroutine gpu_htod_memcpy_async



  subroutine gpu_dtoh_memcpy_async(p_h, p_d, sz, idx) bind(c, name='gpu_dtoh_memcpy_async')

    use cudafor, only: cudaMemcpyAsync, cudaMemcpyDeviceToHost, c_devptr
    use iso_c_binding, only: c_ptr, c_size_t

    implicit none

    type(c_ptr), value :: p_h
    type(c_devptr), value :: p_d
    integer(c_size_t) :: sz
    integer :: idx

    integer :: s
    integer :: cudaResult

    if (idx < 0) then

       cudaResult = cudaMemcpyAsync(p_h, p_d, sz)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, cuda_streams(s))

    endif

  end subroutine gpu_dtoh_memcpy_async



  subroutine gpu_htod_memprefetch_async(p, sz, idx) bind(c, name='gpu_htod_memprefetch_async')

    use cudafor, only: cudaMemPrefetchAsync, c_devptr
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: p
    integer(c_size_t) :: sz
    integer :: idx

    integer :: s
    integer :: cudaResult

    if (idx < 0) then

       cudaResult = cudaMemPrefetchAsync(p, sz, cuda_device_id, 0)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemPrefetchAsync(p, sz, cuda_device_id, cuda_streams(s))

    endif

  end subroutine gpu_htod_memprefetch_async



  subroutine gpu_dtoh_memprefetch_async(p, sz, idx) bind(c, name='gpu_dtoh_memprefetch_async')

    use cudafor, only: cudaMemPrefetchAsync, c_devptr, cudaCpuDeviceId
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: p
    integer(c_size_t) :: sz
    integer :: idx

    integer :: s
    integer :: cudaResult

    if (idx < 0) then

       cudaResult = cudaMemPrefetchAsync(p, sz, cudaCpuDeviceId, 0)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemPrefetchAsync(p, sz, cudaCpuDeviceId, cuda_streams(s))

    endif

  end subroutine gpu_dtoh_memprefetch_async



  subroutine gpu_synchronize() bind(c, name='gpu_synchronize')

    use cudafor, only: cudaDeviceSynchronize

    implicit none

    integer :: cudaResult

    cudaResult = cudaDeviceSynchronize();

  end subroutine gpu_synchronize



  subroutine mem_advise_set_preferred(p, sz, device) bind(c, name='mem_advise_set_preferred')

    use cudafor, only: c_devptr, cudaMemAdvise, cudaMemAdviseSetPreferredLocation
    use iso_c_binding, only: c_size_t, c_int

    type(c_devptr) :: p
    integer(c_size_t) :: sz
    integer(c_int) :: device

    integer :: cudaResult

    cudaResult = cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device)

  end subroutine mem_advise_set_preferred



  subroutine mem_advise_set_readonly(p, sz) bind(c, name='mem_advise_set_readonly')

    use cudafor, only: c_devptr, cudaMemAdvise, cudaMemAdviseSetReadMostly, cudaCpuDeviceId
    use iso_c_binding, only: c_size_t

    type(c_devptr) :: p
    integer(c_size_t) :: sz

    integer :: cudaResult

    ! Note: the device argument in this call is ignored, so we arbitrarily pick the CPU.
    cudaResult = cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId)

  end subroutine mem_advise_set_readonly

end module cuda_module
