module cuda_module

  use cudafor, only: cuda_stream_kind, cudaSuccess, cudaGetErrorString

  implicit none

  integer, parameter :: max_cuda_streams = 100
  integer(kind=cuda_stream_kind) :: cuda_streams(max_cuda_streams)

  integer, save :: cuda_device_id

contains

  subroutine initialize_cuda() bind(c, name='initialize_cuda')

    use cudafor, only: cudaStreamCreate
    
    implicit none

    integer :: i, cudaResult
    character(len=32) :: cudaResultStr

    do i = 1, max_cuda_streams
       cudaResult = cudaStreamCreate(cuda_streams(i))

       if (cudaResult /= cudaSuccess) then
          write(cudaResultStr, "(I32)") cudaResult
          call bl_abort("CUDA failure in initialize_cuda(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
       endif
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
    character(len=32) :: cudaResultStr

    cudaResult = cudaMalloc(x, sz)

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_malloc(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_malloc



  subroutine gpu_hostalloc(x, sz) bind(c, name='gpu_hostalloc')

    use cudafor, only: cudaHostAlloc, cudaHostAllocMapped, cudaHostAllocWriteCombined
    use iso_c_binding, only: c_size_t, c_ptr

    implicit none

    type(c_ptr) :: x
    integer :: sz

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    cudaResult = cudaHostAlloc(x, sz, ior(cudaHostAllocMapped, cudaHostAllocWriteCombined))

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_hostalloc(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_hostalloc



  subroutine gpu_malloc_managed(x, sz) bind(c, name='gpu_malloc_managed')

    use cudafor, only: cudaMallocManaged, cudaMemAttachGlobal, c_devptr
    use iso_c_binding, only: c_size_t

    implicit none

    type(c_devptr) :: x
    integer(c_size_t) :: sz

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    cudaResult = cudaMallocManaged(x, sz, cudaMemAttachGlobal)

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_malloc_managed(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_malloc_managed



  subroutine gpu_free(x) bind(c, name='gpu_free')

    use cudafor, only: cudaFree, c_devptr

    implicit none

    type(c_devptr), value :: x

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    cudaResult = cudaFree(x)

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_free(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_free



  subroutine gpu_freehost(x) bind(c, name='gpu_freehost')

    use cudafor, only: cudaFreeHost
    use iso_c_binding, only: c_ptr

    implicit none

    type(c_ptr), value :: x

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    cudaResult = cudaFreeHost(x)

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_freehost(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_freehost



  subroutine gpu_host_device_ptr(x, y) bind(c, name='gpu_host_device_ptr')

    use cudafor, only: cudaHostGetDevicePointer, c_devptr
    use iso_c_binding, only: c_ptr

    type(c_devptr) :: x
    type(c_ptr), value :: y
    integer :: flags

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    flags = 0 ! This argument does nothing at present, but is part of the API.

    cudaResult = cudaHostGetDevicePointer(x, y, flags)

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_host_device_ptr(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_host_device_ptr



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
    character(len=32) :: cudaResultStr

    if (idx < 0) then

       cudaResult = cudaMemcpyAsync(p_d, p_h, sz)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemcpyAsync(p_d, p_h, sz, cudaMemcpyHostToDevice, cuda_streams(s))

    endif

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_htod_memcpy_async(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
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
    character(len=32) :: cudaResultStr

    if (idx < 0) then

       cudaResult = cudaMemcpyAsync(p_h, p_d, sz)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemcpyAsync(p_h, p_d, sz, cudaMemcpyDeviceToHost, cuda_streams(s))

    endif

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_dtoh_memcpy_async(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
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
    character(len=32) :: cudaResultStr

    if (idx < 0) then

       cudaResult = cudaMemPrefetchAsync(p, sz, cuda_device_id, 0)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemPrefetchAsync(p, sz, cuda_device_id, cuda_streams(s))

    endif

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_htod_memprefetch_async(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
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
    character(len=32) :: cudaResultStr

    if (idx < 0) then

       cudaResult = cudaMemPrefetchAsync(p, sz, cudaCpuDeviceId, 0)

    else

       s = stream_from_index(idx)

       cudaResult = cudaMemPrefetchAsync(p, sz, cudaCpuDeviceId, cuda_streams(s))

    endif

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_dtoh_memprefetch_async(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_dtoh_memprefetch_async



  subroutine gpu_synchronize() bind(c, name='gpu_synchronize')

    use cudafor, only: cudaDeviceSynchronize

    implicit none

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    cudaResult = cudaDeviceSynchronize();

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in gpu_synchronize(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine gpu_synchronize



  subroutine mem_advise_set_preferred(p, sz, device) bind(c, name='mem_advise_set_preferred')

    use cudafor, only: c_devptr, cudaMemAdvise, cudaMemAdviseSetPreferredLocation
    use iso_c_binding, only: c_size_t, c_int

    type(c_devptr), value :: p
    integer(c_size_t) :: sz
    integer(c_int) :: device

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    cudaResult = cudaMemAdvise(p, sz, cudaMemAdviseSetPreferredLocation, device)

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in mem_advise_set_preferred(), Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine mem_advise_set_preferred



  subroutine mem_advise_set_readonly(p, sz) bind(c, name='mem_advise_set_readonly')

    use cudafor, only: c_devptr, cudaMemAdvise, cudaMemAdviseSetReadMostly, cudaCpuDeviceId
    use iso_c_binding, only: c_size_t

    type(c_devptr), value :: p
    integer(c_size_t) :: sz

    integer :: cudaResult
    character(len=32) :: cudaResultStr

    ! Note: the device argument in this call is ignored, so we arbitrarily pick the CPU.
    cudaResult = cudaMemAdvise(p, sz, cudaMemAdviseSetReadMostly, cudaCpuDeviceId)

    if (cudaResult /= cudaSuccess) then
       write(cudaResultStr, "(I32)") cudaResult
       call bl_abort("CUDA failure in mem_advise_set_readonly(): Error " // trim(adjustl(cudaResultStr)) // ", " // cudaGetErrorString(cudaResult))
    endif

  end subroutine mem_advise_set_readonly

end module cuda_module
