module basefab_nd_module

  use amrex_fort_module, only: amrex_real, get_loop_bounds
#ifdef CUDA
  use cuda_module, only: stream_from_index, cuda_streams
  use cudafor, only: cuda_stream_kind
#endif

  implicit none

contains

  ! dst = src
#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_copy_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    integer, intent(in), value :: ncomp

    integer :: i,j,k,n,off(3)
    integer :: blo(3), bhi(3)

    off = sblo - lo

    call get_loop_bounds(blo, bhi, lo, hi)

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                dst(i,j,k,n) = src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_copy_doit


  ! copy from multi-d array to 1d array
  function fort_fab_copytomem (lo, hi, dst, src, slo, shi, ncomp) result(nelems) &
       bind(c,name='fort_fab_copytomem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real)             :: dst(*)
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)

    integer :: i, j, k, n, nx
    integer(c_long) :: offset

    nx = hi(1)-lo(1)+1
    offset = 1-lo(1)
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(offset+i) = src(i,j,k,n) 
             end do
             offset = offset + nx
          end do
       end do
    end do    

    nelems = offset - (1-lo(1))
  end function fort_fab_copytomem


  ! copy from 1d array to multi-d array
  function fort_fab_copyfrommem (lo, hi, dst, dlo, dhi, ncomp, src) result(nelems) &
       bind(c,name='fort_fab_copyfrommem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(amrex_real), intent(in   ) :: src(*)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n, nx
    integer(c_long) :: offset

    nx = hi(1)-lo(1)+1
    offset = 1-lo(1)
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n)  = src(offset+i)
             end do
             offset = offset + nx
          end do
       end do
    end do    

    nelems = offset - (1-lo(1))
  end function fort_fab_copyfrommem



#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_setval_doit(lo, hi, dst, dlo, dhi, ncomp, val)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in), value :: val
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                dst(i,j,k,n) = val
             end do
          end do
       end do
    end do

  end subroutine fort_fab_setval_doit


  function fort_fab_norm_doit (lo, hi, src, slo, shi, ncomp, p) result(nrm)
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp, p
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real) :: nrm

    integer :: i,j,k,n

#ifdef CUDA
    attributes(device) :: src
#endif

    nrm = 0.0_amrex_real
    if (p .eq. 0) then ! max norm
       !$cuf kernel do(4) <<<*, 256>>>
       do n = 1, ncomp
          do       k = lo(3), hi(3)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   nrm = max(nrm, abs(src(i,j,k,n)))
                end do
             end do
          end do
       end do
    else if (p .eq. 1) then
       !$cuf kernel do(4) <<<*, 256>>>
       do n = 1, ncomp
          do       k = lo(3), hi(3)
             do    j = lo(2), hi(2)
                do i = lo(1), hi(1)
                   nrm = nrm + abs(src(i,j,k,n))
                end do
             end do
          end do
       end do
    end if
  end function fort_fab_norm_doit


  function fort_fab_sum_doit(lo, hi, src, slo, shi, ncomp) result(sm)
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real) :: sm

    integer :: i,j,k,n

#ifdef CUDA
    attributes(device) :: src
#endif

    sm = 0.0_amrex_real
    !$cuf kernel do(4) <<<*, 256>>>
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                sm = sm + src(i,j,k,n)
             end do
          end do
       end do
    end do
  end function fort_fab_sum_doit


#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_plus_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n,off(3)
    integer :: blo(3), bhi(3)

    off = sblo - lo

    call get_loop_bounds(blo, bhi, lo, hi)

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_plus_doit


  subroutine fort_fab_minus_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp, index
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n,off(3)

#ifdef CUDA
    attributes(device) :: src, dst, off
    integer(cuda_stream_kind) :: stream
    stream = cuda_streams(stream_from_index(index))
#endif

    off = sblo - lo

    !$cuf kernel do(4) <<<*, 256, 0, stream>>>
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) - src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_minus_doit


  subroutine fort_fab_mult_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp, index
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n,off(3)

#ifdef CUDA
    attributes(device) :: src, dst, off
    integer(cuda_stream_kind) :: stream
    stream = cuda_streams(stream_from_index(index))
#endif

    off = sblo - lo

    !$cuf kernel do(4) <<<*, 256, 0, stream>>>
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) * src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_mult_doit


  subroutine fort_fab_divide_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp, index
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n,off(3)

#ifdef CUDA
    attributes(device) :: src, dst, off
    integer(cuda_stream_kind) :: stream
    stream = cuda_streams(stream_from_index(index))
#endif

    off = sblo - lo

    !$cuf kernel do(4) <<<*, 256, 0, stream>>>
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) / src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_divide_doit


  subroutine fort_fab_protdivide(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='fort_fab_protdivide')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                if (src(i+off(1),j+off(2),k+off(3),n) .ne. 0._amrex_real) then
                   dst(i,j,k,n) = dst(i,j,k,n) / src(i+off(1),j+off(2),k+off(3),n)
                end if
             end do
          end do
       end do
    end do
  end subroutine fort_fab_protdivide


  ! dst = a/src
  subroutine fort_fab_invert(lo, hi, dst, dlo, dhi, ncomp, a) &
       bind(c,name='fort_fab_invert')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(amrex_real), intent(in   ) :: a
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = a / dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_invert


  ! dst += a*src
#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_saxpy_doit(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n,off(3)
    integer :: blo(3), bhi(3)

    off = sblo - lo

    call get_loop_bounds(blo, bhi, lo, hi)

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + a * src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_saxpy_doit


  ! dst = src + a*dst
  subroutine fort_fab_xpay(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp) &
       bind(c,name='fort_fab_xpay')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(amrex_real), intent(in   ) :: a
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n,off(3)

    off = sblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = src(i+off(1),j+off(2),k+off(3),n) + a * dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_xpay


  ! dst = a*x + b*y
  subroutine fort_fab_lincomb(lo, hi, dst, dlo, dhi, a, x, xlo, xhi, xblo, &
       b, y, ylo, yhi, yblo, ncomp) bind(c,name='fort_fab_lincomb')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), xlo(3), xhi(3), xblo(3), &
         ylo(3), yhi(3), yblo(3), ncomp
    real(amrex_real), intent(in   ) :: a, b
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)

    integer :: i,j,k,n,xoff(3),yoff(3)

    xoff = xblo - lo
    yoff = yblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = a * x(i+xoff(1),j+xoff(2),k+xoff(3),n) &
                     +         b * y(i+yoff(1),j+yoff(2),k+yoff(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_lincomb

  ! dst = dst + src1*src2
  subroutine fort_fab_addproduct(lo, hi, dst, dlo, dhi, src1, s1lo, s1hi, src2, s2lo, s2hi,ncomp) &
       bind(c,name='fort_fab_addproduct')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), s1lo(3), s1hi(3), s2lo(3), s2hi(3), ncomp
    real(amrex_real), intent(in   ) :: src1(s1lo(1):s1hi(1),s1lo(2):s1hi(2),s1lo(3):s1hi(3),ncomp)
    real(amrex_real), intent(in   ) :: src2(s2lo(1):s2hi(1),s2lo(2):s2hi(2),s2lo(3):s2hi(3),ncomp)
    real(amrex_real), intent(inout) ::  dst( dlo(1): dhi(1), dlo(2): dhi(2), dlo(3): dhi(3),ncomp)

    integer :: i,j,k,n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = src1(i,j,k,n) * src2(i,j,k,n) + dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_addproduct

  ! dot_product
  function fort_fab_dot(lo, hi, x, xlo, xhi, y, ylo, yhi, yblo, ncomp) result(dp) &
       bind(c,name='fort_fab_dot')
    integer, intent(in) :: lo(3), hi(3), xlo(3), xhi(3), ylo(3), yhi(3), yblo(3), ncomp
    real(amrex_real), intent(in) :: x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in) :: y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)
    real(amrex_real) :: dp

    integer :: i,j,k,n, off(3)

    dp = 0.0_amrex_real

    off = yblo - lo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dp = dp + x(i,j,k,n)*y(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end function fort_fab_dot

end module basefab_nd_module


  subroutine fort_fab_copy(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index) &
                           bind(c,name='fort_fab_copy')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_copy_doit
#ifdef CUDA
  use cuda_module, only: stream_from_index, cuda_streams, threads_and_blocks
  use cudafor, only: cuda_stream_kind, dim3
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), index
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    integer, intent(in), value :: ncomp

#ifdef CUDA
    attributes(managed) :: src, dst, lo, hi, dlo, dhi, slo, shi, sblo

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(stream_from_index(index))

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call fort_fab_copy_doit &
#ifdef CUDA
    <<<numBlocks, numThreads, 0, stream>>> &
#endif
    (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_copy



  subroutine fort_fab_setval(lo, hi, dst, dlo, dhi, ncomp, val, index) &
                             bind(c,name='fort_fab_setval')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_setval_doit
#ifdef CUDA
  use cuda_module, only: stream_from_index, cuda_streams, threads_and_blocks
  use cudafor, only: cuda_stream_kind, dim3
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), index
    real(amrex_real), intent(in), value :: val
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    integer, intent(in), value :: ncomp

#ifdef CUDA
    attributes(managed) :: dst, lo, hi, dlo, dhi

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(stream_from_index(index))

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call fort_fab_setval_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, ncomp, val)

  end subroutine fort_fab_setval




  function fort_fab_norm (lo, hi, src, slo, shi, ncomp, p) result(nrm) &
       bind(c,name='fort_fab_norm')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_norm_doit

    implicit none

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp, p
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real) :: nrm

#ifdef CUDA
    attributes(device) :: src
#endif

    nrm = fort_fab_norm_doit(lo, hi, src, slo, shi, ncomp, p)

  end function fort_fab_norm



  subroutine fort_fab_saxpy(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp, index) bind(c,name='fort_fab_saxpy')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_saxpy_doit
#ifdef CUDA
  use cuda_module, only: stream_from_index, cuda_streams, threads_and_blocks
  use cudafor, only: cuda_stream_kind, dim3
#endif

    implicit none
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), index
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(managed) :: src, dst, lo, hi, dlo, dhi, slo, shi, sblo

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(stream_from_index(index))

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call fort_fab_saxpy_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_saxpy



  subroutine fort_fab_plus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index) &
                           bind(c,name='fort_fab_plus')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_plus_doit
#ifdef CUDA
    use cuda_module, only: stream_from_index, cuda_streams, threads_and_blocks
    use cudafor, only: cuda_stream_kind, dim3
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), index
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(managed) :: src, dst, lo, hi, dlo, dhi, slo, shi, sblo

    integer :: cuda_result
    integer(kind=cuda_stream_kind) :: stream
    type(dim3) :: numThreads, numBlocks

    stream = cuda_streams(stream_from_index(index))

    call threads_and_blocks(lo, hi, numBlocks, numThreads)
#endif

    call fort_fab_plus_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_plus



  function fort_fab_sum(lo, hi, src, slo, shi, ncomp) result(sm) &
                        bind(c,name='fort_fab_sum')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_sum_doit

    implicit none

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real) :: sm

#ifdef CUDA
    attributes(device) :: src
#endif

    sm = fort_fab_sum_doit(lo, hi, src, slo, shi, ncomp)

  end function fort_fab_sum



  subroutine fort_fab_minus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index) &
                            bind(c,name='fort_fab_minus')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_minus_doit

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp, index
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: src, dst
#endif

    call fort_fab_minus_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index)

  end subroutine fort_fab_minus



  subroutine fort_fab_mult(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index) &
                           bind(c,name='fort_fab_mult')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_mult_doit

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp, index
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: src, dst
#endif

    call fort_fab_mult_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index)

  end subroutine fort_fab_mult


  subroutine fort_fab_divide(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index) &
                             bind(c,name='fort_fab_divide')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_divide_doit

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp, index
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: src, dst
#endif

    call fort_fab_divide_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp, index)

  end subroutine fort_fab_divide
