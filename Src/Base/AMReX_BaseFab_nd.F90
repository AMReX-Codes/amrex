module basefab_nd_module

  use amrex_fort_module, only: amrex_real, get_loop_bounds, amrex_add, amrex_max

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
  subroutine fort_fab_copytomem (lo, hi, dst, src, slo, shi, ncomp) &
       bind(c,name='fort_fab_copytomem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real)             :: dst(*)
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)

    integer :: i, j, k, n
    integer(c_long) :: offset, tile_size(3)

    tile_size(1:3) = hi - lo + 1

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             offset = tile_size(1) * (j - lo(2)) + tile_size(1) * tile_size(2) * (k - lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1)
             do i = lo(1), hi(1)
                dst(offset+i) = src(i,j,k,n)
             end do
          end do
       end do
    end do

  end subroutine fort_fab_copytomem


  ! copy from 1d array to multi-d array
  subroutine fort_fab_copyfrommem (lo, hi, dst, dlo, dhi, ncomp, src) &
       bind(c,name='fort_fab_copyfrommem')
    use iso_c_binding, only : c_long
    integer(c_long) :: nelems
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(amrex_real), intent(in   ) :: src(*)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n, nx
    integer(c_long) :: offset, tile_size(3)

    tile_size = hi - lo + 1

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             offset = tile_size(1) * (j - lo(2)) + tile_size(1) * tile_size(2) * (k - lo(3)) + &
                      tile_size(1) * tile_size(2) * tile_size(3) * (n - 1)
             do i = lo(1), hi(1)
                dst(i,j,k,n)  = src(offset+i)
             end do
          end do
       end do
    end do

  end subroutine fort_fab_copyfrommem



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



#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_norm_doit (lo, hi, src, slo, shi, ncomp, p, nrm)
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp, p
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: nrm

    integer :: i,j,k,n
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    nrm = 0.0_amrex_real
    if (p .eq. 0) then ! max norm
       do n = 1, ncomp
          do       k = blo(3), bhi(3)
             do    j = blo(2), bhi(2)
                do i = blo(1), bhi(1)
                   call amrex_max(nrm, abs(src(i,j,k,n)))
                end do
             end do
          end do
       end do
    else if (p .eq. 1) then
       do n = 1, ncomp
          do       k = blo(3), bhi(3)
             do    j = blo(2), bhi(2)
                do i = blo(1), bhi(1)
                   call amrex_add(nrm, abs(src(i,j,k,n)))
                end do
             end do
          end do
       end do
    end if
  end subroutine fort_fab_norm_doit



#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_sum_doit(lo, hi, src, slo, shi, ncomp, sm)
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: sm

    integer :: i,j,k,n
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    sm = 0.0_amrex_real
    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                call amrex_add(sm, src(i,j,k,n))
             end do
          end do
       end do
    end do
  end subroutine fort_fab_sum_doit



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


#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_minus_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)
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
                dst(i,j,k,n) = dst(i,j,k,n) - src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_minus_doit


#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_mult_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)
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
                dst(i,j,k,n) = dst(i,j,k,n) * src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_mult_doit


#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_divide_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)
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
                dst(i,j,k,n) = dst(i,j,k,n) / src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_divide_doit


#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_protdivide_doit(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)
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
                if (src(i+off(1),j+off(2),k+off(3),n) .ne. 0._amrex_real) then
                   dst(i,j,k,n) = dst(i,j,k,n) / src(i+off(1),j+off(2),k+off(3),n)
                end if
             end do
          end do
       end do
    end do
  end subroutine fort_fab_protdivide_doit


  ! dst = a/src
#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_invert_doit(lo, hi, dst, dlo, dhi, ncomp, a)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i,j,k,n
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                dst(i,j,k,n) = a / dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_invert_doit


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
#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_xpay_doit(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp)
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
                dst(i,j,k,n) = src(i+off(1),j+off(2),k+off(3),n) + a * dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_xpay_doit


  ! dst = a*x + b*y
#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_lincomb_doit(lo, hi, dst, dlo, dhi, a, x, xlo, xhi, xblo, &
                                   b, y, ylo, yhi, yblo, ncomp)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), xlo(3), xhi(3), xblo(3), &
                           ylo(3), yhi(3), yblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a, b
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)

    integer :: i,j,k,n,xoff(3),yoff(3)
    integer :: blo(3), bhi(3)

    xoff = xblo - lo
    yoff = yblo - lo

    call get_loop_bounds(blo, bhi, lo, hi)

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                dst(i,j,k,n) = a * x(i+xoff(1),j+xoff(2),k+xoff(3),n) &
                     +         b * y(i+yoff(1),j+yoff(2),k+yoff(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_lincomb_doit

  ! dst = dst + src1*src2
#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_addproduct_doit(lo, hi, dst, dlo, dhi, src1, s1lo, s1hi, src2, s2lo, s2hi, ncomp)
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), s1lo(3), s1hi(3), s2lo(3), s2hi(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src1(s1lo(1):s1hi(1),s1lo(2):s1hi(2),s1lo(3):s1hi(3),ncomp)
    real(amrex_real), intent(in   ) :: src2(s2lo(1):s2hi(1),s2lo(2):s2hi(2),s2lo(3):s2hi(3),ncomp)
    real(amrex_real), intent(inout) ::  dst( dlo(1): dhi(1), dlo(2): dhi(2), dlo(3): dhi(3),ncomp)

    integer :: i,j,k,n
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                dst(i,j,k,n) = src1(i,j,k,n) * src2(i,j,k,n) + dst(i,j,k,n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_addproduct_doit

  ! dot_product
#ifdef CUDA
  attributes(global) &
#endif
  subroutine fort_fab_dot_doit(lo, hi, x, xlo, xhi, y, ylo, yhi, yblo, ncomp, dp)
    integer, intent(in) :: lo(3), hi(3), xlo(3), xhi(3), ylo(3), yhi(3), yblo(3), ncomp
    real(amrex_real), intent(in) :: x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in) :: y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)
    real(amrex_real) :: dp

    integer :: i,j,k,n, off(3)
    integer :: blo(3), bhi(3)

    call get_loop_bounds(blo, bhi, lo, hi)

    dp = 0.0_amrex_real

    off = yblo - lo

    do n = 1, ncomp
       do       k = blo(3), bhi(3)
          do    j = blo(2), bhi(2)
             do i = blo(1), bhi(1)
                call amrex_add(dp, x(i,j,k,n)*y(i+off(1),j+off(2),k+off(3),n))
             end do
          end do
       end do
    end do
  end subroutine fort_fab_dot_doit

end module basefab_nd_module


  subroutine fort_fab_copy(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
                           bind(c,name='fort_fab_copy')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_copy_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    integer, intent(in), value :: ncomp

#ifdef CUDA
    attributes(managed) :: src, dst, lo, hi, dlo, dhi, slo, shi, sblo
#endif

    call fort_fab_copy_doit &
#ifdef CUDA
    <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
    (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_copy



  subroutine fort_fab_setval(lo, hi, dst, dlo, dhi, ncomp, val) &
                             bind(c,name='fort_fab_setval')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_setval_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    real(amrex_real), intent(in), value :: val
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    integer, intent(in), value :: ncomp

#ifdef CUDA
    attributes(managed) :: dst, lo, hi, dlo, dhi
#endif

    call fort_fab_setval_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, ncomp, val)

  end subroutine fort_fab_setval




  subroutine fab_fort_norm (lo, hi, src, slo, shi, ncomp, p, nrm) &
       bind(c,name='fort_fab_norm')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_norm_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp, p
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: nrm

#ifdef CUDA
    attributes(device) :: lo, hi, src, slo, shi, ncomp, p, nrm
#endif

    call fort_fab_norm_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, src, slo, shi, ncomp, p, nrm)

  end subroutine fab_fort_norm



  subroutine fort_fab_saxpy(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp) bind(c,name='fort_fab_saxpy')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_saxpy_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(managed) :: src, dst, lo, hi, dlo, dhi, slo, shi, sblo
#endif

    call fort_fab_saxpy_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_saxpy



  subroutine fort_fab_plus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
                           bind(c,name='fort_fab_plus')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_plus_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(managed) :: src, dst, lo, hi, dlo, dhi, slo, shi, sblo
#endif

    call fort_fab_plus_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_plus



  subroutine fort_fab_sum(lo, hi, src, slo, shi, ncomp, sm) &
                        bind(c,name='fort_fab_sum')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_sum_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(amrex_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: sm

#ifdef CUDA
    attributes(device) :: lo, hi, src, slo, shi, ncomp, sm
#endif

    call fort_fab_sum_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, src, slo, shi, ncomp, sm)

  end subroutine fort_fab_sum



  subroutine fort_fab_minus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
                            bind(c,name='fort_fab_minus')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_minus_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi, src, slo, shi, sblo
#endif

    call fort_fab_minus_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_minus



  subroutine fort_fab_mult(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
                           bind(c,name='fort_fab_mult')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_mult_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi, src, slo, shi, sblo
#endif

    call fort_fab_mult_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_mult


  subroutine fort_fab_divide(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
                             bind(c,name='fort_fab_divide')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_divide_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi, src, slo, shi, sblo
#endif

    call fort_fab_divide_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_divide



  subroutine fort_fab_protdivide(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
                                 bind(c,name='fort_fab_protdivide')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_protdivide_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi, src, slo, shi, sblo
#endif

    call fort_fab_protdivide_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_protdivide



  subroutine fort_fab_invert(lo, hi, dst, dlo, dhi, ncomp, a) bind(c, name='fort_fab_invert')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_invert_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi
#endif

    call fort_fab_invert_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, ncomp, a)

  end subroutine fort_fab_invert



  subroutine fort_fab_xpay(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp) bind(c, name='fort_fab_xpay')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_xpay_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a
    real(amrex_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi, src, slo, shi, sblo
#endif

    call fort_fab_xpay_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp)

  end subroutine fort_fab_xpay



  subroutine fort_fab_lincomb(lo, hi, dst, dlo, dhi, a, x, xlo, xhi, xblo, &
                              b, y, ylo, yhi, yblo, ncomp) bind(c, name='fort_fab_lincomb')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_lincomb_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), xlo(3), xhi(3), xblo(3), &
                           ylo(3), yhi(3), yblo(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ), value :: a, b
    real(amrex_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in   ) ::   y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi, x, xlo, xhi, xblo, y, ylo, yhi, yblo
#endif

    call fort_fab_lincomb_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, a, x, xlo, xhi, xblo, b, y, ylo, yhi, yblo, ncomp)

  end subroutine fort_fab_lincomb



  subroutine fort_fab_addproduct(lo, hi, dst, dlo, dhi, src1, s1lo, s1hi, src2, s2lo, s2hi, ncomp) &
       bind(c, name='fort_fab_addproduct')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_addproduct_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), s1lo(3), s1hi(3), s2lo(3), s2hi(3)
    integer, intent(in), value :: ncomp
    real(amrex_real), intent(in   ) :: src1(s1lo(1):s1hi(1),s1lo(2):s1hi(2),s1lo(3):s1hi(3),ncomp)
    real(amrex_real), intent(in   ) :: src2(s2lo(1):s2hi(1),s2lo(2):s2hi(2),s2lo(3):s2hi(3),ncomp)
    real(amrex_real), intent(inout) ::  dst( dlo(1): dhi(1), dlo(2): dhi(2), dlo(3): dhi(3),ncomp)

#ifdef CUDA
    attributes(device) :: lo, hi, dst, dlo, dhi, src1, s1lo, s1hi, src2, s2lo, s2hi
#endif

    call fort_fab_addproduct_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, dst, dlo, dhi, src1, s1lo, s1hi, src2, s2lo, s2hi, ncomp)

  end subroutine fort_fab_addproduct



  subroutine fort_fab_dot(lo, hi, x, xlo, xhi, y, ylo, yhi, yblo, ncomp, dp) bind(c, name='fort_fab_dot')

    use amrex_fort_module, only: amrex_real
    use basefab_nd_module, only: fort_fab_dot_doit
#ifdef CUDA
    use cuda_module, only: numBlocks, numThreads, cuda_stream
#endif

    implicit none

    integer, intent(in) :: lo(3), hi(3), xlo(3), xhi(3), ylo(3), yhi(3), yblo(3), ncomp
    real(amrex_real), intent(in) :: x(xlo(1):xhi(1),xlo(2):xhi(2),xlo(3):xhi(3),ncomp)
    real(amrex_real), intent(in) :: y(ylo(1):yhi(1),ylo(2):yhi(2),ylo(3):yhi(3),ncomp)
    real(amrex_real), intent(inout) :: dp

#ifdef CUDA
    attributes(device) :: lo, hi, x, xlo, xhi, y, ylo, yhi, yblo, ncomp, dp
#endif

    call fort_fab_dot_doit &
#ifdef CUDA
         <<<numBlocks, numThreads, 0, cuda_stream>>> &
#endif
         (lo, hi, x, xlo, xhi, y, ylo, yhi, yblo, ncomp, dp)

  end subroutine fort_fab_dot
