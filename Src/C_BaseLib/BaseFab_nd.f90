module basefab_nd_module

  use bl_fort_module, only : c_real

  implicit none

contains

  ! dst = src
  subroutine fort_fab_copy(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='fort_fab_copy')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(c_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(c_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = lo - sblo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_copy
    
  
  ! copy from multi-d array to 1d array
  subroutine fort_fab_copytomem (lo, hi, dst, src, slo, shi, ncomp) &
       bind(c,name='fort_fab_copytomem')
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp
    real(c_real)             :: dst(*)
    real(c_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)

    integer :: i, j, k, n, nx, offset

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
  end subroutine fort_fab_copytomem


  ! copy from 1d array to multi-d array
  subroutine fort_fab_copyfrommem (lo, hi, dst, dlo, dhi, ncomp, src) &
       bind(c,name='fort_fab_copyfrommem')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(c_real), intent(in   ) :: src(*)
    real(c_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n, nx, offset

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
  end subroutine fort_fab_copyfrommem


  subroutine fort_fab_setval(lo, hi, dst, dlo, dhi, ncomp, val) &
       bind(c,name='fort_fab_setval')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), ncomp
    real(c_real), intent(in) :: val
    real(c_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)

    integer :: i, j, k, n

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = val
             end do
          end do
       end do
    end do
  end subroutine fort_fab_setval


  function fort_fab_norm (lo, hi, src, slo, shi, ncomp, p) result(nrm) &
       bind(c,name='fort_fab_norm')
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp, p
    real(c_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(c_real) :: nrm

    integer :: i,j,k,n

    nrm = 0.0_c_real
    if (p .eq. 0) then ! max norm
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
  end function fort_fab_norm


  function fort_fab_sum (lo, hi, src, slo, shi, ncomp, p) result(sm) &
       bind(c,name='fort_fab_sum')
    integer, intent(in) :: lo(3), hi(3), slo(3), shi(3), ncomp, p
    real(c_real), intent(in) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(c_real) :: sm

    integer :: i,j,k,n

    sm = 0.0_c_real
    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                sm = sm + src(i,j,k,n)
             end do
          end do
       end do
    end do
  end function fort_fab_sum


  subroutine fort_fab_plus(lo, hi, dst, dlo, dhi, src, slo, shi, sblo, ncomp) &
       bind(c,name='fort_fab_plus')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(c_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(c_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = lo - sblo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_plus


  ! dst += a*src
  subroutine fort_fab_saxpy(lo, hi, dst, dlo, dhi, a, src, slo, shi, sblo, ncomp) &
       bind(c,name='fort_fab_saxpy')
    integer, intent(in) :: lo(3), hi(3), dlo(3), dhi(3), slo(3), shi(3), sblo(3), ncomp
    real(c_real), intent(in   ) :: a
    real(c_real), intent(in   ) :: src(slo(1):shi(1),slo(2):shi(2),slo(3):shi(3),ncomp)
    real(c_real), intent(inout) :: dst(dlo(1):dhi(1),dlo(2):dhi(2),dlo(3):dhi(3),ncomp)
    
    integer :: i,j,k,n,off(3)

    off = lo - sblo

    do n = 1, ncomp
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,n) = dst(i,j,k,n) + a * src(i+off(1),j+off(2),k+off(3),n)
             end do
          end do
       end do
    end do
  end subroutine fort_fab_saxpy

end module basefab_nd_module
