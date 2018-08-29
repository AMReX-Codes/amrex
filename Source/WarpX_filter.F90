module warpx_filter_module

  use iso_c_binding
  use amrex_fort_module, only : rt => amrex_real
  use amrex_constants_module, only : fourth, eighth

  implicit none

contains

  subroutine warpx_filter_3d(lo, hi, src, slo, shi, dst, dlo, dhi, nc) &
     bind(c, name='warpx_filter_3d')
    integer, intent(in), value :: nc
    integer, dimension(3), intent(in) :: lo, hi, dlo, dhi, slo, shi
    real(rt), intent(inout) :: dst(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3),nc)
    real(rt), intent(in   ) :: src(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3),nc)
    
    real(rt), parameter :: sixteenth = 1.d0/16.d0
    real(rt), parameter :: thirtysecond = 1.d0/32.d0
    real(rt), parameter :: sixtyfourth = 1.d0/64.d0
    integer          :: i,j,k,c

    do c = 1, nc
       do       k = lo(3), hi(3)
          do    j = lo(2), hi(2)
             do i = lo(1), hi(1)
                dst(i,j,k,c) = sixtyfourth*(src(i-1,j-1,k-1,c) &
                     &                     +src(i+1,j-1,k-1,c) &
                     &                     +src(i-1,j+1,k-1,c) &
                     &                     +src(i+1,j+1,k-1,c) &
                     &                     +src(i-1,j-1,k+1,c) &
                     &                     +src(i+1,j-1,k+1,c) &
                     &                     +src(i-1,j+1,k+1,c) &
                     &                     +src(i+1,j+1,k+1,c)) &
                     +        thirtysecond*(src(i-1,j-1,k  ,c) &
                     &                     +src(i+1,j-1,k  ,c) &
                     &                     +src(i-1,j+1,k  ,c) &
                     &                     +src(i+1,j+1,k  ,c) &
                     &                     +src(i-1,j  ,k-1,c) &
                     &                     +src(i+1,j  ,k-1,c) &
                     &                     +src(i-1,j  ,k+1,c) &
                     &                     +src(i+1,j  ,k+1,c) &
                     &                     +src(i  ,j-1,k-1,c) &
                     &                     +src(i  ,j+1,k-1,c) &
                     &                     +src(i  ,j-1,k+1,c) &
                     &                     +src(i  ,j+1,k+1,c)) &
                     +        sixteenth   *(src(i-1,j  ,k  ,c) &
                     &                     +src(i+1,j  ,k  ,c) &
                     &                     +src(i  ,j-1,k  ,c) &
                     &                     +src(i  ,j+1,k  ,c) &
                     &                     +src(i  ,j  ,k-1,c) &
                     &                     +src(i  ,j  ,k+1,c)) &
                     +        eighth      * src(i  ,j  ,k  ,c)
             end do
          end do
       end do
    end do
  end subroutine warpx_filter_3d

  subroutine warpx_filter_2d(lo, hi, src, slo, shi, dst, dlo, dhi, nc) &
     bind(c, name='warpx_filter_2d')
    integer, intent(in), value :: nc
    integer, dimension(2), intent(in) :: lo, hi, dlo, dhi, slo, shi
    real(rt), intent(inout) :: dst(dlo(1):dhi(1), dlo(2):dhi(2), nc)
    real(rt), intent(in   ) :: src(slo(1):shi(1), slo(2):shi(2), nc)
    
    real(rt), parameter :: sixteenth = 1.d0/16.d0
    integer          :: i,j,comp

    do comp = 1, nc
       do    j = lo(2), hi(2)
          do i = lo(1), hi(1)
             dst(i, j, comp) = sixteenth*(src(i-1,j-1,comp) &
                  &                      +src(i+1,j-1,comp) &
                  &                      +src(i-1,j+1,comp) &
                  &                      +src(i+1,j+1,comp)) &
                  +            eighth   *(src(i  ,j-1,comp) &
                  &                      +src(i  ,j+1,comp) &
                  &                      +src(i-1,j  ,comp) &
                  &                      +src(i+1,j  ,comp)) &
                  +            fourth   * src(i  ,j  ,comp)
          end do
       end do
    end do
  end subroutine warpx_filter_2d

end module warpx_filter_module
