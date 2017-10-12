module warpx_filter_module

  use iso_c_binding
  use amrex_fort_module, only : rt => amrex_real

  implicit none

contains

  subroutine warpx_filter_3d(src_fab, slo, shi, &
                             dst_fab, dlo, dhi, nc) &
     bind(c, name='warpx_filter_3d')

    integer, value       :: nc
    integer              :: dlo(3), dhi(3), slo(3), shi(3)
    real(rt)             :: dst_fab(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3),nc)
    real(rt)             :: src_fab(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), nc)
    
    integer          :: i,j,k,comp

    do comp = 1, nc

       ! pass in the x-direction
       do k = dlo(3), dhi(3)
          do j = dlo(2), dhi(2)
             do i = dlo(1), dhi(1)
                dst_fab(i, j, k, comp) = 0.25d0 * src_fab(i-1, j, k, comp) + & 
                                         0.5d0  * src_fab(i  , j, k, comp) + &
                                         0.25d0 * src_fab(i+1, j, k, comp)
             end do
          end do
       end do

       ! pass in the y-direction
       do k = dlo(3), dhi(3)
          do j = dlo(2), dhi(2)
             do i = dlo(1), dhi(1)
                dst_fab(i, j, k, comp) = 0.25d0 * src_fab(i, j-1, k, comp) + & 
                                         0.5d0  * src_fab(i, j  , k, comp) + &
                                         0.25d0 * src_fab(i, j+1, k, comp)
             end do
          end do
       end do

       ! pass in the z-direction
       do k = dlo(3), dhi(3)
          do j = dlo(2), dhi(2)
             do i = dlo(1), dhi(1)
                dst_fab(i, j, k, comp) = 0.25d0 * src_fab(i, j, k-1, comp) + & 
                                         0.5d0  * src_fab(i, j, k  , comp) + &
                                         0.25d0 * src_fab(i, j, k+1, comp)
             end do
          end do
       end do

    end do

  end subroutine warpx_filter_3d

  subroutine warpx_filter_2d(src_fab, slo, shi, &
                             dst_fab, dlo, dhi, nc) &
     bind(c, name='warpx_filter_2d')

    integer, value       :: nc
    integer              :: dlo(2), dhi(2), slo(2), shi(2)
    real(rt)             :: dst_fab(dlo(1):dhi(1), dlo(2):dhi(2), nc)
    real(rt)             :: src_fab(slo(1):shi(1), slo(2):shi(2), nc)
    
    integer          :: i,j,comp

    do comp = 1, nc

       ! pass in the x-direction
       do j = dlo(2), dhi(2)
          do i = dlo(1), dhi(1)
             dst_fab(i, j, comp) = 0.25d0 * src_fab(i-1, j, comp) + & 
                                   0.5d0  * src_fab(i  , j, comp) + &
                                   0.25d0 * src_fab(i+1, j, comp)
          end do
       end do

       ! pass in the z-direction
       do j = dlo(2), dhi(2)
          do i = dlo(1), dhi(1)
             dst_fab(i, j, comp) = 0.25d0 * src_fab(i, j-1, comp) + & 
                                   0.5d0  * src_fab(i, j  , comp) + &
                                   0.25d0 * src_fab(i, j+1, comp)
          end do
       end do

    end do

  end subroutine warpx_filter_2d

end module warpx_filter_module
