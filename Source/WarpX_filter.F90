module warpx_filter_module

  use iso_c_binding
  use amrex_fort_module, only : rt => amrex_real

  implicit none

  ! 1D, 3-point filter:
  real(rt), dimension(-1:1), parameter :: ff1_1D = &
       (/ 0.25e0_rt, 0.5e0_rt, 0.25e0_rt/)

  ! 2D, 3-point filter:
  real(rt), dimension(-1:1, -1:1), parameter :: ff1_2D = reshape( &
       [ 0.0625e0_rt, 0.125e0_rt, 0.0625e0_rt, &
         0.125e0_rt,  0.25e0_rt,  0.125e0_rt,  &
         0.0625e0_rt, 0.125e0_rt, 0.0625e0_rt], &
       shape(ff1_2D) )

  ! 3D, 3-point filter:
  real(rt), dimension(-1:1, -1:1, -1:1), parameter :: ff1_3D = reshape( &
       [ 0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt,  &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt,  &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.0625e0_rt,   0.125e0_rt,    0.0625e0_rt,    &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt,  &
         0.03125e0_rt,  0.0625e0_rt,   0.03125e0_rt,   &
         0.015625e0_rt, 0.03125e0_rt,  0.015625e0_rt], &
       shape(ff1_3D) )

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

  subroutine warpx_filter_and_accumulate_3d(src_fab, slo, shi, &
                                            dst_fab, dlo, dhi, nc) &
     bind(c, name='warpx_filter_and_accumulate_3d')

    integer, value       :: nc
    integer              :: dlo(3), dhi(3), slo(3), shi(3)
    real(rt)             :: dst_fab(dlo(1):dhi(1), dlo(2):dhi(2), dlo(3):dhi(3),nc)
    real(rt)             :: src_fab(slo(1):shi(1), slo(2):shi(2), slo(3):shi(3), nc)
    
    integer          :: i,j,k,comp
    integer          :: ii, jj, kk
    
    do comp = 1, nc
       do k = dlo(3), dhi(3)
          do j = dlo(2), dhi(2)
             do i = dlo(1), dhi(1)

                do kk = -1, 1
                   do jj = -1, 1
                      do ii = -1, 1

                         dst_fab(i,j,k,comp) = dst_fab(i,j,k,comp) + ff1_3D(ii, jj, kk) * src_fab(i+ii,j+jj,k+kk,comp)

                      end do
                   end do
                end do
                
             end do
          end do
       end do
    end do
    
  end subroutine warpx_filter_and_accumulate_3d

  subroutine warpx_filter_and_accumulate_2d(local_fab, tile_lo, tile_hi, &
                                            global_fab, lo, hi, nc) &
       bind(c, name='warpx_filter_and_accumulate_2d')

    integer, value       :: nc
    integer              :: lo(2)
    integer              :: hi(2)
    real(rt)             :: global_fab(lo(1):hi(1), lo(2):hi(2), nc)
    integer              :: tile_lo(2)
    integer              :: tile_hi(2)
    real(rt)             :: local_fab(tile_lo(1):tile_hi(1), tile_lo(2):tile_hi(2), nc)

    integer              :: i,k,comp
    integer              :: ii, kk

    do comp = 1, nc
       do k = tile_lo(2), tile_hi(2)
          do i = tile_lo(1), tile_hi(1)

                do kk = -1, 1
                   do ii = -1, 1
#ifdef _OPENMP
                      !$omp atomic
#endif
                      
                      global_fab(i+ii,k+kk,comp) = global_fab(i+ii,k+kk,comp) + ff1_2D(ii, kk) * local_fab(i,k,comp)
#ifdef _OPENMP
                      !$omp end atomic
#endif
                   end do
                end do

          end do
       end do
    end do

  end subroutine warpx_filter_and_accumulate_2d

end module warpx_filter_module
