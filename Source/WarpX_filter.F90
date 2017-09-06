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

  subroutine warpx_filter_and_accumulate_3d(local_fab, tile_lo, tile_hi, &
                                            global_fab, lo, hi, nc) &
     bind(c, name='warpx_filter_and_accumulate_3d')

    integer, value       :: nc
    integer              :: lo(3)
    integer              :: hi(3)
    real(rt)             :: global_fab(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
    integer              :: tile_lo(3)
    integer              :: tile_hi(3)
    real(rt)             :: local_fab(tile_lo(1):tile_hi(1), tile_lo(2):tile_hi(2), tile_lo(3):tile_hi(3), nc)
    
    integer          :: i,j,k,comp
    integer          :: ii, jj, kk
    
    do comp = 1, nc
       do k = tile_lo(3), tile_hi(3)
          do j = tile_lo(2), tile_hi(2)
             do i = tile_lo(1), tile_hi(1)

                do kk = -1, 1
                   do jj = -1, 1
                      do ii = -1, 1
#ifdef _OPENMP
                         !$omp atomic
#endif

                         global_fab(i+ii,j+jj,k+kk,comp) = global_fab(i+ii,j+jj,k+kk,comp) + ff1_3D(ii, jj, kk) * local_fab(i,j,k,comp)
#ifdef _OPENMP
                         !$omp end atomic
#endif
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
