subroutine amrex_atomic_accumulate_fab(local_fab, tile_lo, tile_hi, &
                                       global_fab, lo, hi, nc) &
     bind(c, name='amrex_atomic_accumulate_fab')

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  implicit none
  
  integer, value       :: nc
  integer              :: lo(3)
  integer              :: hi(3)
  real(amrex_real)     :: global_fab(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3),nc)
  integer              :: tile_lo(3)
  integer              :: tile_hi(3)
  real(amrex_real)     :: local_fab(tile_lo(1):tile_hi(1), tile_lo(2):tile_hi(2), tile_lo(3):tile_hi(3), nc)

  integer          :: i,j,k,comp

  do comp = 1, nc
     do k = tile_lo(3), tile_hi(3)
        do j = tile_lo(2), tile_hi(2)
           do i = tile_lo(1), tile_hi(1)
#ifdef _OPENMP
              !$omp atomic
#endif
              global_fab(i,j,k,comp) = global_fab(i,j,k,comp) + local_fab(i,j,k,comp)
#ifdef _OPENMP
              !$omp end atomic
#endif
           end do
        end do
     end do
  end do

end subroutine amrex_atomic_accumulate_fab
