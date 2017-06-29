subroutine amrex_atomic_accumulate_fab(local_fab, tile_lo, tile_hi, &
                                       global_fab, lo, hi, nc) &
     bind(c, name='amrex_atomic_accumulate_fab')

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  implicit none

  integer, value       :: nc
  integer              :: lo(1)
  integer              :: hi(1)
  real(amrex_real)     :: global_fab(lo(1):hi(1),nc)
  integer              :: tile_lo(1)
  integer              :: tile_hi(1)
  real(amrex_real)     :: local_fab(tile_lo(1):tile_hi(1),nc)

  integer          :: i, comp

  do comp = 1, nc
     do i = tile_lo(1), tile_hi(1)
        !$omp atomic
        global_fab(i, comp) = global_fab(i, comp) + local_fab(i, comp)
        !$omp end atomic
     end do
  end do

end subroutine amrex_atomic_accumulate_fab
