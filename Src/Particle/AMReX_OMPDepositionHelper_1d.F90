subroutine amrex_atomic_accumulate_fab(local_fab, tile_lo, tile_hi, &
                                       global_fab, lo, hi) &
     bind(c, name='amrex_atomic_accumulate_fab')

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  implicit none
  
  integer              :: lo(1)
  integer              :: hi(1)
  real(amrex_real)     :: global_fab(lo(1):hi(1))
  integer              :: tile_lo(1)
  integer              :: tile_hi(1)
  real(amrex_real)     :: local_fab(tile_lo(1):tile_hi(1))

  integer          :: i

  do i = tile_lo(1), tile_hi(1)
     !$omp atomic
     global_fab(i) = global_fab(i) + local_fab(i)
     !$omp end atomic
  end do

end subroutine amrex_atomic_accumulate_fab
