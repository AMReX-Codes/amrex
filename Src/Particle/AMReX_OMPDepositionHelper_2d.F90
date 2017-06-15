subroutine amrex_atomic_accumulate_fab(local_fab, tile_lo, tile_hi, &
                                       global_fab, lo, hi) &
     bind(c, name='amrex_atomic_accumulate_fab')

  use iso_c_binding
  use amrex_fort_module, only : amrex_real
  
  implicit none
  
  integer              :: lo(2)
  integer              :: hi(2)
  real(amrex_real)     :: global_fab(lo(1):hi(1), lo(2):hi(2))
  integer              :: tile_lo(2)
  integer              :: tile_hi(2)
  real(amrex_real)     :: local_fab(tile_lo(1):tile_hi(1), tile_lo(2):tile_hi(2))

  integer          :: i,j

  do j = tile_lo(2), tile_hi(2)
     do i = tile_lo(1), tile_hi(1)
        !$omp atomic
        global_fab(i,j) = global_fab(i,j) + local_fab(i,j)
        !$omp end atomic
     end do
  end do

end subroutine amrex_atomic_accumulate_fab
